# -*- coding: utf-8 -*-
# @Time    : 2025/4/19


r"""
FreqRepFusion
################################################
Reference:

"""

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
import math

from recbole.model.abstract_recommender import SequentialRecommender
from recbole.model.loss import BPRLoss
from recbole.data.interaction import Interaction


class EnhancedFrequencyDomainDecoupler(nn.Module):
    """将用户兴趣分解为高频(重复兴趣)和低频(次要兴趣)两个独立表示，着重高频表示
    
    Args:
        embedding_dim (int): 嵌入维度
        cutoff_ratio (float): 高低频分界点比例, 默认0.5 
        high_freq_scale (float): 高频增强比例, 默认1.5
    """
    def __init__(self, embedding_dim, cutoff_ratio=0.5, high_freq_scale=1.5):
        super(EnhancedFrequencyDomainDecoupler, self).__init__()
        self.embedding_dim = embedding_dim
        self.cutoff_ratio = cutoff_ratio
        self.high_freq_scale = high_freq_scale
        
        # 特征提取器 
        self.encoder = nn.Sequential(
            nn.Linear(embedding_dim, embedding_dim*2),
            nn.LayerNorm(embedding_dim*2),  
            nn.GELU(),  
            nn.Dropout(0.1),  
            nn.Linear(embedding_dim*2, embedding_dim)
        )
        
        # 高频投影层
        self.high_freq_projector = nn.Sequential(
            nn.Linear(embedding_dim, embedding_dim*2),
            nn.LayerNorm(embedding_dim*2),
            nn.ReLU(),
            nn.Linear(embedding_dim*2, embedding_dim),
            nn.LayerNorm(embedding_dim),
        )
        
        # 低频投影层
        self.low_freq_projector = nn.Sequential(
            nn.Linear(embedding_dim, embedding_dim),
            nn.LayerNorm(embedding_dim),
            nn.ReLU()
        )
        
    def forward(self, seq_output):
        """增强高频信息的提取和表示
        
        Args:
            seq_output (torch.Tensor): 序列表示, shape: [B, H]
            
        Returns:
            tuple: 高频表示, 低频表示, 频域表示
        """
        if torch.isnan(seq_output).any():
            print("输入seq_output中检测到NaN！")
            seq_output = torch.nan_to_num(seq_output, nan=0.0, posinf=1.0, neginf=-1.0)

        seq_output = seq_output.unsqueeze(1)  # [B, 1, H]
        
        # 编码用户行为序列
        encoded = self.encoder(seq_output)  # [B, 1, H]

        if torch.isnan(encoded).any():
            print("encoded中检测到NaN！")
            encoded = torch.nan_to_num(encoded, nan=0.0, posinf=1.0, neginf=-1.0)

        encoded_mean = encoded.mean(dim=2, keepdim=True)
        encoded_std = encoded.std(dim=2, keepdim=True) + 1e-5
        encoded_norm = (encoded - encoded_mean) / encoded_std
        
        # 转换到频域 (对嵌入维度进行FFT)
        freq_domain = torch.fft.rfft(encoded_norm, dim=2, norm="ortho")  # [B, 1, H//2+1]
        
        # 高低频划分
        cutoff = int(freq_domain.shape[2] * self.cutoff_ratio)
        high_freq_mask = torch.zeros_like(freq_domain)
        low_freq_mask = torch.zeros_like(freq_domain)
        
        # 高频成分 (重复兴趣)
        high_freq_mask[:, :, :cutoff] = 1.0 
        
        # 低频成分 (次要兴趣)
        low_freq_mask[:, :, cutoff:] = 1.0
        
        # 高频增强处理
        high_freq_components = freq_domain * high_freq_mask * self.high_freq_scale
        low_freq_components = freq_domain * low_freq_mask
        
        # 逆变换回时域
        high_freq_repr = torch.fft.irfft(high_freq_components, n=self.embedding_dim, dim=2, norm="ortho")  # [B, 1, H]
        low_freq_repr = torch.fft.irfft(low_freq_components, n=self.embedding_dim, dim=2, norm="ortho")  # [B, 1, H]

        if torch.isnan(high_freq_repr).any() or torch.isnan(low_freq_repr).any():
            print("逆FFT后检测到NaN！")
        
        # 维度压缩
        high_freq_repr = high_freq_repr.squeeze(1)  # [B, H]
        low_freq_repr = low_freq_repr.squeeze(1)  # [B, H]
        
        # 投影到表示空间
        high_freq_repr = self.high_freq_projector(high_freq_repr)  # [B, H]
        low_freq_repr = self.low_freq_projector(low_freq_repr)  # [B, H]

        if torch.isnan(high_freq_repr).any() or torch.isnan(low_freq_repr).any():
            print("投影后检测到NaN！")
            high_freq_repr = torch.nan_to_num(high_freq_repr, nan=0.0, posinf=1.0, neginf=-1.0)
            low_freq_repr = torch.nan_to_num(low_freq_repr, nan=0.0, posinf=1.0, neginf=-1.0)
        
        return high_freq_repr, low_freq_repr, freq_domain


class RepetitiveInterestTracker(nn.Module):
    """跟踪用户重复兴趣模式及其强度变化
    
    Args:
        embedding_dim (int): 嵌入维度
        window_size (int): 频谱历史窗口大小, 默认5
        rep_threshold (float): 重复兴趣阈值, 默认0.6
    """
    def __init__(self, embedding_dim, window_size=5, rep_threshold=0.6):
        super(RepetitiveInterestTracker, self).__init__()
        self.embedding_dim = embedding_dim
        self.window_size = window_size
        self.rep_threshold = rep_threshold
        
        # 重复模式检测网络 
        self.pattern_detector = nn.Sequential(
                nn.Linear(embedding_dim//4, embedding_dim//2),
                nn.LayerNorm(embedding_dim//2, eps=1e-5),
                nn.GELU(),
                nn.Dropout(0.05),  
                nn.Linear(embedding_dim//2, embedding_dim//4),
                nn.LayerNorm(embedding_dim//4, eps=1e-5),
                nn.GELU()
        )
        
        # 重复强度分类器
        self.repetition_classifier = nn.Sequential(
            nn.Linear(embedding_dim//4, embedding_dim),
            nn.LayerNorm(embedding_dim),
            nn.GELU(),
            nn.Linear(embedding_dim, embedding_dim//2),
            nn.GELU(),
            nn.Linear(embedding_dim//2, 1),
            nn.Sigmoid()
        )
        
    def forward(self, freq_history):
        """分析用户重复兴趣强度
        
        Args:
            freq_history (torch.Tensor): 频谱历史序列, shape: [B, T, F]
                T: 历史长度, F: 频域维度 (embedding_dim//2+1)
            
        Returns:
            tuple: 重复强度分数, 重复兴趣模式
        """
        batch_size = freq_history.shape[0]
       
        if torch.isnan(freq_history).any():
            print("freq_history中检测到NaN！")
            freq_history = torch.nan_to_num(freq_history, nan=0.0, posinf=1.0, neginf=-1.0)
        
        # 历史处理
        if freq_history.shape[1] < self.window_size:
            padding = torch.zeros(
                (batch_size, self.window_size - freq_history.shape[1], freq_history.shape[2]),
                dtype=freq_history.dtype,
                device=freq_history.device
            )
            freq_history = torch.cat([padding, freq_history], dim=1)
        
        # 提取高频部分幅度 (只关注高频部分)
        cutoff = int(freq_history.shape[2] * 0.5)
        high_freq_magnitude = torch.abs(freq_history[:, :, :cutoff])
        
        high_freq_magnitude = torch.clamp(high_freq_magnitude, min=1e-5, max=10.0)
        mean = high_freq_magnitude.mean(dim=(1,2), keepdim=True)
        std = high_freq_magnitude.std(dim=(1,2), keepdim=True) + 1e-5
        high_freq_magnitude = (high_freq_magnitude - mean) / std
        high_freq_magnitude = torch.clamp(high_freq_magnitude, min=-3, max=3)
        
        latest_features = high_freq_magnitude[:, -1, :]
        pattern_features = self.pattern_detector(latest_features)

        if torch.isnan(pattern_features).any():
            print("pattern_features中检测到NaN！")
            pattern_features = torch.nan_to_num(pattern_features, nan=0.0, posinf=1.0, neginf=-1.0)
        
        
        # 计算重复强度分数 (0-1)
        repetition_scores = self.repetition_classifier(pattern_features)

        if torch.isnan(repetition_scores).any():
            print("repetition_scores中检测到NaN！")
            repetition_scores = torch.nan_to_num(repetition_scores, nan=0.5, posinf=1.0, neginf=0.0)
        
        # 分析重复模式变化
        repetition_patterns = None
        if freq_history.shape[1] >= 2:
            try:
                # 计算高频幅度变化率 
                first_window = torch.abs(freq_history[:, 0, :cutoff])
                last_window = torch.abs(freq_history[:, -1, :cutoff])
                first_window = torch.clamp(first_window, min=1e-5)  # 避免除零
                amplitude_change = (last_window - first_window) / (first_window + 1e-5)
            
                # 提取增长和衰减模式
                growth_mask = amplitude_change > self.rep_threshold
                decay_mask = amplitude_change < -self.rep_threshold
            
                repetition_patterns = {
                    'growth_patterns': growth_mask,  # 增强的重复兴趣
                    'decay_patterns': decay_mask,    # 衰减的重复兴趣
                    'repetition_scores': repetition_scores
                }
            except Exception as e:
                print(f"计算重复模式时出错: {e}")
                repetition_patterns = {
                    'repetition_scores': repetition_scores
                }
        
        return repetition_scores, repetition_patterns


class RepetitionAwareRecommender(nn.Module):
    """基于重复兴趣的融合推荐器
    
    Args:
        embedding_dim (int): 嵌入维度
        rep_weight (float): 重复兴趣权重基准, 默认0.7
    """
    def __init__(self, embedding_dim, rep_weight=0.7):
        super(RepetitionAwareRecommender, self).__init__()
        self.embedding_dim = embedding_dim
        self.rep_weight = rep_weight
        
        # 融合网络 - 动态调整重复兴趣权重
        self.fusion_network = nn.Sequential(
            nn.Linear(embedding_dim * 2 + 1, embedding_dim),
            nn.LayerNorm(embedding_dim),
            nn.GELU(),
            nn.Dropout(0.1),
            nn.Linear(embedding_dim, embedding_dim//2),
            nn.GELU(),
            nn.Linear(embedding_dim//2, 1),
            nn.Sigmoid()
        )
        
        # 重复兴趣注意力
        self.repetition_attention = nn.MultiheadAttention(
            embed_dim=embedding_dim,
            num_heads=4,
            dropout=0.1,
            batch_first=True
        )
        
        # 最终融合投影
        self.final_projection = nn.Sequential(
            nn.Linear(embedding_dim * 2, embedding_dim),
            nn.LayerNorm(embedding_dim),
            nn.GELU()
        )
        
    def forward(self, high_freq_repr, low_freq_repr, repetition_scores):
        """融合高频和低频表示
        
        Args:
            high_freq_repr (torch.Tensor): 高频表示 (重复兴趣), shape: [B, H]
            low_freq_repr (torch.Tensor): 低频表示 (次要兴趣), shape: [B, H]
            repetition_scores (torch.Tensor): 重复强度分数, shape: [B, 1]
            
        Returns:
            tuple: 融合表示, 高频权重, 低频权重
        """
        batch_size = high_freq_repr.shape[0]
        
        # 动态计算重复兴趣权重
        fusion_input = torch.cat([
            high_freq_repr, 
            low_freq_repr, 
            repetition_scores
        ], dim=1)
        
        # 调整基准权重
        dynamic_weight = self.fusion_network(fusion_input)
        repetition_weight = self.rep_weight + (1 - self.rep_weight) * dynamic_weight
        novelty_weight = 1 - repetition_weight
        
        # 高频特征自注意力增强
        high_freq_enhanced = high_freq_repr.unsqueeze(1)
        high_freq_attn, _ = self.repetition_attention(
            high_freq_enhanced, 
            high_freq_enhanced, 
            high_freq_enhanced
        )
        high_freq_attn = high_freq_attn.squeeze(1)
        
        # 基于权重融合
        weighted_high = repetition_weight * high_freq_attn
        weighted_low = novelty_weight * low_freq_repr
        
        # 连接并最终投影
        concat_repr = torch.cat([weighted_high, weighted_low], dim=1)
        final_repr = self.final_projection(concat_repr)
        
        return final_repr, repetition_weight, novelty_weight


class FreqRepFusion(SequentialRecommender):
    """频域重复兴趣融合推荐模型
    
    基于频域分析的推荐模型，将用户兴趣分解为高频(重复兴趣)和低频(次要兴趣)，
    重点关注重复兴趣模式变化，以提高推荐准确性。

    Args:
        config (Config): 全局配置对象
        dataset (Dataset): 数据集
    """
    def __init__(self, config, dataset):
        super(FreqRepFusion, self).__init__(config, dataset)
        
        # 加载参数信息
        self.embedding_size = config['embedding_size']  # 嵌入维度
        self.hidden_size = config['hidden_size']  # 隐层维度
        self.n_layers = config['n_layers']  # 层数
        self.n_heads = config['n_heads']  # 注意力头数
        self.inner_size = config['inner_size']  # 前馈网络内部维度
        self.hidden_dropout_prob = config['hidden_dropout_prob']  # 隐层dropout比例
        self.attn_dropout_prob = config['attn_dropout_prob']  # 注意力dropout比例
        self.hidden_act = config['hidden_act']  # 激活函数
        self.layer_norm_eps = config['layer_norm_eps']  # 层归一化epsilon
        
        # 模型特定参数
        self.cutoff_ratio = config['cutoff_ratio'] if 'cutoff_ratio' in config else 0.5  # 默认提高高频比例
        self.window_size = config['window_size'] if 'window_size' in config else 5  # 频谱历史窗口大小
        self.rep_threshold = config['rep_threshold'] if 'rep_threshold' in config else 0.6  # 重复兴趣阈值
        self.rep_weight = config['rep_weight'] if 'rep_weight' in config else 0.7  # 重复兴趣基准权重
        self.high_freq_scale = config['high_freq_scale'] if 'high_freq_scale' in config else 1.5  # 高频增强系数
        self.contrast_weight = config['contrast_weight'] if 'contrast_weight' in config else 0.5  # 对比学习权重调整
        self.freq_history_length = config['freq_history_length'] if 'freq_history_length' in config else 10  # 频谱历史长度
        
        # 初始化嵌入层
        self.item_embedding = nn.Embedding(
            self.n_items, self.embedding_size, padding_idx=0
        )
        self.position_embedding = nn.Embedding(
            self.max_seq_length, self.embedding_size
        )
        
        # 序列编码器 
        try:
            from recbole.model.sequential_recommender.fearec import FEAEncoder
            self.encoder = FEAEncoder(
                n_layers=self.n_layers,
                n_heads=self.n_heads,
                hidden_size=self.hidden_size,
                inner_size=self.inner_size,
                hidden_dropout_prob=self.hidden_dropout_prob,
                attn_dropout_prob=self.attn_dropout_prob,
                hidden_act=self.hidden_act,
                layer_norm_eps=self.layer_norm_eps,
                config=config
            )
        except ImportError:
            # 如果没有FEARec，使用标准Transformer编码器
            from recbole.model.layers import TransformerEncoder
            self.encoder = TransformerEncoder(
                n_layers=self.n_layers,
                n_heads=self.n_heads,
                hidden_size=self.hidden_size,
                inner_size=self.inner_size,
                hidden_dropout_prob=self.hidden_dropout_prob,
                attn_dropout_prob=self.attn_dropout_prob,
                hidden_act=self.hidden_act,
                layer_norm_eps=self.layer_norm_eps
            )
        
        # 初始化频域相关模块
        self.freq_decoupler = EnhancedFrequencyDomainDecoupler(
            self.embedding_size, 
            self.cutoff_ratio,
            self.high_freq_scale
        )
        self.repetition_tracker = RepetitiveInterestTracker(
            self.embedding_size, 
            self.window_size, 
            self.rep_threshold
        )
        self.recommender = RepetitionAwareRecommender(
            self.embedding_size,
            self.rep_weight
        )
        
        # 层归一化与dropout
        self.LayerNorm = nn.LayerNorm(self.embedding_size, eps=self.layer_norm_eps)
        self.dropout = nn.Dropout(self.hidden_dropout_prob)
        
        # 损失函数
        self.loss_type = config['loss_type']
        if config['loss_type'] == 'BPR':
            self.loss_fct = BPRLoss()
        elif config['loss_type'] == 'CE':
            self.loss_fct = nn.CrossEntropyLoss()
        else:
            raise NotImplementedError("Make sure 'loss_type' in ['BPR', 'CE']!")
        
        # 对比学习损失
        self.contrast_loss = nn.CrossEntropyLoss()
        
        # 频谱历史存储
        self.register_buffer('freq_history', None)
        
        # 参数初始化
        self.apply(self._init_weights)
    
    def _init_weights(self, module):
        """初始化权重"""
        if isinstance(module, (nn.Linear, nn.Embedding)):
            module.weight.data.normal_(mean=0.0, std=0.02)
        elif isinstance(module, nn.LayerNorm):
            module.bias.data.zero_()
            module.weight.data.fill_(1.0)
        if isinstance(module, nn.Linear) and module.bias is not None:
            module.bias.data.zero_()
    
    def get_attention_mask(self, item_seq):
        """生成注意力掩码"""
        attention_mask = (item_seq > 0).long()
        extended_attention_mask = attention_mask.unsqueeze(1).unsqueeze(2)
        max_len = attention_mask.size(-1)
        attn_shape = (1, max_len, max_len)
        subsequent_mask = torch.triu(torch.ones(attn_shape), diagonal=1)
        subsequent_mask = (subsequent_mask == 0).unsqueeze(1)
        subsequent_mask = subsequent_mask.long().to(item_seq.device)
        
        extended_attention_mask = extended_attention_mask * subsequent_mask
        extended_attention_mask = extended_attention_mask.to(
            dtype=next(self.parameters()).dtype
        )
        extended_attention_mask = (1.0 - extended_attention_mask) * -10000.0
        return extended_attention_mask
    
    def forward(self, item_seq, item_seq_len):
        """前向传播，重点关注重复兴趣
        
        Args:
            item_seq (torch.LongTensor): 物品序列, shape: [B, L]
            item_seq_len (torch.LongTensor): 序列长度, shape: [B]
            
        Returns:
            torch.Tensor: 用户表示, shape: [B, H]
        """
        # 序列编码部分
        position_ids = torch.arange(
            item_seq.size(1), dtype=torch.long, device=item_seq.device
        )
        position_ids = position_ids.unsqueeze(0).expand_as(item_seq)
        position_embedding = self.position_embedding(position_ids)
        
        item_emb = self.item_embedding(item_seq)
        input_emb = item_emb + position_embedding
        input_emb = self.LayerNorm(input_emb)
        input_emb = self.dropout(input_emb)
        
        extended_attention_mask = self.get_attention_mask(item_seq)
        
        # 序列编码
        if hasattr(self.encoder, 'output_all_encoded_layers'):
            # FEAEncoder返回所有层的输出
            encoder_outputs = self.encoder(
                input_emb, 
                extended_attention_mask,
                output_all_encoded_layers=True
            )
            seq_output = encoder_outputs[-1]
        else:
            # 标准Transformer编码器直接返回最后层输出
            encoder_output = self.encoder(input_emb, extended_attention_mask)
            if isinstance(encoder_output, list):
                seq_output = encoder_output[-1]
            else:
                 seq_output = encoder_output
        
        output = self.gather_indexes(seq_output, item_seq_len - 1)

        if torch.isnan(output).any():
            print("output中检测到NaN！")
            output = torch.nan_to_num(output, nan=0.0, posinf=1.0, neginf=-1.0)
            
        # 频域解耦 
        high_freq_repr, low_freq_repr, freq_domain = self.freq_decoupler(output)
        
        # 更新频谱历史
        if self.freq_history is None or self.freq_history.shape[0] != freq_domain.shape[0]:
            # 首次运行或批次大小变化时初始化
            self.freq_history = freq_domain.detach()
        else:
            if len(freq_domain.shape) != len(self.freq_history.shape):
                if len(freq_domain.shape) > len(self.freq_history.shape):
                    while len(freq_domain.shape) > len(self.freq_history.shape):
                        freq_domain = freq_domain.squeeze(1)
                else:
                    while len(freq_domain.shape) < len(self.freq_history.shape):
                        freq_domain = freq_domain.unsqueeze(1)
                        
            # 截断历史以防止NaN传播
            if torch.isnan(self.freq_history).any():
                print("freq_history中检测到NaN，重置历史!")
                self.freq_history = freq_domain.detach()
            else:
                self.freq_history = torch.cat([self.freq_history, freq_domain.detach()], dim=1)
                if self.freq_history.shape[1] > self.freq_history_length:
                    self.freq_history = self.freq_history[:, 1:]
        
        # 跟踪重复兴趣 (使用新的重复兴趣跟踪器)
        # 添加安全检查
        try:
            repetition_scores, _ = self.repetition_tracker(self.freq_history)
        except Exception as e:
            print(f"重复兴趣跟踪器出错: {e}")
            # 出错时提供默认值
            repetition_scores = torch.ones((output.shape[0], 1), device=output.device) * 0.5
    
        # 确保repetition_scores没有NaN
        if torch.isnan(repetition_scores).any():
            print("repetition_scores中检测到NaN，使用默认值替代")
            repetition_scores = torch.ones((output.shape[0], 1), device=output.device) * 0.5
    
        # 重复兴趣感知的融合推荐
        fused_repr, _, _ = self.recommender(high_freq_repr, low_freq_repr, repetition_scores)
        
        return fused_repr
    
    def calculate_rep_contrastive_loss(self, high_freq_repr, low_freq_repr, pos_items_emb, neg_items_emb=None):
        """计算重复兴趣对比学习损失 (调整对比学习策略)
        
        Args:
            high_freq_repr (torch.Tensor): 高频表示 (重复兴趣), shape: [B, H]
            low_freq_repr (torch.Tensor): 低频表示 (次要兴趣), shape: [B, H]
            pos_items_emb (torch.Tensor): 正样本嵌入, shape: [B, H]
            neg_items_emb (torch.Tensor, optional): 负样本嵌入, shape: [B, H]
            
        Returns:
            torch.Tensor: 对比学习损失
        """
        batch_size = high_freq_repr.shape[0]
        
        # 重点增强高频表示与正样本的相似度
        high_freq_pos_sim = torch.sum(high_freq_repr * pos_items_emb, dim=1)
        high_freq_target = torch.ones_like(high_freq_pos_sim)
        high_freq_loss = F.mse_loss(torch.sigmoid(high_freq_pos_sim), high_freq_target)
        
        # 如果有负样本，增强高频表示对正负样本的区分
        if neg_items_emb is not None:
            high_freq_neg_sim = torch.sum(high_freq_repr * neg_items_emb, dim=1)
            high_freq_logits = torch.stack([high_freq_pos_sim, high_freq_neg_sim], dim=1)
            high_freq_labels = torch.zeros(batch_size, dtype=torch.long, device=high_freq_repr.device)
            high_freq_contrast = self.contrast_loss(high_freq_logits, high_freq_labels)
            
            # 低频表示对比 (较弱的对比信号)
            low_freq_pos_sim = torch.sum(low_freq_repr * pos_items_emb, dim=1)
            low_freq_neg_sim = torch.sum(low_freq_repr * neg_items_emb, dim=1)
            low_freq_logits = torch.stack([low_freq_pos_sim, low_freq_neg_sim], dim=1)
            low_freq_contrast = self.contrast_loss(low_freq_logits, high_freq_labels)
            
            # 合并损失 (高频损失权重更大)
            contrastive_loss = high_freq_loss + 0.8 * high_freq_contrast + 0.2 * low_freq_contrast
        else:
            # 无负样本时
            contrastive_loss = high_freq_loss
        
        return contrastive_loss
    
    def calculate_loss(self, interaction):
        """计算模型损失 (调整损失计算)
        
        Args:
            interaction (Interaction): 交互数据
            
        Returns:
            torch.Tensor: 损失值
        """
        item_seq = interaction[self.ITEM_SEQ]
        item_seq_len = interaction[self.ITEM_SEQ_LEN]
        pos_items = interaction[self.POS_ITEM_ID]
        
        seq_output = self.forward(item_seq, item_seq_len)
        
        # 频域解耦表示 --重新计算是为了训练频域模块
        high_freq_repr, low_freq_repr, _ = self.freq_decoupler(seq_output)
        
        # 计算主要推荐损失
        if self.loss_type == 'BPR':
            neg_items = interaction[self.NEG_ITEM_ID]
            pos_items_emb = self.item_embedding(pos_items)
            neg_items_emb = self.item_embedding(neg_items)
            pos_score = torch.sum(seq_output * pos_items_emb, dim=-1)
            neg_score = torch.sum(seq_output * neg_items_emb, dim=-1)
            main_loss = self.loss_fct(pos_score, neg_score)
            
            # 重复兴趣对比学习损失
            contrastive_loss = self.calculate_rep_contrastive_loss(
                high_freq_repr, low_freq_repr, pos_items_emb, neg_items_emb
            )
        else:  # CE损失
            test_item_emb = self.item_embedding.weight
            logits = torch.matmul(seq_output, test_item_emb.transpose(0, 1))
            main_loss = self.loss_fct(logits, pos_items)
            
            # 重复兴趣对比学习损失
            pos_items_emb = self.item_embedding(pos_items)
            contrastive_loss = self.calculate_rep_contrastive_loss(
                high_freq_repr, low_freq_repr, pos_items_emb
            )
        
        # 组合损失 
        loss = main_loss + self.contrast_weight * contrastive_loss
        
        return loss
    
    def predict(self, interaction):
        """预测单个物品的分数
        
        Args:
            interaction (Interaction): 交互数据
            
        Returns:
            torch.Tensor: 预测分数
        """
        item_seq = interaction[self.ITEM_SEQ]
        item_seq_len = interaction[self.ITEM_SEQ_LEN]
        test_item = interaction[self.ITEM_ID]
        
        seq_output = self.forward(item_seq, item_seq_len)
        test_item_emb = self.item_embedding(test_item)
        scores = torch.mul(seq_output, test_item_emb).sum(dim=1)
        
        return scores
    
    def full_sort_predict(self, interaction):
        """预测所有物品的分数
        
        Args:
            interaction (Interaction): 交互数据
            
        Returns:
            torch.Tensor: 所有物品的预测分数
        """
        item_seq = interaction[self.ITEM_SEQ]
        item_seq_len = interaction[self.ITEM_SEQ_LEN]
        
        seq_output = self.forward(item_seq, item_seq_len)
        test_items_emb = self.item_embedding.weight
        scores = torch.matmul(seq_output, test_items_emb.transpose(0, 1))
        
        return scores
