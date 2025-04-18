# -*- coding: utf-8 -*-
# @Time    : 2025/4/10


r"""
FreqDualDynamic
################################################
Reference:
    OuyangNingtai et al. "Frequency Domain Dual-Channel Dynamic Recommendation for Reducing Repetition."
    In RecSys 2025.
"""

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
import math

from recbole.model.abstract_recommender import SequentialRecommender
from recbole.model.loss import BPRLoss
from recbole.data.interaction import Interaction


class FrequencyDomainDecoupler(nn.Module):
    """将用户兴趣分解为高频(常见兴趣)和低频(独特兴趣)两个独立表示
    
    Args:
        embedding_dim (int): 嵌入维度
        cutoff_ratio (float): 高低频分界点比例, 默认0.3
    """
    def __init__(self, embedding_dim, cutoff_ratio=0.3):
        super(FrequencyDomainDecoupler, self).__init__()
        self.embedding_dim = embedding_dim
        self.cutoff_ratio = cutoff_ratio
        
        # 特征提取器
        self.encoder = nn.Sequential(
            nn.Linear(embedding_dim, embedding_dim*2),
            nn.LayerNorm(embedding_dim*2),  
            nn.ReLU(),
            nn.Dropout(0.2), 
            nn.Linear(embedding_dim*2, embedding_dim)
        )
        
        # 频域投影层
        self.high_freq_projector = nn.Sequential(
            nn.Linear(embedding_dim, embedding_dim),
            nn.LayerNorm(embedding_dim),
            nn.ReLU()
        )
        self.low_freq_projector = nn.Sequential(
            nn.Linear(embedding_dim, embedding_dim),
            nn.LayerNorm(embedding_dim),
            nn.ReLU()
        )
        
    def forward(self, seq_output):
        """
        Args:
            seq_output (torch.Tensor): 序列表示, shape: [B, H]
            
        Returns:
            tuple: 高频表示, 低频表示, 频域表示
        """
        # 扩展维度以适应FFT (FFT需要在序列维度上操作)
        seq_output = seq_output.unsqueeze(1)  # [B, 1, H]
        
        # 编码用户行为序列
        encoded = self.encoder(seq_output)  # [B, 1, H]
        
        # 转换到频域 (对嵌入维度进行FFT)
        freq_domain = torch.fft.rfft(encoded, dim=2, norm="ortho")  # [B, 1, H//2+1]
        
        # 划分高频和低频区域
        cutoff = int(freq_domain.shape[2] * self.cutoff_ratio)
        high_freq_mask = torch.zeros_like(freq_domain)
        low_freq_mask = torch.zeros_like(freq_domain)
        
        # 高频成分 (一般兴趣/重复内容倾向)
        high_freq_mask[:, :, :cutoff] = 1.0 
        
        # 低频成分 (独特偏好/新颖内容倾向)
        low_freq_mask[:, :, cutoff:] = 1.0
        
        # 分离高低频表示
        high_freq_components = freq_domain * high_freq_mask
        low_freq_components = freq_domain * low_freq_mask
        
        # 逆变换回时域
        high_freq_repr = torch.fft.irfft(high_freq_components, n=self.embedding_dim, dim=2, norm="ortho")  # [B, 1, H]
        low_freq_repr = torch.fft.irfft(low_freq_components, n=self.embedding_dim, dim=2, norm="ortho")  # [B, 1, H]
        
        # 维度压缩
        high_freq_repr = high_freq_repr.squeeze(1)  # [B, H]
        low_freq_repr = low_freq_repr.squeeze(1)  # [B, H]
        
        # 投影到表示空间
        high_freq_repr = self.high_freq_projector(high_freq_repr)  # [B, H]
        low_freq_repr = self.low_freq_projector(low_freq_repr)  # [B, H]
        
        return high_freq_repr, low_freq_repr, freq_domain


class FrequencyEvolutionTracker(nn.Module):
    """跟踪用户兴趣频谱演化，检测用户对特定内容类型的疲劳
    
    Args:
        embedding_dim (int): 嵌入维度
        window_size (int): 频谱历史窗口大小
        decay_threshold (float): 频谱衰减阈值, 默认0.3
    """
    def __init__(self, embedding_dim, window_size=5, decay_threshold=0.3):
        super(FrequencyEvolutionTracker, self).__init__()
        self.embedding_dim = embedding_dim
        self.window_size = window_size
        self.decay_threshold = decay_threshold
        
        # 频谱变化检测网络
        self.change_detector = nn.GRU(
            input_size=embedding_dim//2+1,  # FFT结果的维度
            hidden_size=embedding_dim,
            num_layers=2,
            batch_first=True
        )
        
        # 疲劳分类器
        self.fatigue_classifier = nn.Sequential(
            nn.Linear(embedding_dim, 128),
            nn.ReLU(),
            nn.Linear(128, 1),
            nn.Sigmoid()
        )
        
    def forward(self, freq_history):
        """
        Args:
            freq_history (torch.Tensor): 频谱历史序列, shape: [B, T, F]
                T: 历史长度, F: 频域维度 (embedding_dim//2+1)
            
        Returns:
            tuple: 疲劳分数, 疲劳内容类型
        """
        batch_size = freq_history.shape[0]
        
        # 频谱历史序列处理
        # 对于新用户或历史不足的情况，填充历史
        if freq_history.shape[1] < self.window_size:
            padding = torch.zeros(
                (batch_size, self.window_size - freq_history.shape[1], freq_history.shape[2]),
                dtype=freq_history.dtype,
                device=freq_history.device
            )
            freq_history = torch.cat([padding, freq_history], dim=1)
        
        # 提取频谱幅度信息
        freq_magnitude = torch.abs(freq_history)
        
        # GRU处理序列
        self.change_detector.flatten_parameters()
        lstm_out, _ = self.change_detector(freq_magnitude)
        
        # 提取最后时间步特征
        latest_features = lstm_out[:, -1, :]
        
        # 计算疲劳分数 (0-1)
        fatigue_scores = self.fatigue_classifier(latest_features)
        
        # 检测频谱衰减 (疲劳信号)
        if freq_history.shape[1] >= 2:
            first_window = torch.abs(freq_history[:, 0, :])
            last_window = torch.abs(freq_history[:, -1, :])
            amplitude_change = (last_window - first_window) / (first_window + 1e-8)
            decay_mask = amplitude_change < -self.decay_threshold
            
            # 计算疲劳内容类型
            fatigue_content_types = {
                'decay_patterns': decay_mask,
                'fatigue_confidence': fatigue_scores
            }
        else:
            # 历史不足时的默认值
            fatigue_content_types = None
        
        return fatigue_scores, fatigue_content_types


class DiversityEnhancer(nn.Module):
    """增强用户低频(长尾)兴趣，提供更多样化的推荐
    
    Args:
        embedding_dim (int): 嵌入维度
        novelty_boost (float): 新颖性提升系数, 默认2.0
    """
    def __init__(self, embedding_dim, novelty_boost=2.0):
        super(DiversityEnhancer, self).__init__()
        self.embedding_dim = embedding_dim
        self.novelty_boost = novelty_boost
        
        # 长尾兴趣增强网络
        self.interest_enhancer = nn.Sequential(
            nn.Linear(embedding_dim, embedding_dim),
            nn.LeakyReLU(0.2),
            nn.Linear(embedding_dim, embedding_dim)
        )
        
    def forward(self, low_freq_repr, all_item_emb=None):
        """
        Args:
            low_freq_repr (torch.Tensor): 低频兴趣表示, shape: [B, H]
            all_item_emb (torch.Tensor, optional): 所有物品的嵌入, shape: [N, H]
                N: 物品数量, H: 嵌入维度
                
        Returns:
            torch.Tensor: 增强后的低频表示, shape: [B, H]
        """
        # 增强低频兴趣表示
        enhanced_repr = self.interest_enhancer(low_freq_repr)
        
        return enhanced_repr


class DynamicFrequencyFusion(nn.Module):
    """动态融合高频和低频兴趣表示
    
    Args:
        embedding_dim (int): 嵌入维度
    """
    def __init__(self, embedding_dim):
        super(DynamicFrequencyFusion, self).__init__()
        self.embedding_dim = embedding_dim
        
        # 融合网络
        self.fusion_network = nn.Sequential(
            nn.Linear(embedding_dim * 2 + 1, embedding_dim),
            nn.LayerNorm(embedding_dim),
            nn.ReLU(),
            nn.Dropout(0.1),
            nn.Linear(embedding_dim, 64),
            nn.ReLU(),
            nn.Linear(64, 1),
            nn.Sigmoid()
        )
        self.residual_projection = nn.Linear(embedding_dim * 2, embedding_dim)
    
        
    def forward(self, high_freq_repr, low_freq_repr, fatigue_scores):
        """
        Args:
            high_freq_repr (torch.Tensor): 高频兴趣表示, shape: [B, H]
            low_freq_repr (torch.Tensor): 低频兴趣表示, shape: [B, H]
            fatigue_scores (torch.Tensor): 疲劳度分数, shape: [B, 1]
            
        Returns:
            tuple: 融合表示, 高频权重, 低频权重
        """
        batch_size = high_freq_repr.shape[0]
        
        # 拼接高低频表示和疲劳度
        fusion_input = torch.cat([
            high_freq_repr, 
            low_freq_repr, 
            fatigue_scores
        ], dim=1)
        
        # 计算动态权重 (高频权重)
        high_freq_weight = self.fusion_network(fusion_input)
        low_freq_weight = 1 - high_freq_weight

        weighted_repr = high_freq_weight * high_freq_repr + low_freq_weight * low_freq_repr
        combined_repr = torch.cat([high_freq_repr, low_freq_repr], dim=1)
        residual_repr = self.residual_projection(combined_repr)
        
        # 融合表示
        fused_repr = weighted_repr + residual_repr
        
        return fused_repr, high_freq_weight, low_freq_weight


class FreqDualDynamic(SequentialRecommender):
    """频域双通道动态推荐模型
    
    基于频域分析的推荐模型，将用户兴趣分解为高频(常见兴趣)和低频(独特兴趣)，
    动态调整它们的权重，减少重复推荐。

    Args:
        config (Config): 全局配置对象
        dataset (Dataset): 数据集
    """
    def __init__(self, config, dataset):
        super(FreqDualDynamic, self).__init__(config, dataset)
        
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
        self.cutoff_ratio = config['cutoff_ratio']  # 高低频分界点比例
        self.window_size = config['window_size']  # 频谱历史窗口大小
        self.decay_threshold = config['decay_threshold']  # 频谱衰减阈值
        self.novelty_boost = config['novelty_boost']  # 新颖性提升系数
        self.contrast_weight = config['contrast_weight']  # 对比学习权重
        self.freq_history_length = config['freq_history_length']  # 频谱历史长度
        
        # 初始化嵌入层
        self.item_embedding = nn.Embedding(
            self.n_items, self.embedding_size, padding_idx=0
        )
        self.position_embedding = nn.Embedding(
            self.max_seq_length, self.embedding_size
        )
        

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
        self.freq_decoupler = FrequencyDomainDecoupler(
            self.embedding_size, 
            self.cutoff_ratio
        )
        self.evolution_tracker = FrequencyEvolutionTracker(
            self.embedding_size, 
            self.window_size, 
            self.decay_threshold
        )
        self.diversity_enhancer = DiversityEnhancer(
            self.embedding_size, 
            self.novelty_boost
        )
        self.fusion_layer = DynamicFrequencyFusion(self.embedding_size)
        
        # 层归一化与dropout
        self.LayerNorm = nn.LayerNorm(self.embedding_size, eps=self.layer_norm_eps)
        self.dropout = nn.Dropout(self.hidden_dropout_prob)

        
        self.loss_type = config['loss_type']
        # 损失函数
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
        """
        Args:
            item_seq (torch.LongTensor): 物品序列, shape: [B, L]
            item_seq_len (torch.LongTensor): 序列长度, shape: [B]
            
        Returns:
            torch.Tensor: 用户表示, shape: [B, H]
        """
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
            self.freq_history = torch.cat([self.freq_history, freq_domain.detach()], dim=1)
            if self.freq_history.shape[1] > self.freq_history_length:
                self.freq_history = self.freq_history[:, 1:]
        
        # 跟踪兴趣演化
        fatigue_scores, _ = self.evolution_tracker(self.freq_history)
        
        # 增强低频兴趣多样性
        enhanced_low_freq = self.diversity_enhancer(low_freq_repr)
        
        # 动态融合高低频表示
        fused_repr, _, _ = self.fusion_layer(high_freq_repr, enhanced_low_freq, fatigue_scores)
        
        return fused_repr
    
    def calculate_freq_contrastive_loss(self, high_freq_repr, low_freq_repr, pos_items_emb, neg_items_emb=None):
        """计算频域对比学习损失
        
        Args:
            high_freq_repr (torch.Tensor): 高频表示, shape: [B, H]
            low_freq_repr (torch.Tensor): 低频表示, shape: [B, H]
            pos_items_emb (torch.Tensor): 正样本嵌入, shape: [B, H]
            neg_items_emb (torch.Tensor, optional): 负样本嵌入, shape: [B, H]
            
        Returns:
            torch.Tensor: 对比学习损失
        """
        batch_size = high_freq_repr.shape[0]
        
        # 计算高频表示和低频表示与正样本的相似度
        high_freq_pos_sim = torch.sum(high_freq_repr * pos_items_emb, dim=1)
        low_freq_pos_sim = torch.sum(low_freq_repr * pos_items_emb, dim=1)
        
        # 鼓励高频表示与正样本更相似
        high_freq_target = torch.ones_like(high_freq_pos_sim)
        high_freq_loss = F.mse_loss(torch.sigmoid(high_freq_pos_sim), high_freq_target)
        
        # 如果有负样本，进一步增强对比
        if neg_items_emb is not None:
            high_freq_neg_sim = torch.sum(high_freq_repr * neg_items_emb, dim=1)
            low_freq_neg_sim = torch.sum(low_freq_repr * neg_items_emb, dim=1)
            
            # 鼓励低频表示更好地区分正负样本
            low_freq_logits = torch.stack([low_freq_pos_sim, low_freq_neg_sim], dim=1)
            low_freq_labels = torch.zeros(batch_size, dtype=torch.long, device=low_freq_repr.device)
            low_freq_loss = self.contrast_loss(low_freq_logits, low_freq_labels)
            
            # 合并两种损失
            contrastive_loss = high_freq_loss + low_freq_loss
        else:
            # 无负样本时，简化为仅使用高频损失
            contrastive_loss = high_freq_loss
        
        return contrastive_loss
    
    def calculate_loss(self, interaction):
        """计算模型损失
        
        Args:
            interaction (Interaction): 交互数据
            
        Returns:
            torch.Tensor: 损失值
        """
        item_seq = interaction[self.ITEM_SEQ]
        item_seq_len = interaction[self.ITEM_SEQ_LEN]
        pos_items = interaction[self.POS_ITEM_ID]
        
        seq_output = self.forward(item_seq, item_seq_len)
        
        # 频域解耦表示 (重新计算是为了训练频域模块)
        high_freq_repr, low_freq_repr, _ = self.freq_decoupler(seq_output)
        
        # 计算主要推荐损失
        if self.loss_type == 'BPR':
            neg_items = interaction[self.NEG_ITEM_ID]
            pos_items_emb = self.item_embedding(pos_items)
            neg_items_emb = self.item_embedding(neg_items)
            pos_score = torch.sum(seq_output * pos_items_emb, dim=-1)
            neg_score = torch.sum(seq_output * neg_items_emb, dim=-1)
            main_loss = self.loss_fct(pos_score, neg_score)
            
            # 频域对比学习损失
            contrastive_loss = self.calculate_freq_contrastive_loss(
                high_freq_repr, low_freq_repr, pos_items_emb, neg_items_emb
            )
        else:  # CE损失
            test_item_emb = self.item_embedding.weight
            logits = torch.matmul(seq_output, test_item_emb.transpose(0, 1))
            main_loss = self.loss_fct(logits, pos_items)
            
            # 频域对比学习损失
            pos_items_emb = self.item_embedding(pos_items)
            contrastive_loss = self.calculate_freq_contrastive_loss(
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
        
        return scores
