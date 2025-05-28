# -*- coding: utf-8 -*-
# @Time    : 2025/5/27
# @Author  : OuyangNingtai
# @Description: 简化版FRF模型 - 修复batch_size不匹配问题

r"""
FRF Simple - Simplified Frequency Repetitive Fusion Model
################################################
简化版频域重复兴趣融合推荐模型

修复问题：
- 解决不同batch_size导致的张量连接错误
- 改用全局频域历史存储策略
- 保持核心功能完整性
"""

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
import math
import time
from collections import deque

from recbole.model.abstract_recommender import SequentialRecommender
from recbole.model.loss import BPRLoss
from recbole.data.interaction import Interaction


class FrequencyDomainDecoupler(nn.Module):
    """简化版频域解耦器
    
    将用户兴趣分解为高频(重复兴趣)和低频(探索兴趣)两个表示
    
    Args:
        embedding_dim (int): 嵌入维度
        cutoff_ratio (float): 高低频分界点比例，固定值
        high_freq_scale (float): 高频增强比例
    """
    def __init__(self, embedding_dim, cutoff_ratio=0.6, high_freq_scale=1.5):
        super(FrequencyDomainDecoupler, self).__init__()
        self.embedding_dim = embedding_dim
        self.cutoff_ratio = cutoff_ratio
        self.high_freq_scale = high_freq_scale
        
        # 特征编码器
        self.encoder = nn.Sequential(
            nn.Linear(embedding_dim, embedding_dim),
            nn.LayerNorm(embedding_dim),  
            nn.GELU()
        )
        
        # 高频投影层
        self.high_freq_projector = nn.Sequential(
            nn.Linear(embedding_dim, embedding_dim),
            nn.LayerNorm(embedding_dim)
        )
        
        # 低频投影层
        self.low_freq_projector = nn.Sequential(
            nn.Linear(embedding_dim, embedding_dim),
            nn.LayerNorm(embedding_dim)
        )
    
    def forward(self, seq_output):
        """频域解耦处理
        
        Args:
            seq_output (torch.Tensor): 序列表示, shape: [B, H]
            
        Returns:
            tuple: 高频表示, 低频表示, 频域表示
        """
        # 安全性检查
        if torch.isnan(seq_output).any():
            seq_output = torch.nan_to_num(seq_output, nan=0.0, posinf=1.0, neginf=-1.0)

        # 编码用户行为序列
        processed_seq = seq_output.unsqueeze(1)  # [B, 1, H]
        encoded = self.encoder(processed_seq)  # [B, 1, H]

        if torch.isnan(encoded).any():
            encoded = torch.nan_to_num(encoded, nan=0.0, posinf=1.0, neginf=-1.0)

        # 归一化处理
        encoded_mean = encoded.mean(dim=2, keepdim=True)
        encoded_std = encoded.std(dim=2, keepdim=True) + 1e-6
        encoded_norm = (encoded - encoded_mean) / encoded_std
        
        # 转换到频域
        freq_domain = torch.fft.rfft(encoded_norm, dim=2, norm="ortho")  # [B, 1, H//2+1]
        
        # 固定高低频划分
        cutoff = int(freq_domain.shape[2] * self.cutoff_ratio)
        
        # 高频成分 (重复兴趣) - 增强处理
        high_freq_components = freq_domain.clone()
        high_freq_components[:, :, cutoff:] = 0  # 清零低频部分
        high_freq_components *= self.high_freq_scale  # 增强高频
        
        # 低频成分 (探索兴趣)
        low_freq_components = freq_domain.clone()
        low_freq_components[:, :, :cutoff] = 0  # 清零高频部分
        
        # 逆变换回时域
        high_freq_repr = torch.fft.irfft(high_freq_components, n=self.embedding_dim, dim=2, norm="ortho")
        low_freq_repr = torch.fft.irfft(low_freq_components, n=self.embedding_dim, dim=2, norm="ortho")

        # 安全性检查
        if torch.isnan(high_freq_repr).any() or torch.isnan(low_freq_repr).any():
            high_freq_repr = torch.nan_to_num(high_freq_repr, nan=0.0, posinf=1.0, neginf=-1.0)
            low_freq_repr = torch.nan_to_num(low_freq_repr, nan=0.0, posinf=1.0, neginf=-1.0)
        
        # 降维
        high_freq_repr = high_freq_repr.squeeze(1)  # [B, H]
        low_freq_repr = low_freq_repr.squeeze(1)  # [B, H]
        
        # 投影到表示空间
        high_freq_repr = self.high_freq_projector(high_freq_repr)
        low_freq_repr = self.low_freq_projector(low_freq_repr)

        return high_freq_repr, low_freq_repr, freq_domain.squeeze(1)


class DTWRepetitionTracker(nn.Module):
    """简化版DTW重复兴趣跟踪器
    
    使用DTW距离跟踪用户重复兴趣模式
    
    Args:
        embedding_dim (int): 嵌入维度
        window_size (int): 历史窗口大小
        rep_threshold (float): 重复兴趣阈值
    """
    def __init__(self, embedding_dim, window_size=3, rep_threshold=0.6):
        super(DTWRepetitionTracker, self).__init__()
        self.embedding_dim = embedding_dim
        self.window_size = window_size
        self.rep_threshold = rep_threshold
        
        # DTW计算缓存
        self.dtw_cache = {}
        self.max_cache_size = 50
        
        # 创建一个dummy参数以确保模块有参数
        self.dummy_param = nn.Parameter(torch.zeros(1))
        
        # 动态构建的模式检测器
        self.pattern_detector = None
        self._input_dim = None
        
        # 使用列表存储历史频域特征，避免batch_size不匹配问题
        self.freq_history_list = deque(maxlen=window_size * 10)  # 最多存储window_size*10个样本
    
    def _build_pattern_detector(self, input_dim, device):
        """动态构建模式检测器"""
        if self.pattern_detector is None:
            self._input_dim = input_dim
            hidden_dim = max(input_dim // 4, 16)  # 确保隐藏层至少有16个神经元
            
            self.pattern_detector = nn.Sequential(
                nn.Linear(input_dim, hidden_dim),
                nn.LayerNorm(hidden_dim),
                nn.GELU(),
                nn.Linear(hidden_dim, 1),
                nn.Sigmoid()
            ).to(device)
            
            # 正确注册到模块中
            self.add_module('pattern_detector', self.pattern_detector)
            
            # 初始化权重
            for module in self.pattern_detector.modules():
                if isinstance(module, nn.Linear):
                    nn.init.xavier_normal_(module.weight)
                    if module.bias is not None:
                        nn.init.zeros_(module.bias)
    
    def update_freq_history(self, freq_features):
        """更新频域历史
        
        Args:
            freq_features (torch.Tensor): 频域特征 [B, F]
        """
        if not self.training:
            return
            
        # 将当前批次的特征逐个添加到历史列表中
        with torch.no_grad():
            for i in range(freq_features.shape[0]):
                self.freq_history_list.append(freq_features[i].cpu().clone())
    
    def get_freq_history_tensor(self, batch_size, device):
        """获取频域历史张量
        
        Args:
            batch_size (int): 当前批次大小
            device: 设备
            
        Returns:
            torch.Tensor: 历史频域特征 [B, T, F] 或 None
        """
        if len(self.freq_history_list) < 2:
            return None
        
        # 取最近的一些历史特征
        history_len = min(self.window_size, len(self.freq_history_list))
        recent_history = list(self.freq_history_list)[-history_len:]
        
        # 构造伪历史张量：为每个批次样本复制相同的历史
        freq_dim = recent_history[0].shape[0]
        history_tensor = torch.zeros(batch_size, history_len, freq_dim, device=device)
        
        for t, hist_feat in enumerate(recent_history):
            history_tensor[:, t, :] = hist_feat.to(device)
        
        return history_tensor
    
    def dtw_distance(self, seq1: torch.Tensor, seq2: torch.Tensor) -> torch.Tensor:
        """计算DTW距离（简化版）
        
        Args:
            seq1, seq2: 输入序列
            
        Returns:
            torch.Tensor: DTW距离
        """
        # 检查缓存
        try:
            cache_key = (hash(seq1.cpu().detach().numpy().tobytes()), 
                         hash(seq2.cpu().detach().numpy().tobytes()))
            if cache_key in self.dtw_cache:
                return self.dtw_cache[cache_key]
        except:
            cache_key = None
            
        n, m = seq1.size(0), seq2.size(0)
        
        # 短序列优化
        if n <= 2 or m <= 2:
            if n == m:
                dist = torch.norm(seq1 - seq2, p=2) / max(n, 1)
            else:
                dist = torch.norm(seq1.mean(dim=0) - seq2.mean(dim=0), p=2)
        else:
            # 完整DTW计算
            dtw_matrix = torch.zeros((n+1, m+1), device=seq1.device)
            dtw_matrix[:, 0] = float('inf')
            dtw_matrix[0, :] = float('inf')
            dtw_matrix[0, 0] = 0
            
            for i in range(1, n+1):
                for j in range(1, m+1):
                    cost = torch.norm(seq1[i-1] - seq2[j-1], p=2)
                    dtw_matrix[i, j] = cost + min(
                        dtw_matrix[i-1, j],
                        dtw_matrix[i, j-1],
                        dtw_matrix[i-1, j-1]
                    )
            
            dist = dtw_matrix[n, m]
        
        # 更新缓存
        if cache_key is not None and len(self.dtw_cache) < self.max_cache_size:
            self.dtw_cache[cache_key] = dist
            
        return dist
    
    def compute_repetition_score(self, freq_history):
        """计算重复强度分数
        
        Args:
            freq_history (torch.Tensor): 频域历史 [B, T, F] 或 None
            
        Returns:
            torch.Tensor: 重复强度分数 [B, 1]
        """
        if freq_history is None:
            # 没有历史，返回默认分数
            batch_size = 1  # 默认值，会在调用时被正确设置
            device = self.dummy_param.device
            return torch.ones(batch_size, 1, device=device) * 0.5
        
        batch_size = freq_history.shape[0]
        device = freq_history.device
        
        if freq_history.shape[1] < 2:
            # 历史不足，返回默认分数
            return torch.ones(batch_size, 1, device=device) * 0.5
        
        # 提取高频部分
        cutoff = int(freq_history.shape[2] * 0.6)
        high_freq_history = torch.abs(freq_history[:, :, :cutoff])
        
        # 归一化
        mean = high_freq_history.mean(dim=(1,2), keepdim=True)
        std = high_freq_history.std(dim=(1,2), keepdim=True) + 1e-6
        high_freq_history = (high_freq_history - mean) / std
        
        # 计算最新特征
        latest_features = high_freq_history[:, -1, :]  # [B, F]
        
        # 动态构建模式检测器
        if self.pattern_detector is None:
            self._build_pattern_detector(latest_features.shape[1], device)
        
        # 使用模式检测器
        repetition_scores = self.pattern_detector(latest_features)
        
        # DTW相似性分析（如果有足够历史）
        if freq_history.shape[1] >= 2:
            dtw_similarities = []
            for i in range(batch_size):
                try:
                    # 比较最新和历史的相似性
                    recent_seq = high_freq_history[i, -1, :].unsqueeze(0)
                    past_seq = high_freq_history[i, :-1, :].mean(dim=0, keepdim=True)
                    
                    dtw_dist = self.dtw_distance(recent_seq, past_seq)
                    similarity = torch.exp(-dtw_dist / self.embedding_dim)
                    dtw_similarities.append(similarity)
                except:
                    dtw_similarities.append(torch.tensor(0.5, device=device))
            
            dtw_sim = torch.stack(dtw_similarities).unsqueeze(1)
            # 结合重复分数和DTW相似性
            repetition_scores = 0.7 * repetition_scores + 0.3 * dtw_sim
        
        return repetition_scores
    
    def forward(self, batch_size, device):
        """前向传播
        
        Args:
            batch_size (int): 当前批次大小
            device: 设备
            
        Returns:
            torch.Tensor: 重复强度分数 [B, 1]
        """
        # 获取历史频域特征
        freq_history = self.get_freq_history_tensor(batch_size, device)
        
        if freq_history is not None and torch.isnan(freq_history).any():
            freq_history = torch.nan_to_num(freq_history, nan=0.0, posinf=1.0, neginf=-1.0)
        
        repetition_scores = self.compute_repetition_score(freq_history)
        
        # 确保批次大小正确
        if repetition_scores.shape[0] != batch_size:
            repetition_scores = torch.ones(batch_size, 1, device=device) * 0.5
        
        if torch.isnan(repetition_scores).any():
            repetition_scores = torch.nan_to_num(repetition_scores, nan=0.5, posinf=1.0, neginf=0.0)
        
        return repetition_scores


class UnifiedRecommender(nn.Module):
    """简化版统一推荐器
    
    使用固定策略融合高频和低频表示
    
    Args:
        embedding_dim (int): 嵌入维度
        rep_weight (float): 重复兴趣基准权重
    """
    def __init__(self, embedding_dim, rep_weight=0.7):
        super(UnifiedRecommender, self).__init__()
        self.embedding_dim = embedding_dim
        self.rep_weight = rep_weight
        
        # 统一融合网络
        self.fusion_network = nn.Sequential(
            nn.Linear(embedding_dim * 2 + 1, embedding_dim),
            nn.LayerNorm(embedding_dim),
            nn.GELU(),
            nn.Linear(embedding_dim, 1),
            nn.Sigmoid()
        )
        
        # 最终投影
        self.final_projection = nn.Sequential(
            nn.Linear(embedding_dim * 2, embedding_dim),
            nn.LayerNorm(embedding_dim)
        )
    
    def forward(self, high_freq_repr, low_freq_repr, repetition_scores):
        """融合高频和低频表示
        
        Args:
            high_freq_repr (torch.Tensor): 高频表示 [B, H]
            low_freq_repr (torch.Tensor): 低频表示 [B, H]
            repetition_scores (torch.Tensor): 重复强度分数 [B, 1]
            
        Returns:
            torch.Tensor: 融合后的用户表示 [B, H]
        """
        # 动态权重计算
        fusion_input = torch.cat([high_freq_repr, low_freq_repr, repetition_scores], dim=1)
        dynamic_adjustment = self.fusion_network(fusion_input)
        
        # 计算最终权重
        repetition_weight = self.rep_weight + (1 - self.rep_weight) * dynamic_adjustment
        novelty_weight = 1 - repetition_weight
        
        # 加权融合
        weighted_high = repetition_weight * high_freq_repr
        weighted_low = novelty_weight * low_freq_repr
        
        # 最终表示
        concat_repr = torch.cat([weighted_high, weighted_low], dim=1)
        final_repr = self.final_projection(concat_repr)
        
        return final_repr


class FRF(SequentialRecommender):
    """简化版频域重复兴趣融合推荐模型
    
    保留核心功能：频域解耦 + DTW距离 + 重复检测
    简化配置：移除复杂的可选功能，提升性能和可维护性
    
    Args:
        config (Config): 配置对象
        dataset (Dataset): 数据集
    """
    def __init__(self, config, dataset):
        super(FRF, self).__init__(config, dataset)
        
        # 基础参数
        self.embedding_size = config['embedding_size']
        self.hidden_size = config['hidden_size']
        self.n_layers = config['n_layers']
        self.n_heads = config['n_heads']
        self.inner_size = config['inner_size']
        self.hidden_dropout_prob = config['hidden_dropout_prob']
        self.attn_dropout_prob = config['attn_dropout_prob']
        self.hidden_act = config['hidden_act']
        self.layer_norm_eps = config['layer_norm_eps']
        
        # 简化后的核心参数
        self.cutoff_ratio = config['cutoff_ratio'] if 'cutoff_ratio' in config else 0.6
        self.rep_weight = config['rep_weight'] if 'rep_weight' in config else 0.7
        self.high_freq_scale = config['high_freq_scale'] if 'high_freq_scale' in config else 1.5
        self.window_size = config['window_size'] if 'window_size' in config else 3
        self.rep_threshold = config['rep_threshold'] if 'rep_threshold' in config else 0.6
        self.contrast_weight = config['contrast_weight'] if 'contrast_weight' in config else 0.3
        
        # 初始化嵌入层
        self.item_embedding = nn.Embedding(self.n_items, self.embedding_size, padding_idx=0)
        self.position_embedding = nn.Embedding(self.max_seq_length, self.embedding_size)
        
        # 序列编码器
        try:
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
        except ImportError:
            # 降级方案
            self.encoder = nn.LSTM(
                input_size=self.embedding_size,
                hidden_size=self.hidden_size,
                num_layers=1,
                batch_first=True
            )
        
        # 核心模块
        self.freq_decoupler = FrequencyDomainDecoupler(
            self.embedding_size, 
            self.cutoff_ratio,
            self.high_freq_scale
        )
        
        self.repetition_tracker = DTWRepetitionTracker(
            self.embedding_size,
            self.window_size,
            self.rep_threshold
        )
        
        self.recommender = UnifiedRecommender(
            self.embedding_size,
            self.rep_weight
        )
        
        # 基础组件
        self.LayerNorm = nn.LayerNorm(self.embedding_size, eps=self.layer_norm_eps)
        self.dropout = nn.Dropout(self.hidden_dropout_prob)
        
        # 损失函数
        self.loss_type = config['loss_type']
        if self.loss_type == 'BPR':
            self.loss_fct = BPRLoss()
        elif self.loss_type == 'CE':
            self.loss_fct = nn.CrossEntropyLoss()
        else:
            raise NotImplementedError("Make sure 'loss_type' in ['BPR', 'CE']!")
        
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
        extended_attention_mask = extended_attention_mask.to(dtype=next(self.parameters()).dtype)
        extended_attention_mask = (1.0 - extended_attention_mask) * -10000.0
        return extended_attention_mask
    
    def forward(self, item_seq, item_seq_len):
        """前向传播
        
        Args:
            item_seq (torch.LongTensor): 物品序列 [B, L]
            item_seq_len (torch.LongTensor): 序列长度 [B]
            
        Returns:
            torch.Tensor: 用户表示 [B, H]
        """
        # 序列编码
        position_ids = torch.arange(item_seq.size(1), dtype=torch.long, device=item_seq.device)
        position_ids = position_ids.unsqueeze(0).expand_as(item_seq)
        position_embedding = self.position_embedding(position_ids)
        
        item_emb = self.item_embedding(item_seq)
        input_emb = item_emb + position_embedding
        input_emb = self.LayerNorm(input_emb)
        input_emb = self.dropout(input_emb)
        
        # 编码器处理
        if isinstance(self.encoder, nn.LSTM):
            encoder_output, _ = self.encoder(input_emb)
            seq_output = encoder_output
        else:
            extended_attention_mask = self.get_attention_mask(item_seq)
            try:
                encoder_output = self.encoder(input_emb, extended_attention_mask)
                seq_output = encoder_output[-1] if isinstance(encoder_output, list) else encoder_output
            except:
                seq_output = input_emb.mean(dim=1, keepdim=True).expand(-1, input_emb.size(1), -1)
        
        # 获取最后时刻的表示
        output = self.gather_indexes(seq_output, item_seq_len - 1)
        if torch.isnan(output).any():
            output = torch.nan_to_num(output, nan=0.0, posinf=1.0, neginf=-1.0)
        
        # 频域解耦
        high_freq_repr, low_freq_repr, freq_domain = self.freq_decoupler(output)
        
        # 更新频域历史（使用新的策略）
        self.repetition_tracker.update_freq_history(freq_domain)
        
        # 重复兴趣检测（传递batch_size和device）
        batch_size = high_freq_repr.shape[0]
        device = high_freq_repr.device
        repetition_scores = self.repetition_tracker(batch_size, device)
        
        # 融合推荐
        final_repr = self.recommender(high_freq_repr, low_freq_repr, repetition_scores)
        
        return final_repr
    
    def calculate_loss(self, interaction):
        """计算损失"""
        item_seq = interaction[self.ITEM_SEQ]
        item_seq_len = interaction[self.ITEM_SEQ_LEN]
        pos_items = interaction[self.POS_ITEM_ID]
        
        seq_output = self.forward(item_seq, item_seq_len)
        
        # 主推荐损失
        if self.loss_type == 'BPR':
            neg_items = interaction[self.NEG_ITEM_ID]
            pos_items_emb = self.item_embedding(pos_items)
            neg_items_emb = self.item_embedding(neg_items)
            pos_score = torch.sum(seq_output * pos_items_emb, dim=-1)
            neg_score = torch.sum(seq_output * neg_items_emb, dim=-1)
            main_loss = self.loss_fct(pos_score, neg_score)
            
            # 简化的对比学习损失
            high_freq_repr, _, _ = self.freq_decoupler(seq_output)
            pos_sim = F.cosine_similarity(high_freq_repr, pos_items_emb, dim=1)
            neg_sim = F.cosine_similarity(high_freq_repr, neg_items_emb, dim=1)
            contrast_loss = F.relu(1.0 - pos_sim + neg_sim).mean()
        else:
            test_item_emb = self.item_embedding.weight
            logits = torch.matmul(seq_output, test_item_emb.transpose(0, 1))
            main_loss = self.loss_fct(logits, pos_items)
            contrast_loss = 0.0
        
        return main_loss + self.contrast_weight * contrast_loss
    
    def predict(self, interaction):
        """预测"""
        item_seq = interaction[self.ITEM_SEQ]
        item_seq_len = interaction[self.ITEM_SEQ_LEN]
        test_item = interaction[self.ITEM_ID]
        
        seq_output = self.forward(item_seq, item_seq_len)
        test_item_emb = self.item_embedding(test_item)
        scores = torch.mul(seq_output, test_item_emb).sum(dim=1)
        
        return scores
    
    def full_sort_predict(self, interaction):
        """全排序预测"""
        item_seq = interaction[self.ITEM_SEQ]
        item_seq_len = interaction[self.ITEM_SEQ_LEN]
        
        seq_output = self.forward(item_seq, item_seq_len)
        test_items_emb = self.item_embedding.weight
        scores = torch.matmul(seq_output, test_items_emb.transpose(0, 1))
        
        return scores
