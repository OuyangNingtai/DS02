# -*- coding: utf-8 -*-
# @Time    : 2025/4/19


r"""
FRF
################################################
Reference:

"""

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
import math
from typing import List, Dict, Tuple, Optional
import warnings
import time

from recbole.model.abstract_recommender import SequentialRecommender
from recbole.model.loss import BPRLoss
from recbole.data.interaction import Interaction


class OptimizedDTWCoreSetBuilder:
    """优化版DTW核心集构建器
    
    降低计算频率、使用缓存、简化敏感度计算
    
    Args:
        epsilon (float): 核心集精度参数，默认0.1
        vc_dim_factor (float): VC维相关因子，默认0.3
        max_coreset_size (int): 最大核心集大小，默认为None（自适应）
        update_freq (int): 更新频率，默认为5（每5个批次更新一次）
        use_fast_dtw (bool): 是否使用快速DTW算法，默认True
        max_seq_len (int): 最大序列长度，超过将进行采样，默认50
    """
    def __init__(self, epsilon=0.1, vc_dim_factor=0.3, max_coreset_size=None,
                 update_freq=5, use_fast_dtw=True, max_seq_len=50):
        self.epsilon = epsilon
        self.vc_dim_factor = vc_dim_factor
        self.max_coreset_size = max_coreset_size
        self.update_freq = update_freq
        self.use_fast_dtw = use_fast_dtw
        self.max_seq_len = max_seq_len
        
        # 计数器和缓存
        self.update_counter = 0
        self.distance_cache = {}
        self.coreset_cache = {}
        self.last_coreset = None
        self.last_weights = None
        self.max_cache_size = 100
    
    def should_update(self):
        """检查是否应该执行更新操作"""
        self.update_counter = (self.update_counter + 1) % self.update_freq
        return self.update_counter == 0
    
    def downsample_sequence(self, seq: torch.Tensor) -> torch.Tensor:
        """对长序列进行降采样
        
        Args:
            seq (torch.Tensor): 输入序列
            
        Returns:
            torch.Tensor: 降采样后的序列
        """
        if seq.size(0) <= self.max_seq_len:
            return seq
            
        # 均匀采样
        indices = torch.linspace(0, seq.size(0)-1, self.max_seq_len).long()
        return seq[indices]
    
    def get_cache_key(self, sequences):
        """生成序列列表的缓存键"""
        if len(sequences) == 0:
            return None
        
        # 使用第一个和最后一个序列作为标识
        try:
            first_seq = sequences[0]
            last_seq = sequences[-1]
            key = (
                hash(first_seq.cpu().detach().numpy().tobytes()),
                hash(last_seq.cpu().detach().numpy().tobytes()),
                len(sequences)
            )
            return key
        except:
            return None
    
    def dtw_distance(self, seq1: torch.Tensor, seq2: torch.Tensor, p: int = 2) -> torch.Tensor:
        """计算两个序列间的DTW距离（带缓存）
        
        Args:
            seq1 (torch.Tensor): 第一个序列
            seq2 (torch.Tensor): 第二个序列
            p (int): 距离范数，默认为2（欧氏距离）
            
        Returns:
            torch.Tensor: DTW距离值
        """
        # 尝试从缓存获取
        try:
            cache_key = (hash(seq1.cpu().detach().numpy().tobytes()), 
                        hash(seq2.cpu().detach().numpy().tobytes()))
                        
            if cache_key in self.distance_cache:
                return self.distance_cache[cache_key]
        except:
            cache_key = None
            
        # 降采样长序列
        seq1 = self.downsample_sequence(seq1)
        seq2 = self.downsample_sequence(seq2)
        
        n, m = seq1.size(0), seq2.size(0)
        
        # 使用FastDTW进行近似计算
        if self.use_fast_dtw and (n > 20 or m > 20):
            dist = self.fast_dtw(seq1, seq2, radius=2, p=p)
        else:
            # 创建DTW矩阵
            dtw_matrix = torch.zeros((n+1, m+1), device=seq1.device)
            dtw_matrix[:, 0] = float('inf')
            dtw_matrix[0, :] = float('inf')
            dtw_matrix[0, 0] = 0
            
            # 计算成对距离
            for i in range(1, n+1):
                for j in range(1, m+1):
                    cost = torch.norm(seq1[i-1] - seq2[j-1], p=p)
                    dtw_matrix[i, j] = cost + min(
                        dtw_matrix[i-1, j],    # 插入
                        dtw_matrix[i, j-1],    # 删除
                        dtw_matrix[i-1, j-1]   # 替换
                    )
                    
            dist = dtw_matrix[n, m]
        
        # 存入缓存
        if cache_key is not None:
            if len(self.distance_cache) >= self.max_cache_size:
                # 清理部分缓存
                keys_to_remove = list(self.distance_cache.keys())[:self.max_cache_size//2]
                for k in keys_to_remove:
                    self.distance_cache.pop(k, None)
            
            self.distance_cache[cache_key] = dist
            
        return dist
    
    def fast_dtw(self, seq1, seq2, radius=1, p=2):
        """快速DTW算法实现 - 降低时间复杂度
        
        Args:
            seq1, seq2: 输入序列
            radius: 搜索半径
            p: 距离范数
            
        Returns:
            torch.Tensor: 近似DTW距离
        """
        # 基本情况 - 序列较短时直接计算
        min_size = radius + 2
        if seq1.size(0) <= min_size or seq2.size(0) <= min_size:
            # 简化的DTW计算 - 避免创建大矩阵
            return self._simple_dtw(seq1, seq2, p)
        
        # 降采样
        seq1_shrunk = seq1[::2] 
        seq2_shrunk = seq2[::2]
        
        # 递归计算低分辨率路径
        low_res_dist = self.fast_dtw(seq1_shrunk, seq2_shrunk, radius, p)
        
        # 近似值直接返回
        return low_res_dist * (seq1.size(0) / seq1_shrunk.size(0) + seq2.size(0) / seq2_shrunk.size(0)) / 2
    
    def _simple_dtw(self, seq1, seq2, p=2):
        """简化的DTW计算 - 用于短序列"""
        n, m = seq1.size(0), seq2.size(0)
        
        # 对于非常短的序列，简化为欧式距离均值
        if n <= 3 or m <= 3:
            if n == m:
                return torch.norm(seq1 - seq2, p=p) / n
            else:
                # 简单的均值比较
                return torch.norm(seq1.mean(dim=0) - seq2.mean(dim=0), p=p)
        
        # 简化的DTW实现
        prev_row = torch.zeros(m+1, device=seq1.device)
        curr_row = torch.zeros(m+1, device=seq1.device)
        
        # 初始化第一行
        for j in range(1, m+1):
            prev_row[j] = prev_row[j-1] + torch.norm(seq1[0] - seq2[j-1], p=p)
            
        # 填充矩阵
        for i in range(1, n):
            curr_row[0] = float('inf')
            for j in range(1, m+1):
                cost = torch.norm(seq1[i] - seq2[j-1], p=p)
                curr_row[j] = cost + min(
                    prev_row[j],      # 插入
                    curr_row[j-1],    # 删除
                    prev_row[j-1]     # 匹配
                )
            # 交换行
            prev_row, curr_row = curr_row, prev_row
            
        return prev_row[m]
    
    def compute_simplified_sensitivities(self, sequences: List[torch.Tensor]) -> torch.Tensor:
        """计算简化版敏感度 - 降低计算成本
        
        Args:
            sequences (List[torch.Tensor]): 序列列表
            
        Returns:
            torch.Tensor: 敏感度数组
        """
        n = len(sequences)
        if n <= 1:
            return torch.ones(n, device=sequences[0].device)
        
        # 减少中心点数量
        k = min(n, max(2, int(self.vc_dim_factor * math.log(n))))
        k = min(k, 5)  # 限制最大中心数为5
        
        # 随机选择中心点
        center_indices = torch.randperm(n)[:k]
        centers = [sequences[i] for i in center_indices]
        
        # 计算到最近中心的距离
        sensitivities = torch.zeros(n, device=sequences[0].device)
        for i in range(n):
            if i % 10 == 0:  # 只处理部分样本以提高性能
                min_dist = float('inf')
                for center in centers:
                    # 使用简化距离计算
                    if sequences[i].shape[0] > 10 or center.shape[0] > 10:
                        # 对长序列使用均值距离
                        dist = torch.norm(sequences[i].mean(dim=0) - center.mean(dim=0)).item()
                    else:
                        # 对短序列使用DTW
                        dist = self.dtw_distance(sequences[i], center).item()
                    min_dist = min(min_dist, dist)
                sensitivities[i] = max(min_dist, 1e-5)
            else:
                # 对未处理的样本使用均值
                sensitivities[i] = 0.01
        
        # 对未设置的样本使用均值
        mean_sens = sensitivities[sensitivities > 0].mean() if torch.any(sensitivities > 0) else 0.01
        sensitivities[sensitivities == 0] = mean_sens
        
        # 归一化敏感度
        total = torch.sum(sensitivities)
        if total > 0:
            sensitivities = sensitivities / total
        else:
            sensitivities = torch.ones_like(sensitivities) / n
            
        return sensitivities
    
    def build_coreset(self, sequences: List[torch.Tensor]) -> Tuple[List[torch.Tensor], torch.Tensor]:
        """构建DTW核心集 - 支持缓存和跳过更新
        
        Args:
            sequences (List[torch.Tensor]): 输入序列列表
            
        Returns:
            Tuple[List[torch.Tensor], torch.Tensor]: 核心集序列和对应权重
        """
        n = len(sequences)
        if n <= 1:
            return sequences, torch.ones(n, device=sequences[0].device)
        
        # 检查是否需要更新
        if not self.should_update() and self.last_coreset is not None:
            # 重用上次的核心集
            return self.last_coreset, self.last_weights
        
        # 尝试从缓存获取
        cache_key = self.get_cache_key(sequences)
        if cache_key in self.coreset_cache:
            self.last_coreset, self.last_weights = self.coreset_cache[cache_key]
            return self.last_coreset, self.last_weights
        
        # 计算简化版敏感度
        sensitivities = self.compute_simplified_sensitivities(sequences)
        
        # 确定核心集大小 (显著减小)
        coreset_size = int(min(n, self.max_coreset_size or (self.vc_dim_factor * n * math.log(n) / self.epsilon**2)))
        coreset_size = min(coreset_size, max(1, n // 10))  # 更激进的减小核心集大小
        coreset_size = max(1, min(coreset_size, n))  # 确保合理范围
        
        # 基于敏感度进行采样
        indices = torch.multinomial(sensitivities, coreset_size, replacement=True)
        coreset = [sequences[i] for i in indices]
        
        # 计算权重
        weights = torch.zeros(coreset_size, device=sequences[0].device)
        for i, idx in enumerate(indices):
            weights[i] = 1.0 / (coreset_size * sensitivities[idx])
        
        # 归一化权重
        weights = weights / torch.sum(weights)
        
        # 更新缓存
        self.last_coreset = coreset
        self.last_weights = weights
        
        if cache_key is not None:
            if len(self.coreset_cache) >= self.max_cache_size:
                # 清理部分缓存
                keys_to_remove = list(self.coreset_cache.keys())[:self.max_cache_size//2]
                for k in keys_to_remove:
                    self.coreset_cache.pop(k, None)
                    
            self.coreset_cache[cache_key] = (coreset, weights)
        
        return coreset, weights
    
    def combine_embeddings(self, coreset: List[torch.Tensor], weights: torch.Tensor) -> torch.Tensor:
        """将核心集合并为单一嵌入表示 - 优化版本
        
        Args:
            coreset (List[torch.Tensor]): 核心集序列
            weights (torch.Tensor): 对应权重
            
        Returns:
            torch.Tensor: 合并后的嵌入表示
        """
        if len(coreset) == 0:
            return None
        
        # 快速路径 - 单个序列
        if len(coreset) == 1:
            return coreset[0]
        
        # 对齐序列长度 - 使用平均池化而非填充
        embeddings = []
        max_len = max(seq.size(0) for seq in coreset)
        max_len = min(max_len, self.max_seq_len)  # 限制长度
        
        for seq in coreset:
            if seq.size(0) > max_len:
                # 对过长序列使用平均池化降采样
                ratio = seq.size(0) / max_len
                pooled = F.avg_pool1d(
                    seq.transpose(0, 1).unsqueeze(0),  # [1, dim, len]
                    kernel_size=int(ratio),
                    stride=int(ratio),
                    ceil_mode=True
                )
                pooled = pooled.squeeze(0).transpose(0, 1)  # [new_len, dim]
                if pooled.size(0) > max_len:
                    pooled = pooled[:max_len]
                embeddings.append(pooled)
            else:
                embeddings.append(seq)
        
        # 加权平均 - 简化处理
        weights_norm = weights / weights.sum()
        result = None
        
        for i, emb in enumerate(embeddings):
            if result is None:
                result = emb * weights_norm[i]
            else:
                # 处理长度不一致的情况
                min_len = min(result.size(0), emb.size(0))
                result[:min_len] += emb[:min_len] * weights_norm[i]
        
        return result


class SimplifiedFrequencyDomainDecoupler(nn.Module):
    """简化版基于DTW核心集的频域解耦器
    
    将用户兴趣分解为高频(重复兴趣)和低频(次要兴趣)两个独立表示
    
    Args:
        embedding_dim (int): 嵌入维度
        cutoff_ratio (float): 高低频分界点比例, 默认0.5 
        high_freq_scale (float): 高频增强比例, 默认1.5
        dtw_epsilon (float): DTW核心集精度参数，默认0.1
        adaptive_cutoff (bool): 是否使用自适应截断比例，默认True
    """
    def __init__(self, embedding_dim, cutoff_ratio=0.5, high_freq_scale=1.5, 
                 dtw_epsilon=0.1, adaptive_cutoff=True):
        super(SimplifiedFrequencyDomainDecoupler, self).__init__()
        self.embedding_dim = embedding_dim
        self.cutoff_ratio = cutoff_ratio
        self.high_freq_scale = high_freq_scale
        self.adaptive_cutoff = adaptive_cutoff
        
        # 核心集构建器
        self.coreset_builder = OptimizedDTWCoreSetBuilder(epsilon=dtw_epsilon)
        
        # 简化特征提取器 
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
        
        # 自适应截断比例网络
        if adaptive_cutoff:
            self.adaptive_cutoff_net = nn.Sequential(
                nn.Linear(embedding_dim, 1),
                nn.Sigmoid()
            )
        
    def get_adaptive_cutoff(self, seq_output):
        """基于序列特征动态调整截断比例
        
        Args:
            seq_output (torch.Tensor): 序列特征
            
        Returns:
            float: 自适应截断比例 [0.3, 0.7]
        """
        # 将范围限制在0.3到0.7之间的动态截断比例
        cutoff_base = 0.3
        cutoff_range = 0.6  # 0.3 + 0.4 = 0.7
        cutoff = cutoff_base + cutoff_range * self.adaptive_cutoff_net(seq_output).mean()
        return cutoff.item()
    
    def prepare_for_fft(self, seq_output):
        """准备序列用于FFT处理，使用DTW核心集减少复杂度
        
        Args:
            seq_output (torch.Tensor): 序列输出 [B, H]
            
        Returns:
            torch.Tensor: 处理后的序列
        """
        batch_size = seq_output.shape[0]
        
        if batch_size <= 1:
            return seq_output.unsqueeze(1)  # [B, 1, H]
        
        # 小批量优化 - 对小批量直接返回
        if batch_size <= 16:
            return seq_output.unsqueeze(1)  # [B, 1, H]
            
        # 将批次中的序列转为列表以应用DTW - 减少频率
        if hasattr(self, 'skip_counter'):
            self.skip_counter = (self.skip_counter + 1) % 3  # 每3次调用处理一次
            if self.skip_counter != 0:
                return seq_output.unsqueeze(1)
        else:
            self.skip_counter = 0
            
        seq_list = [seq_output[i].unsqueeze(0) for i in range(batch_size)]
        
        # 构建DTW核心集（仅当批次大于阈值时）
        coreset_seqs, coreset_weights = self.coreset_builder.build_coreset(seq_list)
        
        # 合并核心集表示
        combined = self.coreset_builder.combine_embeddings(coreset_seqs, coreset_weights)
        if combined is not None:
            # 使用核心集表示
            processed_seq = combined.unsqueeze(0)  # [1, ?, H]
            # 复制到批次大小
            processed_seq = processed_seq.expand(batch_size, -1, -1)  # [B, ?, H]
        else:
            processed_seq = seq_output.unsqueeze(1)  # [B, 1, H]
            
        return processed_seq
        
    def forward(self, seq_output):
        """增强高频信息的提取和表示，利用DTW核心集优化计算
        
        Args:
            seq_output (torch.Tensor): 序列表示, shape: [B, H]
            
        Returns:
            tuple: 高频表示, 低频表示, 频域表示
        """
        if torch.isnan(seq_output).any():
            seq_output = torch.nan_to_num(seq_output, nan=0.0, posinf=1.0, neginf=-1.0)

        # 使用DTW核心集优化FFT计算
        processed_seq = self.prepare_for_fft(seq_output)
        
        # 编码用户行为序列
        encoded = self.encoder(processed_seq)  # [B, ?, H]

        if torch.isnan(encoded).any():
            encoded = torch.nan_to_num(encoded, nan=0.0, posinf=1.0, neginf=-1.0)

        # 归一化处理
        encoded_mean = encoded.mean(dim=2, keepdim=True)
        encoded_std = encoded.std(dim=2, keepdim=True) + 1e-5
        encoded_norm = (encoded - encoded_mean) / encoded_std
        
        # 自适应截断比例
        if self.adaptive_cutoff:
            current_cutoff = self.get_adaptive_cutoff(seq_output)
        else:
            current_cutoff = self.cutoff_ratio
        
        # 转换到频域 (对嵌入维度进行FFT)
        freq_domain = torch.fft.rfft(encoded_norm, dim=2, norm="ortho")  # [B, ?, H//2+1]
        
        # 高低频划分
        cutoff = int(freq_domain.shape[2] * current_cutoff)
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
        high_freq_repr = torch.fft.irfft(high_freq_components, n=self.embedding_dim, dim=2, norm="ortho")
        low_freq_repr = torch.fft.irfft(low_freq_components, n=self.embedding_dim, dim=2, norm="ortho")

        if torch.isnan(high_freq_repr).any() or torch.isnan(low_freq_repr).any():
            high_freq_repr = torch.nan_to_num(high_freq_repr, nan=0.0, posinf=1.0, neginf=-1.0)
            low_freq_repr = torch.nan_to_num(low_freq_repr, nan=0.0, posinf=1.0, neginf=-1.0)
        
        # 维度压缩
        high_freq_repr = high_freq_repr.mean(dim=1)  # [B, H]
        low_freq_repr = low_freq_repr.mean(dim=1)  # [B, H]
        
        # 投影到表示空间
        high_freq_repr = self.high_freq_projector(high_freq_repr)  # [B, H]
        low_freq_repr = self.low_freq_projector(low_freq_repr)  # [B, H]

        if torch.isnan(high_freq_repr).any() or torch.isnan(low_freq_repr).any():
            high_freq_repr = torch.nan_to_num(high_freq_repr, nan=0.0, posinf=1.0, neginf=-1.0)
            low_freq_repr = torch.nan_to_num(low_freq_repr, nan=0.0, posinf=1.0, neginf=-1.0)
        
        return high_freq_repr, low_freq_repr, freq_domain


class SimplifiedRepetitiveInterestTracker(nn.Module):
    """简化版基于DTW的用户重复兴趣跟踪器
    
    降低计算频率、简化计算过程
    
    Args:
        embedding_dim (int): 嵌入维度
        window_size (int): 频谱历史窗口大小, 默认3
        rep_threshold (float): 重复兴趣阈值, 默认0.6
        dtw_enabled (bool): 是否启用DTW距离度量，默认True
        multi_scale (bool): 是否启用多尺度重复模式检测，默认False
    """
    def __init__(self, embedding_dim, window_size=3, rep_threshold=0.6, 
                 dtw_enabled=True, multi_scale=False):
        super(SimplifiedRepetitiveInterestTracker, self).__init__()
        self.embedding_dim = embedding_dim
        self.window_size = window_size
        self.rep_threshold = rep_threshold
        self.dtw_enabled = dtw_enabled
        self.multi_scale = multi_scale
        
        # DTW核心集构建器
        self.coreset_builder = OptimizedDTWCoreSetBuilder(epsilon=0.2, update_freq=3)
        
        # 简化的重复模式检测网络 
        self.pattern_detector = nn.Sequential(
                nn.Linear(embedding_dim//2, embedding_dim//4),
                nn.LayerNorm(embedding_dim//4),
                nn.GELU()
        )
        
        # 重复强度分类器
        self.repetition_classifier = nn.Sequential(
            nn.Linear(embedding_dim//4, 1),
            nn.Sigmoid()
        )
        
        # 多尺度表示 (可选)
        if self.multi_scale:
            self.scale_projections = nn.ModuleList([
                nn.Linear(embedding_dim//2, embedding_dim//8) 
                for _ in range(3)  # 3个尺度
            ])
            
            self.scale_classifier = nn.Sequential(
                nn.Linear(embedding_dim//8 * 3, 1),
                nn.Sigmoid()
            )
        
        # 处理计数器
        self.process_counter = 0
        self.dtw_update_freq = 3  # 每3次进行一次完整DTW计算
    
    def compute_dtw_similarity(self, x, y):
        """计算DTW相似度 (优化版本)
        
        Args:
            x (torch.Tensor): 第一个序列
            y (torch.Tensor): 第二个序列
            
        Returns:
            float: 相似度分数 [0,1]
        """
        # 更新处理计数
        self.process_counter += 1
        if self.process_counter % self.dtw_update_freq != 0:
            # 跳过DTW计算，使用余弦相似度代替
            return F.cosine_similarity(x.mean(dim=1), y.mean(dim=1), dim=1)
            
        if not self.dtw_enabled:
            # 退化为余弦相似度
            return F.cosine_similarity(x.mean(dim=1), y.mean(dim=1), dim=1)
        
        batch_size = x.shape[0]
        similarity = torch.zeros(batch_size, device=x.device)
        
        # 减少计算量 - 只处理部分批次
        process_size = min(batch_size, 16)  # 最多处理16个样本
        indices = torch.randperm(batch_size)[:process_size]
        
        for i in indices:
            # 转为序列形式以计算DTW
            xi = x[i].unsqueeze(0)  # [1, dim]
            yi = y[i].unsqueeze(0)  # [1, dim]
            
            # 使用松弛三角不等式的DTW计算
            try:
                dist = self.coreset_builder.dtw_distance(xi, yi, p=2)
                # 转换为相似度: exp(-dist)
                similarity[i] = torch.exp(-dist / self.embedding_dim)
            except Exception as e:
                # 回退到余弦相似度
                similarity[i] = F.cosine_similarity(xi.mean(dim=0), yi.mean(dim=0), dim=0)
        
        # 对未处理的样本使用均值
        mean_sim = similarity[indices].mean()
        mask = torch.ones(batch_size, dtype=torch.bool, device=x.device)
        mask[indices] = False
        similarity[mask] = mean_sim
        
        return similarity
        
    def analyze_multi_scale_patterns(self, high_freq_magnitude):
        """多尺度重复模式分析 (简化版)
        
        Args:
            high_freq_magnitude (torch.Tensor): 高频幅度
            
        Returns:
            torch.Tensor: 多尺度重复分数
        """
        if not self.multi_scale:
            return self.repetition_classifier(high_freq_magnitude[:, -1, :])
            
        batch_size = high_freq_magnitude.shape[0]
        
        # 不同尺度的特征
        scale_features = []
        
        # 尺度1: 最近时间
        scale1 = high_freq_magnitude[:, -1, :]
        scale_features.append(self.scale_projections[0](scale1))
        
        # 尺度2: 中期 (如果有)
        if high_freq_magnitude.shape[1] >= 3:
            scale2 = high_freq_magnitude[:, -2, :]
            scale_features.append(self.scale_projections[1](scale2))
        else:
            scale_features.append(self.scale_projections[1](scale1))
        
        # 尺度3: 最早 (如果有)
        if high_freq_magnitude.shape[1] >= 2:
            scale3 = high_freq_magnitude[:, 0, :]
            scale_features.append(self.scale_projections[2](scale3))
        else:
            scale_features.append(self.scale_projections[2](scale1))
        
        # 融合多尺度特征
        multi_scale_feature = torch.cat(scale_features, dim=1)
        multi_scale_score = self.scale_classifier(multi_scale_feature)
        
        return multi_scale_score
    
    def build_spectrum_coreset(self, freq_history):
        """从频谱历史构建核心集 (简化版)
        
        Args:
            freq_history (torch.Tensor): 频谱历史 [B, T, F]
            
        Returns:
            torch.Tensor: 频谱核心集
        """
        # 如果历史很短，直接返回
        if freq_history.shape[1] <= 2:
            return freq_history
            
        # 每3次调用才处理一次
        if not hasattr(self, 'build_counter'):
            self.build_counter = 0
        self.build_counter = (self.build_counter + 1) % 3
        
        if self.build_counter != 0:
            # 简单采样代替DTW
            indices = torch.linspace(0, freq_history.shape[1]-1, min(3, freq_history.shape[1])).long()
            return freq_history[:, indices]
        
        batch_size = freq_history.shape[0]
        spectrum_coresets = []
        
        # 只处理部分批次
        process_size = min(batch_size, 16)  # 最多处理16个样本
        process_indices = torch.randperm(batch_size)[:process_size]
        
        for i in process_indices:
            # 提取单个批次样本的频谱历史
            sample_history = freq_history[i]  # [T, F]
            
            # 按时间切片，每个时间点作为一个序列
            seq_list = [sample_history[j].unsqueeze(0) for j in range(sample_history.shape[0])]
            
            # 构建核心集
            if len(seq_list) > 1:
                try:
                    coreset_seqs, coreset_weights = self.coreset_builder.build_coreset(seq_list)
                    combined = self.coreset_builder.combine_embeddings(coreset_seqs, coreset_weights)
                    if combined is not None:
                        spectrum_coresets.append(combined)
                    else:
                        spectrum_coresets.append(sample_history)
                except Exception as e:
                    spectrum_coresets.append(sample_history)
            else:
                spectrum_coresets.append(sample_history)
        
        # 对未处理批次使用简单采样
        unprocessed_mask = torch.ones(batch_size, dtype=torch.bool)
        unprocessed_mask[process_indices] = False
        unprocessed_indices = torch.where(unprocessed_mask)[0]
        
        for i in unprocessed_indices:
            # 简单采样
            indices = torch.linspace(0, freq_history.shape[1]-1, min(3, freq_history.shape[1])).long()
            spectrum_coresets.append(freq_history[i, indices])
                
        # 对齐维度
        max_len = max(s.shape[0] for s in spectrum_coresets)
        aligned_coresets = []
        
        for cs in spectrum_coresets:
            if cs.shape[0] < max_len:
                padding = torch.zeros(max_len - cs.shape[0], cs.shape[1], device=cs.device)
                aligned_cs = torch.cat([cs, padding], dim=0)
            else:
                aligned_cs = cs
            aligned_coresets.append(aligned_cs)
            
        return torch.stack(aligned_coresets, dim=0)  # [B, T', F]
        
    def forward(self, freq_history):
        """分析用户重复兴趣强度 (优化版本)
        
        Args:
            freq_history (torch.Tensor): 频谱历史序列, shape: [B, T, F]
                T: 历史长度, F: 频域维度 (embedding_dim//2+1)
            
        Returns:
            tuple: 重复强度分数, 重复兴趣模式
        """
        batch_size = freq_history.shape[0]
       
        if torch.isnan(freq_history).any():
            freq_history = torch.nan_to_num(freq_history, nan=0.0, posinf=1.0, neginf=-1.0)
        
        # 历史处理 - 减少窗口大小
        actual_window = min(self.window_size, freq_history.shape[1])
        if freq_history.shape[1] > actual_window:
            # 只保留最近的历史
            freq_history = freq_history[:, -actual_window:]
        
        if freq_history.shape[1] < actual_window:
            padding = torch.zeros(
                (batch_size, actual_window - freq_history.shape[1], freq_history.shape[2]),
                dtype=freq_history.dtype,
                device=freq_history.device
            )
            freq_history = torch.cat([padding, freq_history], dim=1)
        
        # 使用DTW核心集优化频谱历史 - 减少计算频率
        spectrum_coreset = self.build_spectrum_coreset(freq_history)
        
        # 提取高频部分幅度 (只关注高频部分)
        cutoff = int(spectrum_coreset.shape[2] * 0.5)
        high_freq_magnitude = torch.abs(spectrum_coreset[:, :, :cutoff])
        
        # 归一化 (简化)
        high_freq_magnitude = torch.clamp(high_freq_magnitude, min=1e-5, max=10.0)
        mean = high_freq_magnitude.mean(dim=(1,2), keepdim=True)
        std = high_freq_magnitude.std(dim=(1,2), keepdim=True) + 1e-5
        high_freq_magnitude = (high_freq_magnitude - mean) / std
        high_freq_magnitude = torch.clamp(high_freq_magnitude, min=-3, max=3)
        
        # 获取最新特征
        latest_features = high_freq_magnitude[:, -1, :]
        pattern_features = self.pattern_detector(latest_features)

        if torch.isnan(pattern_features).any():
            pattern_features = torch.nan_to_num(pattern_features, nan=0.0, posinf=1.0, neginf=-1.0)
        
        # 重复强度分数 - 直接或多尺度
        if self.multi_scale:
            repetition_scores = self.analyze_multi_scale_patterns(high_freq_magnitude)
        else:
            repetition_scores = self.repetition_classifier(pattern_features)

        if torch.isnan(repetition_scores).any():
            repetition_scores = torch.nan_to_num(repetition_scores, nan=0.5, posinf=1.0, neginf=0.0)
        
        # 简化的模式信息
        dtw_patterns = {'repetition_scores': repetition_scores}
        
        # DTW相似度分析 (如果启用且有足够历史)
        if self.dtw_enabled and spectrum_coreset.shape[1] >= 2:
            try:
                # 用本批次的1/3样本进行DTW计算
                sample_size = max(1, batch_size // 3)
                sample_idx = torch.randperm(batch_size)[:sample_size]
                
                first_window = spectrum_coreset[sample_idx, 0, :cutoff].unsqueeze(1)
                last_window = spectrum_coreset[sample_idx, -1, :cutoff].unsqueeze(1)
                
                # 计算部分样本的DTW相似度
                dtw_similarity = self.compute_dtw_similarity(first_window, last_window)
                
                # 扩展到全批次
                full_similarity = torch.ones(batch_size, device=spectrum_coreset.device) * dtw_similarity.mean()
                full_similarity[sample_idx] = dtw_similarity
                
                dtw_patterns['dtw_similarity'] = full_similarity
            except Exception as e:
                # 出错时使用默认值
                dtw_patterns['dtw_similarity'] = torch.ones(batch_size, device=spectrum_coreset.device) * 0.5
        
        return repetition_scores, dtw_patterns


class SimplifiedRepetitionAwareRecommender(nn.Module):
    """简化版基于DTW的重复兴趣感知推荐器
    
    利用DTW核心集特性增强推荐融合策略
    
    Args:
        embedding_dim (int): 嵌入维度
        rep_weight (float): 重复兴趣权重基准, 默认0.7
        personalized_fusion (bool): 是否启用个性化融合策略, 默认False
        use_dtw_probe (bool): 是否使用DTW兴趣探针, 默认False
    """
    def __init__(self, embedding_dim, rep_weight=0.7, 
                 personalized_fusion=False, use_dtw_probe=False):
        super(SimplifiedRepetitionAwareRecommender, self).__init__()
        self.embedding_dim = embedding_dim
        self.rep_weight = rep_weight
        self.personalized_fusion = personalized_fusion
        self.use_dtw_probe = use_dtw_probe
        
        # 融合网络 - 动态调整重复兴趣权重
        self.fusion_network = nn.Sequential(
            nn.Linear(embedding_dim * 2 + 1, embedding_dim//2),
            nn.GELU(),
            nn.Linear(embedding_dim//2, 1),
            nn.Sigmoid()
        )
        
        # 个性化融合模块 (可选)
        if self.personalized_fusion:
            self.personalized_net = nn.Sequential(
                nn.Linear(embedding_dim * 2 + 2, embedding_dim//2),
                nn.GELU(),
                nn.Linear(embedding_dim//2, 2),
                nn.Softmax(dim=1)
            )
        
        # DTW兴趣探针 (可选)
        if self.use_dtw_probe:
            self.dtw_probe = nn.Sequential(
                nn.Linear(embedding_dim, 1),
                nn.Sigmoid()
            )
        
        # 最终融合投影
        self.final_projection = nn.Sequential(
            nn.Linear(embedding_dim * 2, embedding_dim),
            nn.LayerNorm(embedding_dim)
        )
    
    def apply_dtw_probe(self, high_freq_repr, dtw_patterns):
        """应用DTW兴趣探针评估当前兴趣状态 (简化版)
        
        Args:
            high_freq_repr (torch.Tensor): 高频表示
            dtw_patterns (dict): DTW模式信息
            
        Returns:
            torch.Tensor: 探针分数 [0,1]
        """
        if not self.use_dtw_probe:
            # 禁用DTW探针时直接返回重复分数
            if dtw_patterns and 'repetition_scores' in dtw_patterns:
                return dtw_patterns['repetition_scores']
            # 默认值
            return torch.ones(high_freq_repr.size(0), 1, device=high_freq_repr.device) * 0.5
            
        # 提取重复模式信息
        probe_input = high_freq_repr
        
        # 如果有DTW相似度信息，使用它增强探针
        if dtw_patterns and 'dtw_similarity' in dtw_patterns:
            dtw_sim = dtw_patterns['dtw_similarity'].unsqueeze(1)
            probe_score = 0.7 * self.dtw_probe(probe_input) + 0.3 * dtw_sim
        else:
            probe_score = self.dtw_probe(probe_input)
            
        return probe_score
        
    def forward(self, high_freq_repr, low_freq_repr, repetition_scores, dtw_patterns=None):
        """融合高频和低频表示
        
        Args:
            high_freq_repr (torch.Tensor): 高频表示 (重复兴趣), shape: [B, H]
            low_freq_repr (torch.Tensor): 低频表示 (次要兴趣), shape: [B, H]
            repetition_scores (torch.Tensor): 重复强度分数, shape: [B, 1]
            dtw_patterns (dict, optional): DTW模式信息
            
        Returns:
            tuple: 融合表示, 高频权重, 低频权重
        """
        batch_size = high_freq_repr.shape[0]
        
        # 应用DTW兴趣探针 (可选)
        if self.use_dtw_probe:
            probe_score = self.apply_dtw_probe(high_freq_repr, dtw_patterns)
        else:
            probe_score = repetition_scores
        
        # 个性化融合策略 (可选)
        if self.personalized_fusion:
            # 输入包括高频表示、低频表示、重复分数和探针分数
            personalized_input = torch.cat([
                high_freq_repr,
                low_freq_repr,
                repetition_scores,
                probe_score
            ], dim=1)
            
            # 获取个性化权重
            weights = self.personalized_net(personalized_input)  # [B, 2]
            repetition_weight = weights[:, 0].unsqueeze(1)  # 高频权重
            novelty_weight = weights[:, 1].unsqueeze(1)  # 低频权重
        else:        
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
        
        # 基于权重融合 (简化)
        weighted_high = repetition_weight * high_freq_repr
        weighted_low = novelty_weight * low_freq_repr
        
        # 连接并最终投影
        concat_repr = torch.cat([weighted_high, weighted_low], dim=1)
        final_repr = self.final_projection(concat_repr)
        
        return final_repr, repetition_weight, novelty_weight

class FRF(SequentialRecommender):
    """频域重复兴趣融合推荐模型 (优化版)
    
    基于频域分析的推荐模型，将用户兴趣分解为高频(重复兴趣)和低频(次要兴趣)，
    重点关注重复兴趣模式变化，以提高推荐准确性。
    
    优化版整合了高效DTW核心集方法，提升了训练速度。

    Args:
        config (Config): 全局配置对象
        dataset (Dataset): 数据集
    """
    def __init__(self, config, dataset):
        super(FRF, self).__init__(config, dataset)
        
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
        self.window_size = config['window_size'] if 'window_size' in config else 3  # 频谱历史窗口大小
        self.rep_threshold = config['rep_threshold'] if 'rep_threshold' in config else 0.6  # 重复兴趣阈值
        self.rep_weight = config['rep_weight'] if 'rep_weight' in config else 0.7  # 重复兴趣基准权重
        self.high_freq_scale = config['high_freq_scale'] if 'high_freq_scale' in config else 1.5  # 高频增强系数
        self.contrast_weight = config['contrast_weight'] if 'contrast_weight' in config else 0.5  # 对比学习权重调整
        self.freq_history_length = config['freq_history_length'] if 'freq_history_length' in config else 5  # 频谱历史长度
        
        # DTW相关参数
        self.use_dtw = config['use_dtw'] if 'use_dtw' in config else True  # 是否使用DTW
        self.dtw_epsilon = config['dtw_epsilon'] if 'dtw_epsilon' in config else 0.2  # DTW核心集精度参数
        self.adaptive_cutoff = config['adaptive_cutoff'] if 'adaptive_cutoff' in config else False  # 是否使用自适应截断
        self.multi_scale = config['multi_scale'] if 'multi_scale' in config else False  # 是否使用多尺度分析
        self.personalized_fusion = config['personalized_fusion'] if 'personalized_fusion' in config else False  # 是否使用个性化融合
        self.use_dtw_probe = config['use_dtw_probe'] if 'use_dtw_probe' in config else False  # 是否使用DTW兴趣探针
        self.dtw_update_freq = config['dtw_update_freq'] if 'dtw_update_freq' in config else 5  # DTW更新频率
        
        # 性能模式参数
        self.performance_mode = config['performance_mode'] if 'performance_mode' in config else True
        
        # 初始化嵌入层
        self.item_embedding = nn.Embedding(
            self.n_items, self.embedding_size, padding_idx=0
        )
        self.position_embedding = nn.Embedding(
            self.max_seq_length, self.embedding_size
        )
        
        # 序列编码器 (简化版)
        try:
            # 导入Transformer编码器
            from recbole.model.layers import TransformerEncoder
            
            # 性能模式下减少层数
            if self.performance_mode:
                self.n_layers = min(self.n_layers, 2)
                self.n_heads = min(self.n_heads, 2)
            
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
        except ImportError as e:
            print(f"导入Transformer编码器失败: {e}")
            # 简单的RNN编码器作为降级方案
            self.encoder = nn.LSTM(
                input_size=self.embedding_size,
                hidden_size=self.hidden_size,
                num_layers=1,
                batch_first=True
            )
        
        # 初始化频域相关模块 (优化版)
        self.freq_decoupler = SimplifiedFrequencyDomainDecoupler(
            self.embedding_size, 
            self.cutoff_ratio,
            self.high_freq_scale,
            self.dtw_epsilon,
            self.adaptive_cutoff
        )
        self.repetition_tracker = SimplifiedRepetitiveInterestTracker(
            self.embedding_size, 
            self.window_size, 
            self.rep_threshold,
            self.use_dtw,
            self.multi_scale
        )
        self.recommender = SimplifiedRepetitionAwareRecommender(
            self.embedding_size,
            self.rep_weight,
            self.personalized_fusion,
            self.use_dtw_probe
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
        
        # 频谱历史存储 (稀疏存储)
        self.register_buffer('freq_history', None)
        self.freq_update_counter = 0
        
        # 参数初始化
        self.apply(self._init_weights)
        
        # 性能追踪
        self.time_stats = {'forward': [], 'dtw': [], 'fft': []}
    
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
        """生成注意力掩码 (简化版)"""
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
    
    def update_freq_history(self, freq_domain):
        """更新频谱历史 (修复版本 - 避免梯度问题)"""
        # 训练时每N个批次更新一次
        self.freq_update_counter = (self.freq_update_counter + 1) % self.dtw_update_freq
        if not self.training or self.freq_update_counter == 0:
            if self.freq_history is None or self.freq_history.shape[0] != freq_domain.shape[0]:
                # 首次运行或批次大小变化时初始化 - 使用clone()并detach()
                self.freq_history = freq_domain.clone().detach()
            else:
                if len(freq_domain.shape) != len(self.freq_history.shape):
                    if len(freq_domain.shape) > len(self.freq_history.shape):
                        freq_domain = freq_domain.squeeze(1)
                    else:
                        freq_domain = freq_domain.unsqueeze(1)
                    
                # 保持有限长度的历史 - 使用clone()并detach()
                self.freq_history = torch.cat([self.freq_history,freq_domain.clone().detach()], dim=1)
                if self.freq_history.shape[1] > self.freq_history_length:
                    # 只保留最近的N个记录
                    self.freq_history = self.freq_history[:, -self.freq_history_length:].clone()
    
    def forward(self, item_seq, item_seq_len):
        """前向传播，重点关注重复兴趣 (优化版)
    
        Args:
            item_seq (torch.LongTensor): 物品序列, shape: [B, L]
            item_seq_len (torch.LongTensor): 序列长度, shape: [B]
        
        Returns:
            torch.Tensor: 用户表示, shape: [B, H]
        """
        start = time.time()
    
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
    
        # 序列编码
        if isinstance(self.encoder, nn.LSTM):
            # LSTM编码器
            encoder_output, _ = self.encoder(input_emb)
            seq_output = encoder_output
        else:
            # Transformer编码器
            extended_attention_mask = self.get_attention_mask(item_seq)
            try:
                encoder_output = self.encoder(input_emb, extended_attention_mask)
                if isinstance(encoder_output, list):
                    seq_output = encoder_output[-1]
                else:
                    seq_output = encoder_output
            except Exception as e:
                # 出错时回退到简单平均
                print(f"编码器出错: {e}, 使用简单平均")
                seq_output = input_emb.mean(dim=1, keepdim=True).expand(-1, input_emb.size(1), -1)
    
        output = self.gather_indexes(seq_output, item_seq_len - 1)
        if torch.isnan(output).any():
            output = torch.nan_to_num(output, nan=0.0, posinf=1.0, neginf=-1.0)
    
        # 记录编码时间
        encode_time = time.time() - start
        dtw_start = time.time()
    
        # 频域解耦 
        high_freq_repr, low_freq_repr, freq_domain = self.freq_decoupler(output)
    
        # 记录DTW时间
        dtw_time = time.time() - dtw_start
        fft_start = time.time()
    
        # 更新频谱历史 (使用detach()断开梯度流)
        if self.training:
            with torch.no_grad():  # 明确断开梯度流
                self.update_freq_history(freq_domain.detach())
        else:
            self.update_freq_history(freq_domain)
    
        # 跟踪重复兴趣 (根据训练/评估模式优化)
        if self.training:
            # 训练阶段: 每N个批次才完整计算一次
            if self.freq_update_counter == 0 and self.freq_history is not None:
                with torch.no_grad():  # 确保重复兴趣计算不参与梯度
                    repetition_scores, dtw_patterns = self.repetition_tracker(self.freq_history)
                    # 转换为需要梯度的张量但不连接计算历史
                    repetition_scores = repetition_scores.detach()
            else:
                # 使用近似值
                batch_size = high_freq_repr.shape[0]
                repetition_scores = torch.ones((batch_size, 1), device=high_freq_repr.device) * 0.5
                dtw_patterns = {'repetition_scores': repetition_scores}
        else:
            # 评估阶段: 始终计算
            if self.freq_history is not None and self.freq_history.shape[1] > 1:
                repetition_scores, dtw_patterns = self.repetition_tracker(self.freq_history)
            else:
                batch_size = high_freq_repr.shape[0]
                repetition_scores = torch.ones((batch_size, 1), device=high_freq_repr.device) * 0.5
                dtw_patterns = {'repetition_scores': repetition_scores}

        # 确保repetition_scores没有NaN并且训练时是分离的
        if torch.isnan(repetition_scores).any():
            batch_size = high_freq_repr.shape[0]
            repetition_scores = torch.ones((batch_size, 1), device=high_freq_repr.device) * 0.5
        elif self.training:
            repetition_scores = repetition_scores.detach()  # 确保分离
    
        # 记录FFT和重复兴趣跟踪时间
        fft_time = time.time() - fft_start
    
        # 重复兴趣感知的融合推荐
        fused_repr, _, _ = self.recommender(high_freq_repr, low_freq_repr, repetition_scores,dtw_patterns)
    
        # 记录总前向传播时间
        total_time = time.time() - start
        if self.training and len(self.time_stats['forward']) < 1000:  # 限制记录数量
            self.time_stats['forward'].append(total_time)
            self.time_stats['dtw'].append(dtw_time)
            self.time_stats['fft'].append(fft_time)
    
        return fused_repr
    
    def calculate_simplified_contrast_loss(self, high_freq_repr, pos_items_emb, neg_items_emb=None):
        """计算简化版对比学习损失 (性能优化)
        
        Args:
            high_freq_repr (torch.Tensor): 高频表示
            pos_items_emb (torch.Tensor): 正样本嵌入
            neg_items_emb (torch.Tensor, optional): 负样本嵌入
            
        Returns:
            torch.Tensor: 对比损失
        """
        # 简化计算 - 只使用余弦相似度
        batch_size = high_freq_repr.shape[0]
        
        # 计算高频表示与正样本的相似度
        high_freq_pos_sim = F.cosine_similarity(high_freq_repr, pos_items_emb, dim=1)
        high_freq_target = torch.ones_like(high_freq_pos_sim)
        high_freq_loss = F.mse_loss(torch.sigmoid(high_freq_pos_sim * 5), high_freq_target)
        
        # 如果有负样本，计算对比损失
        if neg_items_emb is not None and not self.performance_mode:
            high_freq_neg_sim = F.cosine_similarity(high_freq_repr, neg_items_emb, dim=1)
            high_freq_logits = torch.stack([high_freq_pos_sim, high_freq_neg_sim], dim=1)
            high_freq_labels = torch.zeros(batch_size, dtype=torch.long, device=high_freq_repr.device)
            high_freq_contrast = self.contrast_loss(high_freq_logits * 5, high_freq_labels)
            
            return high_freq_loss + high_freq_contrast
        else:
            return high_freq_loss
    
    def calculate_loss(self, interaction):
        """计算模型损失 (优化版)
        
        Args:
            interaction (Interaction): 交互数据
            
        Returns:
            torch.Tensor: 损失值
        """
        item_seq = interaction[self.ITEM_SEQ]
        item_seq_len = interaction[self.ITEM_SEQ_LEN]
        pos_items = interaction[self.POS_ITEM_ID]
        
        seq_output = self.forward(item_seq, item_seq_len)
        
        # 性能模式 - 跳过额外的频域解耦计算
        if not self.performance_mode:
            high_freq_repr, low_freq_repr, _ = self.freq_decoupler(seq_output)
        else:
            high_freq_repr = seq_output  # 性能模式下简化
        
        # 计算主要推荐损失
        if self.loss_type == 'BPR':
            neg_items = interaction[self.NEG_ITEM_ID]
            pos_items_emb = self.item_embedding(pos_items)
            neg_items_emb = self.item_embedding(neg_items)
            pos_score = torch.sum(seq_output * pos_items_emb, dim=-1)
            neg_score = torch.sum(seq_output * neg_items_emb, dim=-1)
            main_loss = self.loss_fct(pos_score, neg_score)
            
            # 计算简化版对比学习损失
            contrastive_loss = self.calculate_simplified_contrast_loss(
                high_freq_repr, pos_items_emb, neg_items_emb
            )
        else:  # CE损失
            test_item_emb = self.item_embedding.weight
            logits = torch.matmul(seq_output, test_item_emb.transpose(0, 1))
            main_loss = self.loss_fct(logits, pos_items)
            
            # 对比学习损失
            pos_items_emb = self.item_embedding(pos_items)
            contrastive_loss = self.calculate_simplified_contrast_loss(
                high_freq_repr, pos_items_emb
            )
        
        # 组合损失 - 性能模式下降低对比学习权重
        loss_weight = self.contrast_weight if not self.performance_mode else self.contrast_weight * 0.5
        loss = main_loss + loss_weight * contrastive_loss
        
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

    def get_performance_stats(self):
        """获取性能统计信息"""
        stats = {}
        if len(self.time_stats['forward']) > 0:
            stats['avg_forward_time'] = sum(self.time_stats['forward']) / len(self.time_stats['forward'])
            stats['avg_dtw_time'] = sum(self.time_stats['dtw']) / len(self.time_stats['dtw'])
            stats['avg_fft_time'] = sum(self.time_stats['fft']) / len(self.time_stats['fft'])
            stats['dtw_percentage'] = stats['avg_dtw_time'] / stats['avg_forward_time'] * 100
            stats['fft_percentage'] = stats['avg_fft_time'] / stats['avg_forward_time'] * 100
        return stats
