# -*- coding: utf-8 -*-
# @Time    : 2025/5/16


r"""
FRF with DTW
################################################
Reference:

"""

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
import math
import time
import logging
import os
from tqdm import tqdm
from collections import defaultdict
from scipy.sparse import csr_matrix

from recbole.model.abstract_recommender import SequentialRecommender
from recbole.model.loss import BPRLoss
from recbole.data.interaction import Interaction
from recbole.data.dataset import Dataset
from recbole.utils import set_color


class DTWCoreSetProcessor:
    """基于DTW的核心集处理器
    
    对序列数据进行DTW核心集构建，减少数据规模并保持表达能力。
    
    Args:
        dataset (Dataset): 数据集对象
        embedding_size (int): 嵌入维度
        sample_ratio (float): 采样比例，默认0.3
        sensitivity_alpha (float): 敏感性计算中的α参数，默认2.0
        sensitivity_beta (float): 敏感性计算中的β参数，默认4.0
        dtw_epsilon (float): DTW路径度量近似参数，默认0.1
        use_metric_closure (bool): 是否使用路径度量，默认True
        cache_dir (str): 缓存目录，默认'./dtw_coreset_cache'
        force_recompute (bool): 是否强制重新计算，默认False
        dataset_name (str): 数据集名称，用于缓存文件命名
    """
    def __init__(self, dataset, embedding_size, sample_ratio=0.3, 
                 sensitivity_alpha=2.0, sensitivity_beta=4.0, dtw_epsilon=0.1,
                 use_metric_closure=True, cache_dir='./dtw_coreset_cache', 
                 force_recompute=False, dataset_name=None):
        self.logger = logging.getLogger()
        self.dataset = dataset
        self.embedding_size = embedding_size
        self.sample_ratio = sample_ratio
        self.sensitivity_alpha = sensitivity_alpha
        self.sensitivity_beta = sensitivity_beta
        self.dtw_epsilon = dtw_epsilon
        self.use_metric_closure = use_metric_closure
        self.cache_dir = cache_dir
        self.force_recompute = force_recompute
        self.dataset_name = dataset_name or 'default_dataset'
        
        # 创建缓存目录
        if not os.path.exists(cache_dir):
            os.makedirs(cache_dir)
            
        # DTW距离计算缓存
        self.dtw_cache = {}
        self.max_cache_size = 10000
        
        # 数据集信息
        self.uid_field = dataset.uid_field
        self.iid_field = dataset.iid_field
        self.time_field = dataset.time_field if hasattr(dataset, 'time_field') else None
        self.item_num = dataset.item_num
        
        # 初始化嵌入
        self.item_embeddings = None
        
    def _init_item_embeddings(self):
        """初始化物品嵌入
        
        使用Xavier初始化生成随机嵌入或从预训练模型加载
        """
        self.logger.info("初始化物品嵌入用于DTW核心集计算...")
        
        # 检查是否存在预训练嵌入
        pretrained_path = os.path.join(self.cache_dir, f"{self.dataset_name}_pretrained_emb.npy")
        if os.path.exists(pretrained_path) and not self.force_recompute:
            self.logger.info(f"加载预训练物品嵌入: {pretrained_path}")
            self.item_embeddings = torch.tensor(np.load(pretrained_path), dtype=torch.float32)
        else:
            # 随机初始化
            self.logger.info("使用Xavier初始化物品嵌入")
            self.item_embeddings = torch.zeros((self.item_num + 1, self.embedding_size))
            nn.init.xavier_normal_(self.item_embeddings)
    
    def dtw_distance(self, seq1, seq2, p=2):
        """计算两个序列间的DTW距离
        
        Args:
            seq1 (torch.Tensor): 第一个序列
            seq2 (torch.Tensor): 第二个序列
            p (int): 距离范数，默认为2（欧氏距离）
            
        Returns:
            float: DTW距离值
        """
        # 将序列转为numpy进行缓存
        try:
            seq1_bytes = seq1.cpu().numpy().tobytes()
            seq2_bytes = seq2.cpu().numpy().tobytes()
            cache_key = (hash(seq1_bytes), hash(seq2_bytes))
            
            if cache_key in self.dtw_cache:
                return self.dtw_cache[cache_key]
        except:
            cache_key = None
        
        # 序列长度
        n, m = len(seq1), len(seq2)
        
        # 对于非常短的序列，简化为欧式距离
        if n <= 2 or m <= 2:
            if n == m:
                dist = np.linalg.norm(seq1 - seq2) / max(n, 1)
            else:
                # 简单的均值比较
                dist = np.linalg.norm(np.mean(seq1, axis=0) - np.mean(seq2, axis=0))
            
            if cache_key is not None:
                self.dtw_cache[cache_key] = dist
            return dist
        
        # 创建DTW矩阵
        dtw_matrix = np.zeros((n+1, m+1)) + np.inf
        dtw_matrix[0, 0] = 0
        
        # 计算成对距离
        for i in range(1, n+1):
            for j in range(1, m+1):
                cost = np.linalg.norm(seq1[i-1] - seq2[j-1], ord=p)
                dtw_matrix[i, j] = cost + min(
                    dtw_matrix[i-1, j],    # 插入
                    dtw_matrix[i, j-1],    # 删除
                    dtw_matrix[i-1, j-1]   # 替换
                )
        
        dist = dtw_matrix[n, m]
        
        # 存入缓存
        if cache_key is not None:
            if len(self.dtw_cache) >= self.max_cache_size:
                # 清理部分缓存
                keys_to_remove = list(self.dtw_cache.keys())[:self.max_cache_size//2]
                for k in keys_to_remove:
                    self.dtw_cache.pop(k, None)
            
            self.dtw_cache[cache_key] = dist
            
        return dist
    
    def compute_path_metric(self, seq1, seq2, intermediates):
        """计算路径度量（满足三角不等式的DTW变种）
        
        Args:
            seq1 (np.ndarray): 第一个序列
            seq2 (np.ndarray): 第二个序列
            intermediates (list): 中间点序列
            
        Returns:
            float: 路径度量距离
        """
        if not intermediates:
            return self.dtw_distance(seq1, seq2)
        
        # 通过中间点计算路径
        path_dist = self.dtw_distance(seq1, intermediates[0])
        
        for i in range(len(intermediates) - 1):
            path_dist += self.dtw_distance(intermediates[i], intermediates[i+1])
            
        path_dist += self.dtw_distance(intermediates[-1], seq2)
        
        return path_dist
    
    def _get_initial_centers(self, sequences, num_centers=5):
        """初始化聚类中心
        
        使用k-means++思想选择初始聚类中心
        
        Args:
            sequences (list): 序列列表
            num_centers (int): 中心数量
            
        Returns:
            list: 初始聚类中心索引
        """
        self.logger.info(f"使用k-means++策略选择{num_centers}个初始聚类中心...")
        n_samples = len(sequences)
        
        # 随机选择第一个中心
        center_indices = [np.random.randint(0, n_samples)]
        centers = [sequences[center_indices[0]]]
        
        # 选择剩余的中心
        for _ in range(1, num_centers):
            # 计算每个点到最近中心的距离
            min_distances = np.array([
                min([self.dtw_distance(seq, center) for center in centers])
                for seq in sequences
            ])
            
            # 按距离加权随机选择下一个中心
            weights = min_distances / min_distances.sum()
            next_center_idx = np.random.choice(n_samples, p=weights)
            
            center_indices.append(next_center_idx)
            centers.append(sequences[next_center_idx])
        
        return center_indices
    
    def _compute_sensitivities(self, sequences, initial_centers_idx):
        """计算每个序列的敏感性
        
        Args:
            sequences (list): 序列列表
            initial_centers_idx (list): 初始中心索引
            
        Returns:
            np.ndarray: 敏感性数组
        """
        self.logger.info("计算序列敏感性...")
        n_samples = len(sequences)
        sensitivities = np.zeros(n_samples)
        
        # 获取初始聚类中心
        centers = [sequences[idx] for idx in initial_centers_idx]
        
        # 计算Voronoi分区
        assignments = np.zeros(n_samples, dtype=int)
        for i, seq in enumerate(sequences):
            min_dist = float('inf')
            min_idx = 0
            
            for c_idx, center in enumerate(centers):
                dist = self.dtw_distance(seq, center)
                if dist < min_dist:
                    min_dist = dist
                    min_idx = c_idx
            
            assignments[i] = min_idx
        
        # 计算总聚类成本和每个分区的成本
        total_cost = 0
        partition_costs = np.zeros(len(centers))
        partition_sizes = np.zeros(len(centers))
        
        for i, seq in enumerate(sequences):
            center_idx = assignments[i]
            dist = self.dtw_distance(seq, centers[center_idx])
            total_cost += dist
            partition_costs[center_idx] += dist
            partition_sizes[center_idx] += 1
        
        # 使用估计公式计算敏感性上界
        for i, seq in enumerate(sequences):
            center_idx = assignments[i]
            seq_dist = self.dtw_distance(seq, centers[center_idx])
            
            # 防止除零
            v_size = max(1, partition_sizes[center_idx])
            
            # 敏感性上界公式
            term1 = (2 * self.sensitivity_alpha * seq_dist) / total_cost
            term2 = 4 / v_size
            term3 = (8 * self.sensitivity_alpha * partition_costs[center_idx]) / (total_cost * v_size)
            
            # 计算敏感性
            seq_len = len(seq)
            d = seq[0].shape[0] if len(seq) > 0 else self.embedding_size
            sensitivity = (seq_len * d)**(1/2) * (term1 + term2 + term3)
            
            # 确保敏感性为正值且不为零
            sensitivities[i] = max(0.01, sensitivity)
        
        return sensitivities
    
    def _simplify_sequences(self, sequences, max_len=None):
        """简化序列，减少复杂度
        
        Args:
            sequences (list): 序列列表
            max_len (int, optional): 最大长度
            
        Returns:
            list: 简化后的序列
        """
        self.logger.info("简化序列复杂度...")
        
        if max_len is None:
            # 确定合适的最大长度
            lengths = [len(seq) for seq in sequences]
            max_len = min(50, np.percentile(lengths, 90))
        
        simplified = []
        for seq in sequences:
            if len(seq) <= max_len:
                simplified.append(seq)
            else:
                # 等间隔采样
                indices = np.linspace(0, len(seq)-1, max_len, dtype=int)
                simplified.append(seq[indices])
        
        return simplified
    
    def _sample_coreset(self, sequences, sensitivities, coreset_size):
        """基于敏感性采样构建核心集
        
        Args:
            sequences (list): 序列列表
            sensitivities (np.ndarray): 敏感性数组
            coreset_size (int): 核心集大小
            
        Returns:
            tuple: 核心集序列，核心集权重，核心集用户ID
        """
        self.logger.info(f"基于敏感性采样构建大小为{coreset_size}的核心集...")
        n_samples = len(sequences)
        
        # 计算采样概率
        total_sensitivity = np.sum(sensitivities)
        probs = sensitivities / total_sensitivity
        
        # 有放回采样
        sampled_indices = np.random.choice(
            n_samples, size=coreset_size, replace=True, p=probs
        )
        
        # 计算权重
        unique_indices, counts = np.unique(sampled_indices, return_counts=True)
        weights = {}
        
        for idx, count in zip(unique_indices, counts):
            weights[idx] = (count / coreset_size) * (total_sensitivity / sensitivities[idx])
        
        # 获取核心集序列和对应权重
        coreset_sequences = [sequences[idx] for idx in unique_indices]
        coreset_weights = np.array([weights[idx] for idx in unique_indices])
        coreset_uids = [self.uid_list[idx] for idx in unique_indices]
        
        return coreset_sequences, coreset_weights, coreset_uids
    
    def _extract_user_sequences(self):
        """从数据集中提取用户交互序列
        
        Returns:
            tuple: 用户序列列表，用户ID列表
        """
        self.logger.info("从数据集提取用户交互序列...")
        
        # 获取交互数据
        inter_feat = self.dataset.inter_feat
        user_ids = inter_feat[self.uid_field].numpy()
        item_ids = inter_feat[self.iid_field].numpy()
        
        # 时间信息（如果有）
        if self.time_field is not None and self.time_field in inter_feat:
            times = inter_feat[self.time_field].numpy()
            # 按用户ID和时间排序
            sort_idx = np.lexsort((times, user_ids))
        else:
            # 仅按用户ID排序
            sort_idx = np.argsort(user_ids)
        
        user_ids = user_ids[sort_idx]
        item_ids = item_ids[sort_idx]
        
        # 分组构建序列
        user_sequences = defaultdict(list)
        for user_id, item_id in zip(user_ids, item_ids):
            user_sequences[user_id].append(item_id)
        
        # 转换为嵌入序列
        sequences = []
        uid_list = []
        
        for user_id, items in user_sequences.items():
            # 转换物品ID为嵌入
            seq_emb = [self.item_embeddings[item_id].numpy() for item_id in items]
            if len(seq_emb) > 1:  # 忽略长度为1的序列
                sequences.append(np.array(seq_emb))
                uid_list.append(user_id)
        
        self.logger.info(f"提取了{len(sequences)}个用户序列")
        return sequences, uid_list
    
    def _save_coreset_data(self, coreset_uids, coreset_weights):
        """保存核心集数据
        
        Args:
            coreset_uids (list): 核心集用户ID
            coreset_weights (np.ndarray): 核心集权重
        """
        cache_path = os.path.join(self.cache_dir, f"{self.dataset_name}_coreset_info.npz")
        np.savez(
            cache_path,
            uids=np.array(coreset_uids),
            weights=coreset_weights
        )
        self.logger.info(f"核心集数据已保存至: {cache_path}")
    
    def build_coreset(self):
        """构建DTW核心集
        
        返回:
            tuple: 核心集用户ID, 核心集权重
        """
        # 检查缓存
        cache_path = os.path.join(self.cache_dir, f"{self.dataset_name}_coreset_info.npz")
        if os.path.exists(cache_path) and not self.force_recompute:
            self.logger.info(f"加载缓存的核心集数据: {cache_path}")
            data = np.load(cache_path)
            return data['uids'].tolist(), data['weights']
        
        # 初始化物品嵌入
        self._init_item_embeddings()
        
        # 提取用户序列
        sequences, self.uid_list = self._extract_user_sequences()
        
        # 简化序列
        sequences = self._simplify_sequences(sequences)
        
        # 确定核心集大小
        total_users = len(sequences)
        coreset_size = max(int(total_users * self.sample_ratio), 100)
        self.logger.info(f"目标核心集大小: {coreset_size} (原始用户数: {total_users})")
        
        # 初始聚类中心
        k = min(10, max(5, int(np.sqrt(coreset_size))))
        initial_centers = self._get_initial_centers(sequences, num_centers=k)
        
        # 计算敏感性
        sensitivities = self._compute_sensitivities(sequences, initial_centers)
        
        # 采样构建核心集
        _, coreset_weights, coreset_uids = self._sample_coreset(
            sequences, sensitivities, coreset_size
        )
        
        # 保存核心集信息
        self._save_coreset_data(coreset_uids, coreset_weights)
        
        return coreset_uids, coreset_weights
    
    def process_dataset(self):
        """处理数据集，构建核心集并应用到数据集
        
        Returns:
            Dataset: 处理后的数据集
        """
        self.logger.info(set_color("开始DTW核心集处理...", "green"))
        
        # 构建核心集
        coreset_uids, coreset_weights = self.build_coreset()
        
        # 创建核心集用户ID到权重的映射
        uid_to_weight = {uid: weight for uid, weight in zip(coreset_uids, coreset_weights)}
        
        # 应用到数据集
        self.logger.info("基于核心集处理数据集...")
        inter_feat = self.dataset.inter_feat
        uid_field = self.dataset.uid_field
        
        # 获取交互数据中的用户ID
        user_ids = inter_feat[uid_field].numpy()
        
        # 创建用户ID到权重的映射数组
        weights = np.ones(len(user_ids))
        for i, uid in enumerate(user_ids):
            if uid in uid_to_weight:
                weights[i] = uid_to_weight[uid]
            else:
                weights[i] = 0.1  # 非核心集用户给予小权重
        
        # 将权重添加到数据集
        weight_field = 'dtw_coreset_weight'
        inter_feat.update(weight_field, torch.FloatTensor(weights))
        
        # 创建核心集掩码
        coreset_mask = np.zeros(len(user_ids), dtype=bool)
        for i, uid in enumerate(user_ids):
            if uid in uid_to_weight:
                coreset_mask[i] = True
        
        # 将核心集掩码添加到数据集
        mask_field = 'dtw_coreset_mask'
        inter_feat.update(mask_field, torch.BoolTensor(coreset_mask))
        
        self.logger.info(f"数据集处理完成！添加了字段 '{weight_field}' 和 '{mask_field}'")
        
        # 统计信息
        n_coreset = np.sum(coreset_mask)
        self.logger.info(set_color(f"核心集大小: {len(coreset_uids)}, 核心集交互数: {n_coreset}", "green"))
        
        return self.dataset


class DTWFrequencyDomainDecoupler(nn.Module):
    """基于DTW的频域解耦器
    
    将用户兴趣分解为高频(重复兴趣)和低频(次要兴趣)两个独立表示
    
    Args:
        embedding_dim (int): 嵌入维度
        cutoff_ratio (float): 高低频分界点比例, 默认0.5 
        high_freq_scale (float): 高频增强比例, 默认1.5
        adaptive_cutoff (bool): 是否使用自适应截断比例，默认False
    """
    def __init__(self, embedding_dim, cutoff_ratio=0.5, high_freq_scale=1.5, 
                 adaptive_cutoff=False):
        super(DTWFrequencyDomainDecoupler, self).__init__()
        self.embedding_dim = embedding_dim
        self.cutoff_ratio = cutoff_ratio
        self.high_freq_scale = high_freq_scale
        self.adaptive_cutoff = adaptive_cutoff
        
        # 特征提取器
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
        
        # 自适应截断比例网络 (可选)
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
            float: 自适应截断比例 [0.3, 0.9】
        """
        cutoff_base = 0.3
        cutoff_range = 0.55  # 最高到0.95
        cutoff = cutoff_base + cutoff_range * self.adaptive_cutoff_net(seq_output).mean()
        return cutoff.item()
        
    def forward(self, seq_output):
        """增强高频信息的提取和表示
        
        Args:
            seq_output (torch.Tensor): 序列表示, shape: [B, H]
            
        Returns:
            tuple: 高频表示, 低频表示, 频域表示
        """
        if torch.isnan(seq_output).any():
            seq_output = torch.nan_to_num(seq_output, nan=0.0, posinf=1.0, neginf=-1.0)

        # 准备序列用于FFT处理
        processed_seq = seq_output.unsqueeze(1)  # [B, 1, H]
        
        # 编码用户行为序列
        encoded = self.encoder(processed_seq)  # [B, 1, H]

        if torch.isnan(encoded).any():
            encoded = torch.nan_to_num(encoded, nan=0.0, posinf=1.0, neginf=-1.0)

        # 归一化处理
        encoded_mean = encoded.mean(dim=2, keepdim=True)
        encoded_std = encoded.std(dim=2, keepdim=True) + 1e-5
        encoded_norm = (encoded - encoded_mean) / encoded_std
        
        # 获取截断比例
        if self.adaptive_cutoff:
            current_cutoff = self.get_adaptive_cutoff(seq_output)
        else:
            current_cutoff = self.cutoff_ratio
        
        # 转换到频域 (对嵌入维度进行FFT)
        freq_domain = torch.fft.rfft(encoded_norm, dim=2, norm="ortho")  # [B, 1, H//2+1]
        
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


class DTWRepetitiveInterestTracker(nn.Module):
    """基于DTW的用户重复兴趣跟踪器
    
    跟踪用户重复兴趣模式及其强度变化，使用DTW进行序列相似度比较
    
    Args:
        embedding_dim (int): 嵌入维度
        window_size (int): 频谱历史窗口大小, 默认3
        rep_threshold (float): 重复兴趣阈值, 默认0.6
        dtw_enabled (bool): 是否启用DTW距离度量，默认True
        multi_scale (bool): 是否启用多尺度重复模式检测，默认False
        input_dim (int, optional): 高频特征输入维度，默认为32
    """
    def __init__(self, embedding_dim, window_size=3, rep_threshold=0.6, 
                 dtw_enabled=True, multi_scale=False, input_dim=32):
        super(DTWRepetitiveInterestTracker, self).__init__()
        self.embedding_dim = embedding_dim
        self.window_size = window_size
        self.rep_threshold = rep_threshold
        self.dtw_enabled = dtw_enabled
        self.multi_scale = multi_scale
        self.input_dim = input_dim
        
        # DTW 缓存
        self.dtw_cache = {}
        self.max_cache_size = 100
        
        # 目标特征维度 - 使用内部一致的维度
        self.feature_dim = embedding_dim // 4
        
        # 维度适配层 - 使用线性变换代替填充/截断
        if self.input_dim != self.feature_dim:
            self.dim_adapter = nn.Linear(self.input_dim, self.feature_dim)
            # 初始化为近似恒等映射（对于缩小维度）或均匀映射（对于扩大维度）
            if self.input_dim > self.feature_dim:
                # 缩小维度，尝试保留更多信息
                nn.init.xavier_uniform_(self.dim_adapter.weight)
            else:
                # 扩大维度，尽量均匀分布
                nn.init.xavier_uniform_(self.dim_adapter.weight)
            nn.init.zeros_(self.dim_adapter.bias)
        else:
            self.dim_adapter = None
        
        # 重复模式检测网络 
        self.pattern_detector = nn.Sequential(
            nn.Linear(self.feature_dim, self.feature_dim),
            nn.LayerNorm(self.feature_dim),
            nn.GELU()
        )
        
        # 重复强度分类器
        self.repetition_classifier = nn.Sequential(
            nn.Linear(self.feature_dim, 1),
            nn.Sigmoid()
        )
        
        # 多尺度表示 (可选)
        if self.multi_scale:
            self.scale_projections = nn.ModuleList([
                nn.Linear(self.feature_dim, self.feature_dim // 2) 
                for _ in range(3)  # 3个尺度
            ])
            
            self.scale_classifier = nn.Sequential(
                nn.Linear(self.feature_dim // 2 * 3, 1),
                nn.Sigmoid()
            )
    
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
                        
            if cache_key in self.dtw_cache:
                return self.dtw_cache[cache_key]
        except:
            cache_key = None
            
        n, m = seq1.size(0), seq2.size(0)
        
        # 对于非常短的序列，简化为欧式距离
        if n <= 2 or m <= 2:
            if n == m:
                dist = torch.norm(seq1 - seq2, p=p) / max(n, 1)
            else:
                # 简单的均值比较
                dist = torch.norm(seq1.mean(dim=0) - seq2.mean(dim=0), p=p)
            
            if cache_key is not None:
                self.dtw_cache[cache_key] = dist
            return dist
            
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
            if len(self.dtw_cache) >= self.max_cache_size:
                # 清理部分缓存
                keys_to_remove = list(self.dtw_cache.keys())[:self.max_cache_size//2]
                for k in keys_to_remove:
                    self.dtw_cache.pop(k, None)
            
            self.dtw_cache[cache_key] = dist
            
        return dist
    
    def compute_dtw_similarity(self, x, y):
        """计算DTW相似度
        
        Args:
            x (torch.Tensor): 第一个序列
            y (torch.Tensor): 第二个序列
            
        Returns:
            torch.Tensor: 相似度分数 [0,1]
        """
        if not self.dtw_enabled:
            # 退化为余弦相似度
            return F.cosine_similarity(x.mean(dim=1), y.mean(dim=1), dim=1)
        
        batch_size = x.shape[0]
        similarity = torch.zeros(batch_size, device=x.device)
        
        for i in range(batch_size):
            # 转为序列形式以计算DTW
            xi = x[i]
            yi = y[i]
            
            # 使用DTW计算
            try:
                dist = self.dtw_distance(xi, yi, p=2)
                # 转换为相似度: exp(-dist)
                similarity[i] = torch.exp(-dist / self.embedding_dim)
            except Exception as e:
                # 回退到余弦相似度
                similarity[i] = F.cosine_similarity(xi.mean(dim=0, keepdim=True), 
                                                   yi.mean(dim=0, keepdim=True), dim=1)
        
        return similarity
        
    def analyze_multi_scale_patterns(self, high_freq_magnitude):
        """多尺度重复模式分析
        
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
            freq_history = torch.nan_to_num(freq_history, nan=0.0, posinf=1.0, neginf=-1.0)
        
        # 历史处理 - 确保历史长度符合窗口大小要求
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
        
        # 提取高频部分幅度 (只关注高频部分)
        # 动态计算截断点而不是硬编码
        cutoff = int(freq_history.shape[2] * 0.5)
        high_freq_magnitude = torch.abs(freq_history[:, :, :cutoff])
        
        # 归一化
        high_freq_magnitude = torch.clamp(high_freq_magnitude, min=1e-5, max=10.0)
        mean = high_freq_magnitude.mean(dim=(1,2), keepdim=True)
        std = high_freq_magnitude.std(dim=(1,2), keepdim=True) + 1e-5
        high_freq_magnitude = (high_freq_magnitude - mean) / std
        high_freq_magnitude = torch.clamp(high_freq_magnitude, min=-3, max=3)
        
        # 获取最新特征
        latest_features = high_freq_magnitude[:, -1, :]
        
        # 确保输入维度与模型期望匹配
        # 使用线性变换调整维度
        if self.dim_adapter is not None:
            latest_features = self.dim_adapter(latest_features)
        
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
        
        # 模式信息
        dtw_patterns = {'repetition_scores': repetition_scores}
        
        # DTW相似度分析 (如果启用且有足够历史)
        if self.dtw_enabled and freq_history.shape[1] >= 2:
            try:
                first_window = freq_history[:, 0, :cutoff]
                last_window = freq_history[:, -1, :cutoff]
                
                # 计算DTW相似度
                dtw_similarity = self.compute_dtw_similarity(
                    first_window.unsqueeze(1), 
                    last_window.unsqueeze(1)
                )
                dtw_patterns['dtw_similarity'] = dtw_similarity.unsqueeze(1)
            except Exception as e:
                # 出错时使用默认值
                dtw_patterns['dtw_similarity'] = torch.ones(batch_size, 1, device=freq_history.device) * 0.5
        
        return repetition_scores, dtw_patterns


class DTWRepetitionAwareRecommender(nn.Module):
    """基于DTW的重复兴趣感知推荐器
    
    Args:
        embedding_dim (int): 嵌入维度
        rep_weight (float): 重复兴趣权重基准, 默认0.7
        personalized_fusion (bool): 是否启用个性化融合策略, 默认False
        use_dtw_probe (bool): 是否使用DTW兴趣探针, 默认False
    """
    def __init__(self, embedding_dim, rep_weight=0.7, 
                 personalized_fusion=False, use_dtw_probe=False):
        super(DTWRepetitionAwareRecommender, self).__init__()
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
        """应用DTW兴趣探针评估当前兴趣状态
        
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
            dtw_sim = dtw_patterns['dtw_similarity']
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
        
        # 基于权重融合
        weighted_high = repetition_weight * high_freq_repr
        weighted_low = novelty_weight * low_freq_repr
        
        # 连接并最终投影
        concat_repr = torch.cat([weighted_high, weighted_low], dim=1)
        final_repr = self.final_projection(concat_repr)
        
        return final_repr, repetition_weight, novelty_weight


class FRF(SequentialRecommender):
    """频域重复兴趣融合推荐模型 (DTW版)
    
    基于频域分析的推荐模型，将用户兴趣分解为高频(重复兴趣)和低频(次要兴趣)，
    重点关注重复兴趣模式变化，利用DTW分析用户兴趣相似性。

    Args:
        config (Config): 全局配置对象
        dataset (Dataset): 数据集
    """
    def __init__(self, config, dataset):
        super(FRF, self).__init__(config, dataset)
        
        # 获取logger
        self.logger = logging.getLogger()
        
        # 获取数据集名称
        if hasattr(dataset, 'dataset_name'):
            self.dataset_name = dataset.dataset_name
        else:
            self.dataset_name = config['dataset']
        
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
        self.adaptive_cutoff = config['adaptive_cutoff'] if 'adaptive_cutoff' in config else False  # 是否使用自适应截断
        self.multi_scale = config['multi_scale'] if 'multi_scale' in config else False  # 是否使用多尺度分析
        self.personalized_fusion = config['personalized_fusion'] if 'personalized_fusion' in config else False  # 是否使用个性化融合
        self.use_dtw_probe = config['use_dtw_probe'] if 'use_dtw_probe' in config else False  # 是否使用DTW兴趣探针
        self.dtw_update_freq = config['dtw_update_freq'] if 'dtw_update_freq' in config else 5  # DTW更新频率
        
        # 核心集处理相关参数
        self.use_dtw_coreset = config['use_dtw_coreset'] if 'use_dtw_coreset' in config else True  # 是否使用DTW核心集
        self.dtw_coreset_ratio = config['dtw_coreset_ratio'] if 'dtw_coreset_ratio' in config else 0.3  # 核心集采样比例
        self.dtw_coreset_cache = config['dtw_coreset_cache'] if 'dtw_coreset_cache' in config else './dtw_coreset_cache'  # 核心集缓存目录
        self.dtw_force_recompute = config['dtw_force_recompute'] if 'dtw_force_recompute' in config else False  # 强制重新计算核心集
        
        # 检查并处理DTW核心集
        self.has_coreset_weight = False
        if self.use_dtw_coreset:
            self.logger.info(set_color("开始DTW核心集处理...", "green"))
            self.process_dtw_coreset()
        
        # 初始化嵌入层
        self.item_embedding = nn.Embedding(
            self.n_items, self.embedding_size, padding_idx=0
        )
        self.position_embedding = nn.Embedding(
            self.max_seq_length, self.embedding_size
        )
        
        # 序列编码器
        try:
            # 尝试导入Transformer编码器
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
        except ImportError as e:
            self.logger.warning(f"导入Transformer编码器失败: {e}")
            # 简单的RNN编码器作为降级方案
            self.encoder = nn.LSTM(
                input_size=self.embedding_size,
                hidden_size=self.hidden_size,
                num_layers=1,
                batch_first=True
            )
        
        # 初始化频域相关模块
        self.freq_decoupler = DTWFrequencyDomainDecoupler(
            self.embedding_size, 
            self.cutoff_ratio,
            self.high_freq_scale,
            self.adaptive_cutoff
        )
        
        self.repetition_tracker = DTWRepetitiveInterestTracker(
            self.embedding_size, 
            self.window_size, 
            self.rep_threshold,
            self.use_dtw,
            self.multi_scale
        )
        self.recommender = DTWRepetitionAwareRecommender(
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
        
        # 频谱历史存储
        self.register_buffer('freq_history', None)
        self.freq_update_counter = 0
        
        # 参数初始化
        self.apply(self._init_weights)
        
        # 性能追踪
        self.time_stats = {'forward': [], 'dtw': [], 'fft': []}
    
    def process_dtw_coreset(self):
        """处理数据集，应用DTW核心集
        """
        try:
            # 创建DTW核心集处理器
            coreset_processor = DTWCoreSetProcessor(
                dataset=self.dataset,
                embedding_size=self.embedding_size,
                sample_ratio=self.dtw_coreset_ratio,
                cache_dir=self.dtw_coreset_cache,
                force_recompute=self.dtw_force_recompute,
                dataset_name=self.dataset_name
            )
            
            # 处理数据集
            self.dataset = coreset_processor.process_dataset()
            
            # 检查权重字段是否已添加到数据集
            if 'dtw_coreset_weight' in self.dataset.inter_feat:
                self.logger.info(set_color("成功添加DTW核心集权重到数据集", "green"))
                self.has_coreset_weight = True
            else:
                self.logger.warning("未能添加DTW核心集权重到数据集")
                self.has_coreset_weight = False
                
        except Exception as e:
            self.logger.error(f"DTW核心集处理失败: {e}")
            self.has_coreset_weight = False
    
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
    
    def update_freq_history(self, freq_domain):
        """更新频谱历史"""
        # 每N个批次更新一次
        self.freq_update_counter = (self.freq_update_counter + 1) % self.dtw_update_freq
        if not self.training or self.freq_update_counter == 0:
            if self.freq_history is None or self.freq_history.shape[0] != freq_domain.shape[0]:
                # 首次运行或批次大小变化时初始化
                self.freq_history = freq_domain.clone().detach()
            else:
                if len(freq_domain.shape) != len(self.freq_history.shape):
                    if len(freq_domain.shape) > len(self.freq_history.shape):
                        freq_domain = freq_domain.squeeze(1)
                    else:
                        freq_domain = freq_domain.unsqueeze(1)
                    
                # 保持有限长度的历史
                self.freq_history = torch.cat([self.freq_history, freq_domain.clone().detach()], dim=1)
                if self.freq_history.shape[1] > self.freq_history_length:
                    # 只保留最近的N个记录
                    self.freq_history = self.freq_history[:, -self.freq_history_length:].clone()
    
    def forward(self, item_seq, item_seq_len):
        """前向传播，重点关注重复兴趣
    
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
                self.logger.warning(f"编码器出错: {e}, 使用简单平均")
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
    
        # 更新频谱历史
        if self.training:
            with torch.no_grad():  # 明确断开梯度流
                self.update_freq_history(freq_domain.detach())
        else:
            self.update_freq_history(freq_domain)
    
        # 跟踪重复兴趣
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

        # 确保repetition_scores没有NaN
        if torch.isnan(repetition_scores).any():
            batch_size = high_freq_repr.shape[0]
            repetition_scores = torch.ones((batch_size, 1), device=high_freq_repr.device) * 0.5
    
        # 记录FFT和重复兴趣跟踪时间
        fft_time = time.time() - fft_start
    
        # 重复兴趣感知的融合推荐
        fused_repr, _, _ = self.recommender(high_freq_repr, low_freq_repr, repetition_scores, dtw_patterns)
    
        # 记录总前向传播时间
        total_time = time.time() - start
        if self.training and len(self.time_stats['forward']) < 1000:  # 限制记录数量
            self.time_stats['forward'].append(total_time)
            self.time_stats['dtw'].append(dtw_time)
            self.time_stats['fft'].append(fft_time)
    
        return fused_repr
    
    def calculate_simplified_contrast_loss(self, high_freq_repr, pos_items_emb, neg_items_emb=None):
        """计算简化版对比学习损失
        
        Args:
            high_freq_repr (torch.Tensor): 高频表示
            pos_items_emb (torch.Tensor): 正样本嵌入
            neg_items_emb (torch.Tensor, optional): 负样本嵌入
            
        Returns:
            torch.Tensor: 对比损失
        """
        # 计算高频表示与正样本的相似度
        batch_size = high_freq_repr.shape[0]
        high_freq_pos_sim = F.cosine_similarity(high_freq_repr, pos_items_emb, dim=1)
        high_freq_target = torch.ones_like(high_freq_pos_sim)
        high_freq_loss = F.mse_loss(torch.sigmoid(high_freq_pos_sim * 5), high_freq_target)
        
        # 如果有负样本，计算对比损失
        if neg_items_emb is not None:
            high_freq_neg_sim = F.cosine_similarity(high_freq_repr, neg_items_emb, dim=1)
            high_freq_logits = torch.stack([high_freq_pos_sim, high_freq_neg_sim], dim=1)
            high_freq_labels = torch.zeros(batch_size, dtype=torch.long, device=high_freq_repr.device)
            high_freq_contrast = self.contrast_loss(high_freq_logits * 5, high_freq_labels)
            return high_freq_loss + high_freq_contrast
        else:
            return high_freq_loss
    
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
        
        # 获取样本权重（如果存在核心集权重）
        if self.has_coreset_weight and 'dtw_coreset_weight' in interaction:
            sample_weights = interaction['dtw_coreset_weight'].float()
            # 标准化权重，确保平均为1
            sample_weights = sample_weights / (sample_weights.mean() + 1e-8)
        else:
            sample_weights = None
        
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
            
            # 应用样本权重（如果有）
            if sample_weights is not None:
                try:
                    main_loss = self.loss_fct(pos_score, neg_score, sample_weights)
                except TypeError:
                    # 如果BPRLoss不接受权重参数，则使用不带权重的版本
                    main_loss = self.loss_fct(pos_score, neg_score)
                    # 手动应用权重
                    if hasattr(main_loss, 'item'):
                        main_loss = (main_loss.item() * sample_weights).mean()
            else:
                main_loss = self.loss_fct(pos_score, neg_score)
            
            # 计算简化版对比学习损失
            contrastive_loss = self.calculate_simplified_contrast_loss(
                high_freq_repr, pos_items_emb, neg_items_emb
            )
        else:  # CE损失
            test_item_emb = self.item_embedding.weight
            logits = torch.matmul(seq_output, test_item_emb.transpose(0, 1))
            
            # 应用样本权重（如果有）
            if sample_weights is not None:
                main_loss = F.cross_entropy(logits, pos_items, reduction='none')
                main_loss = (main_loss * sample_weights).mean()
            else:
                main_loss = self.loss_fct(logits, pos_items)
            
            # 对比学习损失
            pos_items_emb = self.item_embedding(pos_items)
            contrastive_loss = self.calculate_simplified_contrast_loss(
                high_freq_repr, pos_items_emb
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
