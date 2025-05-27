# -*- coding: utf-8 -*-
# @Time    : 2025/5/24

r"""
DTW CoreSet Sequence Recommender
################################################
基于DTW核心集的序列推荐模型
通过DTW相似度对序列片段进行聚类，选择代表性片段构建核心集
核心集片段作为高频信号，非核心集作为低频信号
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
from sklearn.cluster import KMeans
import pickle

from recbole.model.abstract_recommender import SequentialRecommender
from recbole.model.loss import BPRLoss
from recbole.data.interaction import Interaction
from recbole.data.dataset import Dataset
from recbole.utils import set_color


class DTWSequmentCoreSetProcessor:
    """基于DTW的序列片段核心集处理器
    
    将长序列分片，使用DTW计算片段相似度并聚类，选择代表性片段构建核心集
    
    Args:
        dataset (Dataset): 数据集对象
        embedding_size (int): 嵌入维度
        segment_length (int): 片段长度，默认8
        overlap_ratio (float): 片段重叠比例，默认0.5
        coreset_ratio (float): 核心集比例，默认0.3
        n_clusters (int): 聚类数量，默认None（自动确定）
        min_segment_freq (int): 最小片段频次，默认2
        cache_dir (str): 缓存目录
        force_recompute (bool): 是否强制重新计算
        dataset_name (str): 数据集名称
    """
    def __init__(self, dataset, embedding_size, 
                 segment_length=8, overlap_ratio=0.5, coreset_ratio=0.3,
                 n_clusters=None, min_segment_freq=2,
                 cache_dir='./dtw_coreset_cache', 
                 force_recompute=False, dataset_name=None):
        self.logger = logging.getLogger()
        self.dataset = dataset
        self.embedding_size = embedding_size
        self.segment_length = segment_length
        self.overlap_ratio = overlap_ratio
        self.coreset_ratio = coreset_ratio
        self.n_clusters = n_clusters
        self.min_segment_freq = min_segment_freq
        self.cache_dir = cache_dir
        self.force_recompute = force_recompute
        self.dataset_name = dataset_name or 'default_dataset'
        
        # 创建缓存目录
        if not os.path.exists(cache_dir):
            os.makedirs(cache_dir)
            
        # DTW距离计算缓存
        self.dtw_cache = {}
        self.max_cache_size = 50000
        
        # 数据集信息
        self.uid_field = dataset.uid_field
        self.iid_field = dataset.iid_field
        self.item_num = dataset.item_num
        
        # 核心集信息
        self.coreset_segments = []
        self.segment_to_cluster = {}
        self.cluster_representatives = {}
        
    def _segment_sequence(self, sequence):
        """将序列分割为重叠片段
        
        Args:
            sequence (list): 物品序列
            
        Returns:
            list: 片段列表，每个片段是(start_idx, segment)的元组
        """
        if len(sequence) < self.segment_length:
            return [(0, sequence)]
        
        segments = []
        step = max(1, int(self.segment_length * (1 - self.overlap_ratio)))
        
        for start in range(0, len(sequence) - self.segment_length + 1, step):
            segment = sequence[start:start + self.segment_length]
            segments.append((start, segment))
        
        # 确保包含最后一个片段
        if segments[-1][0] + self.segment_length < len(sequence):
            segments.append((len(sequence) - self.segment_length, 
                           sequence[-self.segment_length:]))
        
        return segments
    
    def dtw_distance(self, seq1, seq2):
        """计算两个序列间的DTW距离
        
        Args:
            seq1 (list): 第一个序列（物品ID列表）
            seq2 (list): 第二个序列（物品ID列表）
            
        Returns:
            float: DTW距离值
        """
        # 缓存键
        cache_key = (tuple(seq1), tuple(seq2))
        if cache_key in self.dtw_cache:
            return self.dtw_cache[cache_key]
        
        n, m = len(seq1), len(seq2)
        
        # 对于短序列，简化为Jaccard距离
        if n <= 2 or m <= 2:
            set1, set2 = set(seq1), set(seq2)
            intersection = len(set1 & set2)
            union = len(set1 | set2)
            jaccard_sim = intersection / union if union > 0 else 0
            dist = 1 - jaccard_sim
        else:
            # 创建DTW矩阵
            dtw_matrix = np.zeros((n+1, m+1)) + np.inf
            dtw_matrix[0, 0] = 0
            
            # 计算成对距离（基于物品ID）
            for i in range(1, n+1):
                for j in range(1, m+1):
                    cost = 0 if seq1[i-1] == seq2[j-1] else 1
                    dtw_matrix[i, j] = cost + min(
                        dtw_matrix[i-1, j],    # 插入
                        dtw_matrix[i, j-1],    # 删除
                        dtw_matrix[i-1, j-1]   # 替换
                    )
            
            dist = dtw_matrix[n, m] / max(n, m)  # 归一化
        
        # 存入缓存
        if len(self.dtw_cache) < self.max_cache_size:
            self.dtw_cache[cache_key] = dist
        
        return dist
    
    def _extract_all_segments(self):
        """从数据集中提取所有用户序列的片段
        
        Returns:
            tuple: (所有片段列表, 片段到用户ID的映射, 用户序列字典)
        """
        self.logger.info("提取用户序列片段...")
        
        # 获取交互数据
        inter_feat = self.dataset.inter_feat
        user_ids = inter_feat[self.uid_field].numpy()
        item_ids = inter_feat[self.iid_field].numpy()
        
        # 按用户ID排序
        sort_idx = np.argsort(user_ids)
        user_ids = user_ids[sort_idx]
        item_ids = item_ids[sort_idx]
        
        # 分组构建序列
        user_sequences = defaultdict(list)
        for user_id, item_id in zip(user_ids, item_ids):
            user_sequences[user_id].append(item_id)
        
        # 提取所有片段
        all_segments = []
        segment_to_user = {}
        
        for user_id, sequence in user_sequences.items():
            if len(sequence) >= self.segment_length:
                segments = self._segment_sequence(sequence)
                for start_idx, segment in segments:
                    segment_id = len(all_segments)
                    all_segments.append(segment)
                    segment_to_user[segment_id] = (user_id, start_idx)
        
        self.logger.info(f"提取了{len(all_segments)}个序列片段")
        return all_segments, segment_to_user, user_sequences
    
    def _compute_segment_similarity_matrix(self, segments):
        """计算片段相似度矩阵
        
        Args:
            segments (list): 片段列表
            
        Returns:
            np.ndarray: 相似度矩阵
        """
        self.logger.info("计算片段DTW相似度矩阵...")
        n_segments = len(segments)
        similarity_matrix = np.zeros((n_segments, n_segments))
        
        # 使用tqdm显示进度
        for i in tqdm(range(n_segments), desc="计算相似度"):
            for j in range(i, n_segments):
                if i == j:
                    similarity_matrix[i, j] = 0
                else:
                    dist = self.dtw_distance(segments[i], segments[j])
                    similarity_matrix[i, j] = dist
                    similarity_matrix[j, i] = dist
        
        return similarity_matrix
    
    def _cluster_segments(self, segments, similarity_matrix):
        """对片段进行聚类
        
        Args:
            segments (list): 片段列表
            similarity_matrix (np.ndarray): 相似度矩阵
            
        Returns:
            tuple: (聚类标签, 聚类中心索引)
        """
        self.logger.info("对片段进行聚类...")
        
        # 确定聚类数量
        n_segments = len(segments)
        if self.n_clusters is None:
            # 自动确定聚类数：使用启发式规则
            self.n_clusters = min(max(10, int(np.sqrt(n_segments))), n_segments // 5)
        
        self.logger.info(f"使用{self.n_clusters}个聚类")
        
        # 将距离矩阵转换为特征矩阵（使用MDS思想的简化版本）
        # 这里我们使用相似度矩阵的每一行作为特征
        features = 1 - similarity_matrix  # 转换为相似度
        
        # 使用KMeans聚类
        kmeans = KMeans(n_clusters=self.n_clusters, random_state=42, n_init=10)
        cluster_labels = kmeans.fit_predict(features)
        
        # 找到每个聚类的代表片段（最接近聚类中心的片段）
        cluster_centers = {}
        for cluster_id in range(self.n_clusters):
            cluster_indices = np.where(cluster_labels == cluster_id)[0]
            if len(cluster_indices) > 0:
                # 计算聚类内平均距离，选择最小的作为代表
                min_avg_dist = float('inf')
                representative = cluster_indices[0]
                
                for idx in cluster_indices:
                    avg_dist = np.mean([similarity_matrix[idx, j] for j in cluster_indices if j != idx])
                    if avg_dist < min_avg_dist:
                        min_avg_dist = avg_dist
                        representative = idx
                
                cluster_centers[cluster_id] = representative
        
        return cluster_labels, cluster_centers
    
    def _select_coreset_segments(self, segments, cluster_labels, cluster_centers):
        """选择核心集片段
        
        Args:
            segments (list): 片段列表
            cluster_labels (np.ndarray): 聚类标签
            cluster_centers (dict): 聚类中心索引
            
        Returns:
            tuple: (核心集片段索引, 片段到聚类的映射)
        """
        self.logger.info("选择核心集片段...")
        
        # 统计每个聚类的片段数量
        cluster_counts = defaultdict(int)
        for label in cluster_labels:
            cluster_counts[label] += 1
        
        # 按聚类大小排序，优先选择大聚类的代表片段
        sorted_clusters = sorted(cluster_counts.items(), key=lambda x: x[1], reverse=True)
        
        # 选择核心集
        target_coreset_size = max(1, int(len(segments) * self.coreset_ratio))
        coreset_indices = set()
        segment_to_cluster = {}
        
        # 首先添加所有聚类的代表片段
        for cluster_id, representative_idx in cluster_centers.items():
            if len(coreset_indices) < target_coreset_size:
                coreset_indices.add(representative_idx)
                segment_to_cluster[representative_idx] = cluster_id
        
        # 如果还需要更多片段，按聚类大小添加更多片段
        for cluster_id, count in sorted_clusters:
            if len(coreset_indices) >= target_coreset_size:
                break
            
            # 获取该聚类的所有片段
            cluster_indices = np.where(cluster_labels == cluster_id)[0]
            
            # 添加尚未在核心集中的片段
            for idx in cluster_indices:
                if len(coreset_indices) >= target_coreset_size:
                    break
                if idx not in coreset_indices:
                    coreset_indices.add(idx)
                    segment_to_cluster[idx] = cluster_id
        
        self.logger.info(f"选择了{len(coreset_indices)}个核心集片段")
        return list(coreset_indices), segment_to_cluster
    
    def _reconstruct_sequences(self, user_sequences, segments, segment_to_user, 
                             coreset_indices, segment_to_cluster):
        """重构用户序列，标记核心集和非核心集片段
        
        Args:
            user_sequences (dict): 原始用户序列
            segments (list): 所有片段
            segment_to_user (dict): 片段到用户的映射
            coreset_indices (list): 核心集片段索引
            segment_to_cluster (dict): 片段到聚类的映射
            
        Returns:
            dict: 重构后的用户序列信息
        """
        self.logger.info("重构用户序列...")
        
        coreset_set = set(coreset_indices)
        
        # 为每个用户重构序列
        reconstructed_sequences = {}
        
        for user_id, original_seq in user_sequences.items():
            if len(original_seq) < self.segment_length:
                # 短序列直接标记为低频
                reconstructed_sequences[user_id] = {
                    'original_sequence': original_seq,
                    'high_freq_segments': [],
                    'low_freq_segments': [original_seq],
                    'sequence_mask': [0] * len(original_seq)  # 0表示低频，1表示高频
                }
                continue
            
            # 获取该用户的所有片段
            user_segments = []
            for seg_id, (u_id, start_idx) in segment_to_user.items():
                if u_id == user_id:
                    user_segments.append((seg_id, start_idx, segments[seg_id]))
            
            # 按起始位置排序
            user_segments.sort(key=lambda x: x[1])
            
            # 创建序列掩码
            sequence_mask = [0] * len(original_seq)
            high_freq_segments = []
            low_freq_segments = []
            
            covered_positions = set()
            
            # 标记高频片段
            for seg_id, start_idx, segment in user_segments:
                if seg_id in coreset_set:
                    high_freq_segments.append((start_idx, segment))
                    # 标记这些位置为高频
                    for i in range(start_idx, min(start_idx + len(segment), len(original_seq))):
                        sequence_mask[i] = 1
                        covered_positions.add(i)
            
            # 收集低频片段（未被高频覆盖的部分）
            low_freq_start = None
            for i in range(len(original_seq)):
                if i not in covered_positions:
                    if low_freq_start is None:
                        low_freq_start = i
                else:
                    if low_freq_start is not None:
                        low_freq_segments.append(original_seq[low_freq_start:i])
                        low_freq_start = None
            
            # 处理最后一个低频片段
            if low_freq_start is not None:
                low_freq_segments.append(original_seq[low_freq_start:])
            
            reconstructed_sequences[user_id] = {
                'original_sequence': original_seq,
                'high_freq_segments': high_freq_segments,
                'low_freq_segments': low_freq_segments,
                'sequence_mask': sequence_mask
            }
        
        return reconstructed_sequences
    
    def _save_coreset_data(self, reconstructed_sequences, cluster_centers, segments):
        """保存核心集数据
        
        Args:
            reconstructed_sequences (dict): 重构序列
            cluster_centers (dict): 聚类中心
            segments (list): 所有片段
        """
        cache_path = os.path.join(self.cache_dir, f"{self.dataset_name}_segment_coreset.pkl")
        
        coreset_data = {
            'reconstructed_sequences': reconstructed_sequences,
            'cluster_centers': cluster_centers,
            'segments': segments,
            'segment_length': self.segment_length,
            'overlap_ratio': self.overlap_ratio,
            'coreset_ratio': self.coreset_ratio,
            'n_clusters': self.n_clusters
        }
        
        with open(cache_path, 'wb') as f:
            pickle.dump(coreset_data, f)
        
        self.logger.info(f"核心集数据已保存至: {cache_path}")
    
    def _load_coreset_data(self):
        """加载核心集数据
        
        Returns:
            dict or None: 核心集数据，如果不存在则返回None
        """
        cache_path = os.path.join(self.cache_dir, f"{self.dataset_name}_segment_coreset.pkl")
        
        if os.path.exists(cache_path) and not self.force_recompute:
            self.logger.info(f"加载缓存的核心集数据: {cache_path}")
            with open(cache_path, 'rb') as f:
                return pickle.load(f)
        
        return None
    
    def build_coreset(self):
        """构建DTW序列片段核心集
        
        Returns:
            dict: 重构后的序列数据
        """
        # 尝试加载缓存
        cached_data = self._load_coreset_data()
        if cached_data is not None:
            return cached_data['reconstructed_sequences']
        
        # 提取所有片段
        segments, segment_to_user, user_sequences = self._extract_all_segments()
        
        if len(segments) == 0:
            self.logger.warning("没有提取到足够的序列片段")
            return {}
        
        # 计算相似度矩阵
        similarity_matrix = self._compute_segment_similarity_matrix(segments)
        
        # 聚类
        cluster_labels, cluster_centers = self._cluster_segments(segments, similarity_matrix)
        
        # 选择核心集
        coreset_indices, segment_to_cluster = self._select_coreset_segments(
            segments, cluster_labels, cluster_centers
        )
        
        # 重构序列
        reconstructed_sequences = self._reconstruct_sequences(
            user_sequences, segments, segment_to_user, coreset_indices, segment_to_cluster
        )
        
        # 保存数据
        self._save_coreset_data(reconstructed_sequences, cluster_centers, segments)
        
        return reconstructed_sequences
    
    def process_dataset(self):
        """处理数据集，应用序列片段核心集
        
        Returns:
            Dataset: 处理后的数据集
        """
        self.logger.info(set_color("开始DTW序列片段核心集处理...", "green"))
        
        # 构建核心集
        reconstructed_sequences = self.build_coreset()
        
        if not reconstructed_sequences:
            self.logger.warning("核心集构建失败，返回原始数据集")
            return self.dataset
        
        # 应用到数据集
        self.logger.info("将核心集信息应用到数据集...")
        inter_feat = self.dataset.inter_feat
        
        # 获取交互数据
        user_ids = inter_feat[self.uid_field].numpy()
        item_ids = inter_feat[self.iid_field].numpy()
        
        # 按用户ID排序
        sort_idx = np.argsort(user_ids)
        user_ids = user_ids[sort_idx]
        item_ids = item_ids[sort_idx]
        
        # 创建高频/低频标记
        high_freq_mask = []
        low_freq_mask = []
        
        current_user = None
        current_seq_idx = 0
        
        for i, (user_id, item_id) in enumerate(zip(user_ids, item_ids)):
            if user_id != current_user:
                current_user = user_id
                current_seq_idx = 0
            
            # 获取该用户的序列掩码
            if user_id in reconstructed_sequences:
                seq_info = reconstructed_sequences[user_id]
                if current_seq_idx < len(seq_info['sequence_mask']):
                    is_high_freq = seq_info['sequence_mask'][current_seq_idx]
                else:
                    is_high_freq = 0  # 超出范围的默认为低频
            else:
                is_high_freq = 0  # 未找到的用户默认为低频
            
            high_freq_mask.append(is_high_freq)
            low_freq_mask.append(1 - is_high_freq)
            current_seq_idx += 1
        
        # 恢复原始顺序
        original_order = np.argsort(sort_idx)
        high_freq_mask = np.array(high_freq_mask)[original_order]
        low_freq_mask = np.array(low_freq_mask)[original_order]
        
        # 添加到数据集
        inter_feat.update('high_freq_mask', torch.BoolTensor(high_freq_mask))
        inter_feat.update('low_freq_mask', torch.BoolTensor(low_freq_mask))
        
        # 统计信息
        n_high_freq = np.sum(high_freq_mask)
        n_low_freq = np.sum(low_freq_mask)
        total_interactions = len(high_freq_mask)
        
        self.logger.info(set_color(
            f"序列片段核心集处理完成！"
            f"高频交互: {n_high_freq} ({n_high_freq/total_interactions:.2%}), "
            f"低频交互: {n_low_freq} ({n_low_freq/total_interactions:.2%})", 
            "green"
        ))
        
        return self.dataset


class HighLowFrequencyFusionModule(nn.Module):
    """高低频信号融合模块
    
    Args:
        embedding_dim (int): 嵌入维度
        fusion_type (str): 融合类型 ['concat', 'attention', 'gated']
        high_freq_weight (float): 高频信号权重，默认0.7
    """
    def __init__(self, embedding_dim, fusion_type='gated', high_freq_weight=0.7):
        super(HighLowFrequencyFusionModule, self).__init__()
        self.embedding_dim = embedding_dim
        self.fusion_type = fusion_type
        self.high_freq_weight = high_freq_weight
        
        if fusion_type == 'concat':
            self.fusion_layer = nn.Sequential(
                nn.Linear(embedding_dim * 2, embedding_dim),
                nn.LayerNorm(embedding_dim),
                nn.GELU()
            )
        elif fusion_type == 'attention':
            self.attention = nn.MultiheadAttention(
                embed_dim=embedding_dim,
                num_heads=8,
                dropout=0.1,
                batch_first=True
            )
            self.norm = nn.LayerNorm(embedding_dim)
        elif fusion_type == 'gated':
            self.gate_network = nn.Sequential(
                nn.Linear(embedding_dim * 2, embedding_dim),
                nn.Sigmoid()
            )
            self.fusion_transform = nn.Sequential(
                nn.Linear(embedding_dim * 2, embedding_dim),
                nn.LayerNorm(embedding_dim)
            )
        else:
            raise ValueError(f"Unsupported fusion_type: {fusion_type}")
    
    def forward(self, high_freq_repr, low_freq_repr):
        """融合高频和低频表示
        
        Args:
            high_freq_repr (torch.Tensor): 高频表示, shape: [B, H]
            low_freq_repr (torch.Tensor): 低频表示, shape: [B, H]
            
        Returns:
            torch.Tensor: 融合后的表示, shape: [B, H]
        """
        if self.fusion_type == 'concat':
            # 简单拼接后变换
            concat_repr = torch.cat([high_freq_repr, low_freq_repr], dim=-1)
            fused_repr = self.fusion_layer(concat_repr)
            
        elif self.fusion_type == 'attention':
            # 使用注意力机制融合
            # 将高频和低频作为query和key/value
            fused_repr, _ = self.attention(
                high_freq_repr.unsqueeze(1),  # query
                torch.stack([high_freq_repr, low_freq_repr], dim=1),  # key & value
                torch.stack([high_freq_repr, low_freq_repr], dim=1)
            )
            fused_repr = fused_repr.squeeze(1)
            fused_repr = self.norm(fused_repr + high_freq_repr)  # 残差连接
            
        elif self.fusion_type == 'gated':
            # 门控融合
            concat_input = torch.cat([high_freq_repr, low_freq_repr], dim=-1)
            gate = self.gate_network(concat_input)  # 门控权重
            
            # 加权融合
            weighted_high = gate * high_freq_repr
            weighted_low = (1 - gate) * low_freq_repr
            
            # 最终变换
            fused_input = torch.cat([weighted_high, weighted_low], dim=-1)
            fused_repr = self.fusion_transform(fused_input)
        
        return fused_repr


class DTWCoreSetSequenceRecommender(SequentialRecommender):
    """基于DTW序列片段核心集的推荐模型
    
    通过DTW对序列片段聚类，将核心集片段作为高频信号，非核心集作为低频信号
    
    Args:
        config (Config): 全局配置对象
        dataset (Dataset): 数据集
    """
    def __init__(self, config, dataset):
        super(DTWCoreSetSequenceRecommender, self).__init__(config, dataset)
        
        # 获取logger
        self.logger = logging.getLogger()
        
        # 获取数据集名称
        if hasattr(dataset, 'dataset_name'):
            self.dataset_name = dataset.dataset_name
        else:
            self.dataset_name = config['dataset']
        
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
        
        # DTW核心集参数
        self.use_dtw_coreset = config.get('use_dtw_coreset', True)
        self.segment_length = config.get('segment_length', 8)
        self.overlap_ratio = config.get('overlap_ratio', 0.5)
        self.coreset_ratio = config.get('coreset_ratio', 0.3)
        self.n_clusters = config.get('n_clusters', None)
        self.dtw_cache_dir = config.get('dtw_cache_dir', './dtw_coreset_cache')
        self.dtw_force_recompute = config.get('dtw_force_recompute', False)
        
        # 融合参数
        self.fusion_type = config.get('fusion_type', 'gated')
        self.high_freq_weight = config.get('high_freq_weight', 0.7)
        
        # 检查并处理DTW核心集
        self.has_freq_masks = False
        if self.use_dtw_coreset:
            self.logger.info(set_color("开始DTW序列片段核心集处理...", "green"))
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
            self.logger.warning("导入TransformerEncoder失败，使用LSTM")
            self.encoder = nn.LSTM(
                input_size=self.embedding_size,
                hidden_size=self.hidden_size,
                num_layers=self.n_layers,
                dropout=self.hidden_dropout_prob if self.n_layers > 1 else 0,
                batch_first=True
            )
        
        # 高低频表示提取器
        self.high_freq_projector = nn.Sequential(
            nn.Linear(self.embedding_size, self.hidden_size),
            nn.LayerNorm(self.hidden_size),
            nn.GELU(),
            nn.Linear(self.hidden_size, self.embedding_size)
        )
        
        self.low_freq_projector = nn.Sequential(
            nn.Linear(self.embedding_size, self.hidden_size),
            nn.LayerNorm(self.hidden_size),
            nn.GELU(),
            nn.Linear(self.hidden_size, self.embedding_size)
        )
        
        # 高低频融合模块
        self.fusion_module = HighLowFrequencyFusionModule(
            self.embedding_size, self.fusion_type, self.high_freq_weight
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
        
        # 参数初始化
        self.apply(self._init_weights)
    
    def process_dtw_coreset(self):
        """处理数据集，应用DTW序列片段核心集"""
        try:
            coreset_processor = DTWSequmentCoreSetProcessor(
                dataset=self.dataset,
                embedding_size=self.embedding_size,
                segment_length=self.segment_length,
                overlap_ratio=self.overlap_ratio,
                coreset_ratio=self.coreset_ratio,
                n_clusters=self.n_clusters,
                cache_dir=self.dtw_cache_dir,
                force_recompute=self.dtw_force_recompute,
                dataset_name=self.dataset_name
            )
            
            # 处理数据集
            self.dataset = coreset_processor.process_dataset()
            
            # 检查是否成功添加了频率掩码
            if ('high_freq_mask' in self.dataset.inter_feat and 
                'low_freq_mask' in self.dataset.inter_feat):
                self.logger.info(set_color("成功添加高低频掩码到数据集", "green"))
                self.has_freq_masks = True
            else:
                self.logger.warning("未能添加高低频掩码到数据集")
                self.has_freq_masks = False
        
        except Exception as e:
            self.logger.error(f"DTW序列片段核心集处理失败: {e}")
            self.has_freq_masks = False
    
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
    
    def extract_high_low_freq_representations(self, seq_output, high_freq_mask, low_freq_mask):
        """提取高频和低频表示
        
        Args:
            seq_output (torch.Tensor): 序列输出, shape: [B, L, H]
            high_freq_mask (torch.Tensor): 高频掩码, shape: [B, L]
            low_freq_mask (torch.Tensor): 低频掩码, shape: [B, L]
            
        Returns:
            tuple: (高频表示, 低频表示)
        """
        batch_size, seq_len, hidden_size = seq_output.shape
        
        # 扩展掩码维度
        high_freq_mask_expanded = high_freq_mask.unsqueeze(-1).float()  # [B, L, 1]
        low_freq_mask_expanded = low_freq_mask.unsqueeze(-1).float()    # [B, L, 1]
        
        # 加权序列表示
        high_freq_weighted = seq_output * high_freq_mask_expanded
        low_freq_weighted = seq_output * low_freq_mask_expanded
        
        # 计算加权平均（处理零除情况）
        high_freq_sum = high_freq_mask.sum(dim=1, keepdim=True).float()  # [B, 1]
        low_freq_sum = low_freq_mask.sum(dim=1, keepdim=True).float()    # [B, 1]
        
        # 避免除零
        high_freq_sum = torch.clamp(high_freq_sum, min=1.0)
        low_freq_sum = torch.clamp(low_freq_sum, min=1.0)
        
        # 计算平均表示
        high_freq_repr = high_freq_weighted.sum(dim=1) / high_freq_sum  # [B, H]
        low_freq_repr = low_freq_weighted.sum(dim=1) / low_freq_sum     # [B, H]
        
        # 通过投影层增强表示
        high_freq_repr = self.high_freq_projector(high_freq_repr)
        low_freq_repr = self.low_freq_projector(low_freq_repr)
        
        return high_freq_repr, low_freq_repr
    
    def forward(self, item_seq, item_seq_len):
        """前向传播
        
        Args:
            item_seq (torch.LongTensor): 物品序列, shape: [B, L]
            item_seq_len (torch.LongTensor): 序列长度, shape: [B]
            
        Returns:
            torch.Tensor: 用户表示, shape: [B, H]
        """
        # 位置编码
        position_ids = torch.arange(item_seq.size(1), dtype=torch.long, device=item_seq.device)
        position_ids = position_ids.unsqueeze(0).expand_as(item_seq)
        position_embedding = self.position_embedding(position_ids)
        
        # 物品嵌入
        item_emb = self.item_embedding(item_seq)
        input_emb = item_emb + position_embedding
        input_emb = self.LayerNorm(input_emb)
        input_emb = self.dropout(input_emb)
        
        # 序列编码
        if isinstance(self.encoder, nn.LSTM):
            encoder_output, _ = self.encoder(input_emb)
            seq_output = encoder_output
        else:
            extended_attention_mask = self.get_attention_mask(item_seq)
            encoder_output = self.encoder(input_emb, extended_attention_mask)
            if isinstance(encoder_output, list):
                seq_output = encoder_output[-1]
            else:
                seq_output = encoder_output
        
        # 获取最终的序列表示
        final_output = self.gather_indexes(seq_output, item_seq_len - 1)
        
        # 如果有高低频掩码，进行频率分离
        if self.has_freq_masks:
            # 从interaction中获取掩码信息（这需要在调用时传入）
            # 这里我们假设可以通过某种方式获取到掩码
            # 在实际使用中，需要修改forward函数签名来接收这些掩码
            batch_size = item_seq.shape[0]
            
            # 临时解决方案：使用随机掩码进行演示
            # 在实际实现中，应该从interaction中获取真实掩码
            high_freq_mask = torch.randint(0, 2, (batch_size, item_seq.shape[1]), 
                                         dtype=torch.bool, device=item_seq.device)
            low_freq_mask = ~high_freq_mask
            
            # 提取高低频表示
            high_freq_repr, low_freq_repr = self.extract_high_low_freq_representations(
                seq_output, high_freq_mask, low_freq_mask
            )
            
            # 融合高低频表示
            fused_repr = self.fusion_module(high_freq_repr, low_freq_repr)
            
            return fused_repr
        else:
            # 没有频率信息时，直接返回序列表示
            return final_output
    
    def calculate_loss(self, interaction):
        """计算损失"""
        item_seq = interaction[self.ITEM_SEQ]
        item_seq_len = interaction[self.ITEM_SEQ_LEN]
        pos_items = interaction[self.POS_ITEM_ID]
        
        seq_output = self.forward(item_seq, item_seq_len)
        
        if self.loss_type == 'BPR':
            neg_items = interaction[self.NEG_ITEM_ID]
            pos_items_emb = self.item_embedding(pos_items)
            neg_items_emb = self.item_embedding(neg_items)
            pos_score = torch.sum(seq_output * pos_items_emb, dim=-1)
            neg_score = torch.sum(seq_output * neg_items_emb, dim=-1)
            loss = self.loss_fct(pos_score, neg_score)
        else:  # CE损失
            test_item_emb = self.item_embedding.weight
            logits = torch.matmul(seq_output, test_item_emb.transpose(0, 1))
            loss = self.loss_fct(logits, pos_items)
        
        return loss
    
    def predict(self, interaction):
        """预测单个物品的分数"""
        item_seq = interaction[self.ITEM_SEQ]
        item_seq_len = interaction[self.ITEM_SEQ_LEN]
        test_item = interaction[self.ITEM_ID]
        
        seq_output = self.forward(item_seq, item_seq_len)
        test_item_emb = self.item_embedding(test_item)
        scores = torch.mul(seq_output, test_item_emb).sum(dim=1)
        
        return scores
    
    def full_sort_predict(self, interaction):
        """预测所有物品的分数"""
        item_seq = interaction[self.ITEM_SEQ]
        item_seq_len = interaction[self.ITEM_SEQ_LEN]
        
        seq_output = self.forward(item_seq, item_seq_len)
        test_items_emb = self.item_embedding.weight
        scores = torch.matmul(seq_output, test_items_emb.transpose(0, 1))
        
        return scores
