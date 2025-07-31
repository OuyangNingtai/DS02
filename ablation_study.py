# -*- coding: utf-8 -*-
# @Time    : 2025/1/31
# @Author  : Generated for FRF Ablation Study
# @Description: Ablation study for FRF model components

r"""
FRF Ablation Study
################################################

This script provides ablation study variants of the FRF model to evaluate
the effectiveness of different components:

1. FRF-w/o-Freq: Remove frequency domain decoupling
2. FRF-w/o-DTW: Remove DTW repetition tracker
3. FRF-w/o-Fusion: Remove unified recommender fusion
4. FRF-w/o-HighFreq: Remove high frequency enhancement
5. FRF-Fixed-Weight: Use fixed repetition weight

Usage:
    python ablation_study.py --variant FRF-w/o-Freq --dataset ml-1m
    python ablation_study.py --variant FRF-w/o-DTW --dataset ml-100k
"""

import argparse
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
from recbole.quick_start import run_recbole


class SimpleSequenceEncoder(nn.Module):
    """Simple sequence encoder without frequency domain processing"""
    def __init__(self, embedding_dim):
        super(SimpleSequenceEncoder, self).__init__()
        self.embedding_dim = embedding_dim
        self.encoder = nn.Sequential(
            nn.Linear(embedding_dim, embedding_dim),
            nn.LayerNorm(embedding_dim),
            nn.GELU()
        )
    
    def forward(self, seq_output):
        """Simple sequence encoding without frequency analysis"""
        return self.encoder(seq_output), seq_output, None


class FrequencyDomainDecouplerAblated(nn.Module):
    """Frequency domain decoupler with ablation options"""
    def __init__(self, embedding_dim, cutoff_ratio=0.6, high_freq_scale=1.5, 
                 enable_freq_domain=True, enable_high_freq_enhancement=True):
        super(FrequencyDomainDecouplerAblated, self).__init__()
        self.embedding_dim = embedding_dim
        self.cutoff_ratio = cutoff_ratio
        self.high_freq_scale = high_freq_scale if enable_high_freq_enhancement else 1.0
        self.enable_freq_domain = enable_freq_domain
        
        # Feature encoder
        self.encoder = nn.Sequential(
            nn.Linear(embedding_dim, embedding_dim),
            nn.LayerNorm(embedding_dim),
            nn.GELU()
        )
        
        # High frequency projector
        self.high_freq_projector = nn.Sequential(
            nn.Linear(embedding_dim, embedding_dim),
            nn.LayerNorm(embedding_dim)
        )
        
        # Low frequency projector
        self.low_freq_projector = nn.Sequential(
            nn.Linear(embedding_dim, embedding_dim),
            nn.LayerNorm(embedding_dim)
        )
    
    def forward(self, seq_output):
        """Forward pass with ablation options"""
        # Safety check
        if torch.isnan(seq_output).any():
            seq_output = torch.nan_to_num(seq_output, nan=0.0, posinf=1.0, neginf=-1.0)

        # If frequency domain is disabled, return simple representations
        if not self.enable_freq_domain:
            encoded = self.encoder(seq_output)
            # Return the same representation for both high and low freq
            return encoded, encoded, None

        # Standard frequency domain processing
        processed_seq = seq_output.unsqueeze(1)  # [B, 1, H]
        encoded = self.encoder(processed_seq)  # [B, 1, H]

        if torch.isnan(encoded).any():
            encoded = torch.nan_to_num(encoded, nan=0.0, posinf=1.0, neginf=-1.0)

        # Normalization
        encoded_mean = encoded.mean(dim=2, keepdim=True)
        encoded_std = encoded.std(dim=2, keepdim=True) + 1e-6
        encoded_norm = (encoded - encoded_mean) / encoded_std
        
        # Convert to frequency domain
        freq_domain = torch.fft.rfft(encoded_norm, dim=2, norm="ortho")  # [B, 1, H//2+1]
        
        # High-low frequency split
        cutoff = int(freq_domain.shape[2] * self.cutoff_ratio)
        
        # High frequency components (repetitive interests) - with enhancement
        high_freq_components = freq_domain.clone()
        high_freq_components[:, :, cutoff:] = 0  # Zero out low freq
        high_freq_components *= self.high_freq_scale  # Enhancement
        
        # Low frequency components (exploratory interests)
        low_freq_components = freq_domain.clone()
        low_freq_components[:, :, :cutoff] = 0  # Zero out high freq
        
        # Inverse transform back to time domain
        high_freq_repr = torch.fft.irfft(high_freq_components, n=self.embedding_dim, dim=2, norm="ortho")
        low_freq_repr = torch.fft.irfft(low_freq_components, n=self.embedding_dim, dim=2, norm="ortho")

        # Safety check
        if torch.isnan(high_freq_repr).any() or torch.isnan(low_freq_repr).any():
            high_freq_repr = torch.nan_to_num(high_freq_repr, nan=0.0, posinf=1.0, neginf=-1.0)
            low_freq_repr = torch.nan_to_num(low_freq_repr, nan=0.0, posinf=1.0, neginf=-1.0)
        
        # Squeeze dimensions
        high_freq_repr = high_freq_repr.squeeze(1)  # [B, H]
        low_freq_repr = low_freq_repr.squeeze(1)  # [B, H]
        
        # Project to representation space
        high_freq_repr = self.high_freq_projector(high_freq_repr)
        low_freq_repr = self.low_freq_projector(low_freq_repr)

        return high_freq_repr, low_freq_repr, freq_domain.squeeze(1)


class SimpleRepetitionTracker(nn.Module):
    """Simple repetition tracker without DTW"""
    def __init__(self, embedding_dim):
        super(SimpleRepetitionTracker, self).__init__()
        self.embedding_dim = embedding_dim
        
        # Simple pattern detector
        self.pattern_detector = nn.Sequential(
            nn.Linear(embedding_dim, embedding_dim // 2),
            nn.LayerNorm(embedding_dim // 2),
            nn.GELU(),
            nn.Linear(embedding_dim // 2, 1),
            nn.Sigmoid()
        )
    
    def forward(self, batch_size, device):
        """Simple repetition detection"""
        # Return fixed repetition scores without DTW
        return torch.ones(batch_size, 1, device=device) * 0.5


class DTWRepetitionTrackerAblated(nn.Module):
    """DTW repetition tracker with ablation options"""
    def __init__(self, embedding_dim, window_size=3, rep_threshold=0.6, enable_dtw=True):
        super(DTWRepetitionTrackerAblated, self).__init__()
        self.embedding_dim = embedding_dim
        self.window_size = window_size
        self.rep_threshold = rep_threshold
        self.enable_dtw = enable_dtw
        
        # DTW computation cache
        self.dtw_cache = {}
        self.max_cache_size = 50
        
        # Dummy parameter to ensure module has parameters
        self.dummy_param = nn.Parameter(torch.zeros(1))
        
        # Dynamic pattern detector
        self.pattern_detector = None
        self._input_dim = None
        
        # Use list to store historical frequency features
        self.freq_history_list = deque(maxlen=window_size * 10)
    
    def _build_pattern_detector(self, input_dim, device):
        """Dynamically build pattern detector"""
        if self.pattern_detector is None:
            self._input_dim = input_dim
            hidden_dim = max(input_dim // 4, 16)
            
            self.pattern_detector = nn.Sequential(
                nn.Linear(input_dim, hidden_dim),
                nn.LayerNorm(hidden_dim),
                nn.GELU(),
                nn.Linear(hidden_dim, 1),
                nn.Sigmoid()
            ).to(device)
            
            # Register to module correctly
            self.add_module('pattern_detector', self.pattern_detector)
            
            # Initialize weights
            for module in self.pattern_detector.modules():
                if isinstance(module, nn.Linear):
                    nn.init.xavier_normal_(module.weight)
                    if module.bias is not None:
                        nn.init.zeros_(module.bias)
    
    def update_freq_history(self, freq_features):
        """Update frequency history"""
        if not self.training or not self.enable_dtw:
            return
            
        # Add current batch features to history list
        with torch.no_grad():
            for i in range(freq_features.shape[0]):
                magnitude = torch.abs(freq_features[i])
                self.freq_history_list.append(magnitude.cpu().clone())
    
    def get_freq_history_tensor(self, batch_size, device):
        """Get frequency history tensor"""
        if len(self.freq_history_list) < 2:
            return None
        
        # Take recent history features
        history_len = min(self.window_size, len(self.freq_history_list))
        recent_history = list(self.freq_history_list)[-history_len:]
        
        # Build pseudo history tensor for each batch sample
        freq_dim = recent_history[0].shape[0]
        history_tensor = torch.zeros(batch_size, history_len, freq_dim, device=device)
        
        for t, hist_feat in enumerate(recent_history):
            history_tensor[:, t, :] = hist_feat.to(device)
        
        return history_tensor
    
    def dtw_distance(self, seq1: torch.Tensor, seq2: torch.Tensor) -> torch.Tensor:
        """Compute DTW distance (simplified version)"""
        if not self.enable_dtw:
            # Return simple L2 distance if DTW is disabled
            return torch.norm(seq1.mean(dim=0) - seq2.mean(dim=0), p=2)
            
        # Check cache
        try:
            cache_key = (hash(seq1.cpu().detach().numpy().tobytes()), 
                         hash(seq2.cpu().detach().numpy().tobytes()))
            if cache_key in self.dtw_cache:
                return self.dtw_cache[cache_key]
        except:
            cache_key = None
            
        n, m = seq1.size(0), seq2.size(0)
        
        # Short sequence optimization
        if n <= 2 or m <= 2:
            if n == m:
                dist = torch.norm(seq1 - seq2, p=2) / max(n, 1)
            else:
                dist = torch.norm(seq1.mean(dim=0) - seq2.mean(dim=0), p=2)
        else:
            # Full DTW computation
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
        
        # Update cache
        if cache_key is not None and len(self.dtw_cache) < self.max_cache_size:
            self.dtw_cache[cache_key] = dist
            
        return dist
    
    def compute_repetition_score(self, freq_history):
        """Compute repetition intensity score"""
        if freq_history is None:
            batch_size = 1
            device = self.dummy_param.device
            return torch.ones(batch_size, 1, device=device) * 0.5
        
        batch_size = freq_history.shape[0]
        device = freq_history.device
        
        if freq_history.shape[1] < 2:
            return torch.ones(batch_size, 1, device=device) * 0.5
        
        # Extract high frequency part
        cutoff = int(freq_history.shape[2] * 0.6)
        high_freq_history = freq_history[:, :, :cutoff]  # Already magnitude spectrum
        
        # Normalization
        mean = high_freq_history.mean(dim=(1,2), keepdim=True)
        std = high_freq_history.std(dim=(1,2), keepdim=True) + 1e-6
        high_freq_history = (high_freq_history - mean) / std
        
        # Compute latest features
        latest_features = high_freq_history[:, -1, :]  # [B, F]
        
        # Dynamically build pattern detector
        if self.pattern_detector is None:
            self._build_pattern_detector(latest_features.shape[1], device)
        
        # Use pattern detector
        repetition_scores = self.pattern_detector(latest_features)
        
        # DTW similarity analysis (only if DTW is enabled)
        if self.enable_dtw and freq_history.shape[1] >= 2:
            dtw_similarities = []
            for i in range(batch_size):
                try:
                    recent_seq = high_freq_history[i, -1, :].unsqueeze(0)
                    past_seq = high_freq_history[i, :-1, :].mean(dim=0, keepdim=True)
                    
                    dtw_dist = self.dtw_distance(recent_seq, past_seq)
                    similarity = torch.exp(-dtw_dist / self.embedding_dim)
                    dtw_similarities.append(similarity)
                except:
                    dtw_similarities.append(torch.tensor(0.5, device=device))
            
            dtw_sim = torch.stack(dtw_similarities).unsqueeze(1)
            repetition_scores = 0.7 * repetition_scores + 0.3 * dtw_sim
        
        return repetition_scores
    
    def forward(self, batch_size, device):
        """Forward pass"""
        freq_history = self.get_freq_history_tensor(batch_size, device)
        
        if freq_history is not None and torch.isnan(freq_history).any():
            freq_history = torch.nan_to_num(freq_history, nan=0.0, posinf=1.0, neginf=-1.0)
        
        repetition_scores = self.compute_repetition_score(freq_history)
        
        # Ensure correct batch size
        if repetition_scores.shape[0] != batch_size:
            repetition_scores = torch.ones(batch_size, 1, device=device) * 0.5
        
        if torch.isnan(repetition_scores).any():
            repetition_scores = torch.nan_to_num(repetition_scores, nan=0.5, posinf=1.0, neginf=0.0)
        
        return repetition_scores


class UnifiedRecommenderAblated(nn.Module):
    """Unified recommender with ablation options"""
    def __init__(self, embedding_dim, rep_weight=0.7, enable_fusion=True, fixed_weight=False):
        super(UnifiedRecommenderAblated, self).__init__()
        self.embedding_dim = embedding_dim
        self.rep_weight = rep_weight
        self.enable_fusion = enable_fusion
        self.fixed_weight = fixed_weight
        
        # Unified fusion network
        if enable_fusion and not fixed_weight:
            self.fusion_network = nn.Sequential(
                nn.Linear(embedding_dim * 2 + 1, embedding_dim),
                nn.LayerNorm(embedding_dim),
                nn.GELU(),
                nn.Linear(embedding_dim, 1),
                nn.Sigmoid()
            )
        
        # Final projection
        self.final_projection = nn.Sequential(
            nn.Linear(embedding_dim * 2, embedding_dim),
            nn.LayerNorm(embedding_dim)
        )
    
    def forward(self, high_freq_repr, low_freq_repr, repetition_scores):
        """Fuse high frequency and low frequency representations"""
        if not self.enable_fusion:
            # Simple concatenation without fusion
            concat_repr = torch.cat([high_freq_repr, low_freq_repr], dim=1)
            return self.final_projection(concat_repr)
        
        if self.fixed_weight:
            # Fixed weight fusion
            repetition_weight = self.rep_weight
            novelty_weight = 1 - repetition_weight
        else:
            # Dynamic weight computation
            fusion_input = torch.cat([high_freq_repr, low_freq_repr, repetition_scores], dim=1)
            dynamic_adjustment = self.fusion_network(fusion_input)
            
            # Compute final weights
            repetition_weight = self.rep_weight + (1 - self.rep_weight) * dynamic_adjustment
            novelty_weight = 1 - repetition_weight
        
        # Weighted fusion
        weighted_high = repetition_weight * high_freq_repr
        weighted_low = novelty_weight * low_freq_repr
        
        # Final representation
        concat_repr = torch.cat([weighted_high, weighted_low], dim=1)
        final_repr = self.final_projection(concat_repr)
        
        return final_repr


class FRFAblation(SequentialRecommender):
    """FRF model with ablation options"""
    
    def __init__(self, config, dataset, variant='full'):
        super(FRFAblation, self).__init__(config, dataset)
        
        # Ablation variant configuration
        self.variant = variant
        self.configure_ablation_variant()
        
        # Basic parameters
        self.embedding_size = config['embedding_size']
        self.hidden_size = config['hidden_size']
        self.n_layers = config['n_layers']
        self.n_heads = config['n_heads']
        self.inner_size = config['inner_size']
        self.hidden_dropout_prob = config['hidden_dropout_prob']
        self.attn_dropout_prob = config['attn_dropout_prob']
        self.hidden_act = config['hidden_act']
        self.layer_norm_eps = config['layer_norm_eps']
        
        # Core parameters
        self.cutoff_ratio = config['cutoff_ratio'] if 'cutoff_ratio' in config else 0.65 
        self.rep_weight = config['rep_weight'] if 'rep_weight' in config else 0.75       
        self.high_freq_scale = config['high_freq_scale'] if 'high_freq_scale' in config else 1.8 
        self.window_size = config['window_size'] if 'window_size' in config else 7     
        self.rep_threshold = config['rep_threshold'] if 'rep_threshold' in config else 0.6
        self.contrast_weight = config['contrast_weight'] if 'contrast_weight' in config else 0.25  
        
        # Training enhancement parameters
        self.label_smoothing = config['label_smoothing'] if 'label_smoothing' in config else 0.1
        self.temperature = config['temperature'] if 'temperature' in config else 0.07
        self.freq_reg_weight = config['freq_reg_weight'] if 'freq_reg_weight' in config else 0.01
        
        # Initialize embedding layers
        self.item_embedding = nn.Embedding(self.n_items, self.embedding_size, padding_idx=0)
        self.position_embedding = nn.Embedding(self.max_seq_length, self.embedding_size)
        
        # Sequence encoder
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
            self.encoder = nn.LSTM(
                input_size=self.embedding_size,
                hidden_size=self.hidden_size,
                num_layers=1,
                batch_first=True
            )
        
        # Initialize ablated components based on variant
        self.init_ablated_components()
        
        # Basic components
        self.LayerNorm = nn.LayerNorm(self.embedding_size, eps=self.layer_norm_eps)
        self.dropout = nn.Dropout(self.hidden_dropout_prob)
        
        # Loss functions
        self.loss_type = config['loss_type']
        if self.loss_type == 'BPR':
            self.loss_fct = BPRLoss()
        elif self.loss_type == 'CE':
            self.loss_fct = nn.CrossEntropyLoss()
        else:
            raise NotImplementedError("Make sure 'loss_type' in ['BPR', 'CE']!")
        
        # Parameter initialization
        self.apply(self._init_weights)
    
    def configure_ablation_variant(self):
        """Configure ablation variant settings"""
        if self.variant == 'FRF-w/o-Freq':
            self.enable_freq_domain = False
            self.enable_dtw = True
            self.enable_fusion = True
            self.enable_high_freq_enhancement = True
            self.fixed_weight = False
        elif self.variant == 'FRF-w/o-DTW':
            self.enable_freq_domain = True
            self.enable_dtw = False
            self.enable_fusion = True
            self.enable_high_freq_enhancement = True
            self.fixed_weight = False
        elif self.variant == 'FRF-w/o-Fusion':
            self.enable_freq_domain = True
            self.enable_dtw = True
            self.enable_fusion = False
            self.enable_high_freq_enhancement = True
            self.fixed_weight = False
        elif self.variant == 'FRF-w/o-HighFreq':
            self.enable_freq_domain = True
            self.enable_dtw = True
            self.enable_fusion = True
            self.enable_high_freq_enhancement = False
            self.fixed_weight = False
        elif self.variant == 'FRF-Fixed-Weight':
            self.enable_freq_domain = True
            self.enable_dtw = True
            self.enable_fusion = True
            self.enable_high_freq_enhancement = True
            self.fixed_weight = True
        else:  # Full model
            self.enable_freq_domain = True
            self.enable_dtw = True
            self.enable_fusion = True
            self.enable_high_freq_enhancement = True
            self.fixed_weight = False
    
    def init_ablated_components(self):
        """Initialize components based on ablation settings"""
        # Frequency domain decoupler
        if self.enable_freq_domain:
            self.freq_decoupler = FrequencyDomainDecouplerAblated(
                self.embedding_size, 
                self.cutoff_ratio,
                self.high_freq_scale,
                enable_freq_domain=True,
                enable_high_freq_enhancement=self.enable_high_freq_enhancement
            )
        else:
            self.freq_decoupler = SimpleSequenceEncoder(self.embedding_size)
        
        # Repetition tracker
        if self.enable_dtw:
            self.repetition_tracker = DTWRepetitionTrackerAblated(
                self.embedding_size,
                self.window_size,
                self.rep_threshold,
                enable_dtw=True
            )
        else:
            self.repetition_tracker = SimpleRepetitionTracker(self.embedding_size)
        
        # Unified recommender
        self.recommender = UnifiedRecommenderAblated(
            self.embedding_size,
            self.rep_weight,
            enable_fusion=self.enable_fusion,
            fixed_weight=self.fixed_weight
        )
    
    def _init_weights(self, module):
        """Initialize weights"""
        if isinstance(module, (nn.Linear, nn.Embedding)):
            module.weight.data.normal_(mean=0.0, std=0.02)
        elif isinstance(module, nn.LayerNorm):
            module.bias.data.zero_()
            module.weight.data.fill_(1.0)
        if isinstance(module, nn.Linear) and module.bias is not None:
            module.bias.data.zero_()
    
    def get_attention_mask(self, item_seq):
        """Generate attention mask"""
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
        """Forward pass"""
        # Sequence encoding
        position_ids = torch.arange(item_seq.size(1), dtype=torch.long, device=item_seq.device)
        position_ids = position_ids.unsqueeze(0).expand_as(item_seq)
        position_embedding = self.position_embedding(position_ids)
        
        item_emb = self.item_embedding(item_seq)
        input_emb = item_emb + position_embedding
        input_emb = self.LayerNorm(input_emb)
        input_emb = self.dropout(input_emb)
        
        # Encoder processing
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
        
        # Get last timestep representation
        output = self.gather_indexes(seq_output, item_seq_len - 1)
        if torch.isnan(output).any():
            output = torch.nan_to_num(output, nan=0.0, posinf=1.0, neginf=-1.0)
        
        # Frequency domain decoupling (or simple encoding)
        high_freq_repr, low_freq_repr, freq_domain = self.freq_decoupler(output)
        
        # Update frequency history (if DTW is enabled)
        if hasattr(self.repetition_tracker, 'update_freq_history') and freq_domain is not None:
            self.repetition_tracker.update_freq_history(freq_domain)
        
        # Repetition interest detection
        batch_size = high_freq_repr.shape[0]
        device = high_freq_repr.device
        
        if hasattr(self.repetition_tracker, 'forward'):
            repetition_scores = self.repetition_tracker(batch_size, device)
        else:
            repetition_scores = torch.ones(batch_size, 1, device=device) * 0.5
        
        # Fusion recommendation
        final_repr = self.recommender(high_freq_repr, low_freq_repr, repetition_scores)
        
        return final_repr
    
    def calculate_loss(self, interaction):
        """Enhanced loss function"""
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
            
            # Label Smoothing for BPR
            margin = 1.0 - self.label_smoothing + self.label_smoothing * torch.rand_like(pos_score)
            bpr_loss = -torch.log(torch.sigmoid(pos_score - neg_score + margin) + 1e-8).mean()
            
            # Add contrast loss only for variants that support it
            if self.enable_freq_domain and hasattr(self.freq_decoupler, 'high_freq_projector'):
                high_freq_repr, _, freq_domain = self.freq_decoupler(seq_output)
                
                # Temperature-scaled contrastive learning
                pos_sim = F.cosine_similarity(high_freq_repr, pos_items_emb, dim=1) / self.temperature
                neg_sim = F.cosine_similarity(high_freq_repr, neg_items_emb, dim=1) / self.temperature
                
                # InfoNCE-style contrastive loss
                contrast_loss = -torch.log(
                    torch.exp(pos_sim) / (torch.exp(pos_sim) + torch.exp(neg_sim) + 1e-8)
                ).mean()
                
                # Frequency domain regularization
                if freq_domain is not None:
                    freq_magnitude = torch.abs(freq_domain)
                    freq_reg = torch.mean(freq_magnitude) * self.freq_reg_weight
                else:
                    freq_reg = 0.0
                
                return bpr_loss + self.contrast_weight * contrast_loss + freq_reg
            else:
                return bpr_loss
            
        else:  # CE loss
            test_item_emb = self.item_embedding.weight
            logits = torch.matmul(seq_output, test_item_emb.transpose(0, 1))
            
            # Label smoothing for CE
            num_classes = logits.size(1)
            smoothed_labels = torch.zeros_like(logits)
            smoothed_labels.scatter_(1, pos_items.unsqueeze(1), 1.0 - self.label_smoothing)
            smoothed_labels += self.label_smoothing / num_classes
            
            ce_loss = F.kl_div(F.log_softmax(logits, dim=1), smoothed_labels, reduction='batchmean')
            
            # Frequency domain regularization
            if self.enable_freq_domain and hasattr(self.freq_decoupler, 'high_freq_projector'):
                _, _, freq_domain = self.freq_decoupler(seq_output)
                if freq_domain is not None:
                    freq_magnitude = torch.abs(freq_domain)
                    freq_reg = torch.mean(freq_magnitude) * self.freq_reg_weight
                else:
                    freq_reg = 0.0
                return ce_loss + freq_reg
            else:
                return ce_loss
    
    def predict(self, interaction):
        """Prediction"""
        item_seq = interaction[self.ITEM_SEQ]
        item_seq_len = interaction[self.ITEM_SEQ_LEN]
        test_item = interaction[self.ITEM_ID]
        
        seq_output = self.forward(item_seq, item_seq_len)
        test_item_emb = self.item_embedding(test_item)
        scores = torch.mul(seq_output, test_item_emb).sum(dim=1)
        
        return scores
    
    def full_sort_predict(self, interaction):
        """Full sort prediction"""
        item_seq = interaction[self.ITEM_SEQ]
        item_seq_len = interaction[self.ITEM_SEQ_LEN]
        
        seq_output = self.forward(item_seq, item_seq_len)
        test_items_emb = self.item_embedding.weight
        scores = torch.matmul(seq_output, test_items_emb.transpose(0, 1))
        
        return scores


def run_ablation_study(variant, dataset, config_file='test.yaml', additional_params=None):
    """Run ablation study for a specific variant"""
    print(f"\n{'='*60}")
    print(f"Running Ablation Study: {variant}")
    print(f"Dataset: {dataset}")
    print(f"{'='*60}")
    
    # Model parameters
    model_params = {
        'variant': variant
    }
    
    # Add any additional parameters
    if additional_params:
        model_params.update(additional_params)
    
    # Configuration files
    config_file_list = [config_file] if config_file else None
    
    try:
        # Run the experiment
        result = run_recbole(
            model='FRFAblation',
            dataset=dataset,
            config_file_list=config_file_list,
            config_dict=model_params
        )
        
        print(f"\nResults for {variant}:")
        print(result)
        return result
        
    except Exception as e:
        print(f"Error running {variant}: {str(e)}")
        return None


def run_all_ablation_variants(dataset='ml-1m', config_file='test.yaml'):
    """Run all ablation study variants"""
    variants = [
        'FRF-Full',
        'FRF-w/o-Freq',
        'FRF-w/o-DTW', 
        'FRF-w/o-Fusion',
        'FRF-w/o-HighFreq',
        'FRF-Fixed-Weight'
    ]
    
    results = {}
    
    for variant in variants:
        print(f"\n\nStarting ablation study for: {variant}")
        result = run_ablation_study(variant, dataset, config_file)
        results[variant] = result
        
        # Add some delay between experiments
        time.sleep(2)
    
    return results


def generate_latex_table(results):
    """Generate LaTeX table from ablation study results"""
    # Extract metrics (assuming RecBole standard metrics)
    metrics = ['Recall@10', 'NDCG@10', 'MRR@10', 'Hit@10', 'Precision@10']
    
    latex_code = """
\\begin{table}[h]
\\centering
\\caption{Ablation Study Results on FRF Model Components}
\\label{tab:ablation_study}
\\begin{tabular}{l|ccccc}
\\hline
\\textbf{Model Variant} & \\textbf{Recall@10} & \\textbf{NDCG@10} & \\textbf{MRR@10} & \\textbf{Hit@10} & \\textbf{Precision@10} \\\\
\\hline
"""
    
    for variant, result in results.items():
        if result is not None:
            row = f"{variant.replace('_', '\\_')}"
            for metric in metrics:
                if metric in result:
                    row += f" & {result[metric]:.4f}"
                else:
                    row += " & -"
            row += " \\\\\n"
            latex_code += row
    
    latex_code += """\\hline
\\end{tabular}
\\end{table}
"""
    
    return latex_code


def main():
    """Main function for ablation study"""
    parser = argparse.ArgumentParser(description='FRF Ablation Study')
    parser.add_argument('--variant', '-v', type=str, default='all', 
                       choices=['all', 'FRF-Full', 'FRF-w/o-Freq', 'FRF-w/o-DTW', 
                               'FRF-w/o-Fusion', 'FRF-w/o-HighFreq', 'FRF-Fixed-Weight'],
                       help='Ablation variant to run')
    parser.add_argument('--dataset', '-d', type=str, default='ml-1m', 
                       choices=['ml-1m', 'ml-100k'],
                       help='Dataset to use')
    parser.add_argument('--config', '-c', type=str, default='test.yaml',
                       help='Configuration file')
    parser.add_argument('--output', '-o', type=str, default='ablation_results.txt',
                       help='Output file for results')
    
    args = parser.parse_args()
    
    print("FRF Model Ablation Study")
    print("=" * 50)
    print(f"Dataset: {args.dataset}")
    print(f"Config: {args.config}")
    print(f"Variant: {args.variant}")
    print("=" * 50)
    
    if args.variant == 'all':
        # Run all variants
        results = run_all_ablation_variants(args.dataset, args.config)
        
        # Generate LaTeX table
        latex_table = generate_latex_table(results)
        
        # Save results
        with open(args.output, 'w') as f:
            f.write("FRF Ablation Study Results\n")
            f.write("=" * 50 + "\n\n")
            
            for variant, result in results.items():
                f.write(f"{variant}:\n")
                f.write(str(result) + "\n\n")
            
            f.write("\nLaTeX Table:\n")
            f.write(latex_table)
        
        print(f"\nResults saved to: {args.output}")
        print("\nLaTeX Table:")
        print(latex_table)
        
    else:
        # Run single variant
        result = run_ablation_study(args.variant, args.dataset, args.config)
        
        # Save result
        with open(args.output, 'w') as f:
            f.write(f"FRF Ablation Study - {args.variant}\n")
            f.write("=" * 50 + "\n")
            f.write(str(result))
        
        print(f"\nResult saved to: {args.output}")


if __name__ == '__main__':
    main()