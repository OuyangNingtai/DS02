#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2025/1/31
# @Author  : Generated for FRF Ablation Study
# @Description: Test script for ablation study components

"""
Test Script for FRF Ablation Study
==================================

This script tests the ablation study implementation without requiring
the full RecBole environment. It validates the component initialization
and basic functionality.
"""

import sys
import os
import torch
import torch.nn as nn
import numpy as np
from datetime import datetime

# Add current directory to path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

# Mock RecBole classes for testing
class MockSequentialRecommender:
    """Mock SequentialRecommender for testing"""
    def __init__(self, config, dataset):
        self.config = config
        self.dataset = dataset
        self.device = torch.device('cpu')
        
    def gather_indexes(self, output, gather_index):
        """Mock gather_indexes method"""
        batch_size = output.shape[0]
        output_dim = output.shape[-1]
        
        # Simple indexing for testing
        gathered = torch.zeros(batch_size, output_dim)
        for i in range(batch_size):
            idx = min(gather_index[i].item(), output.shape[1] - 1)
            gathered[i] = output[i, idx]
        
        return gathered


class MockBPRLoss:
    """Mock BPR loss for testing"""
    def __call__(self, pos_score, neg_score):
        return -torch.log(torch.sigmoid(pos_score - neg_score) + 1e-8).mean()


class MockDataset:
    """Mock dataset for testing"""
    def __init__(self):
        self.field2token_id = {'[PAD]': 0}


class MockConfig(dict):
    """Mock configuration for testing"""
    def __init__(self):
        super().__init__()
        self.update({
            'embedding_size': 64,
            'hidden_size': 64,
            'n_layers': 2,
            'n_heads': 2,
            'inner_size': 128,
            'hidden_dropout_prob': 0.3,
            'attn_dropout_prob': 0.3,
            'hidden_act': 'gelu',
            'layer_norm_eps': 1e-12,
            'loss_type': 'BPR',
            'cutoff_ratio': 0.65,
            'rep_weight': 0.75,
            'high_freq_scale': 1.8,
            'window_size': 7,
            'rep_threshold': 0.6,
            'contrast_weight': 0.25,
            'label_smoothing': 0.1,
            'temperature': 0.07,
            'freq_reg_weight': 0.01
        })


def test_imports():
    """Test if ablation study components can be imported"""
    print("Testing imports...")
    
    try:
        # Replace RecBole imports with mocks
        import sys
        sys.modules['recbole.model.abstract_recommender'] = type('MockModule', (), {
            'SequentialRecommender': MockSequentialRecommender
        })()
        sys.modules['recbole.model.loss'] = type('MockModule', (), {
            'BPRLoss': MockBPRLoss
        })()
        sys.modules['recbole.data.interaction'] = type('MockModule', (), {
            'Interaction': dict
        })()
        
        # Now import ablation components
        from ablation_study import (
            FrequencyDomainDecouplerAblated,
            DTWRepetitionTrackerAblated,
            SimpleRepetitionTracker,
            UnifiedRecommenderAblated,
            FRFAblation
        )
        
        print("✓ All components imported successfully")
        return True
        
    except Exception as e:
        print(f"✗ Import error: {e}")
        return False


def test_frequency_decoupler():
    """Test frequency domain decoupler with different ablation settings"""
    print("\nTesting Frequency Domain Decoupler...")
    
    try:
        embedding_dim = 64
        batch_size = 4
        
        # Test with frequency domain enabled
        decoupler_enabled = FrequencyDomainDecouplerAblated(
            embedding_dim, enable_freq_domain=True, enable_high_freq_enhancement=True
        )
        
        # Test with frequency domain disabled
        decoupler_disabled = FrequencyDomainDecouplerAblated(
            embedding_dim, enable_freq_domain=False, enable_high_freq_enhancement=False
        )
        
        # Create test input
        test_input = torch.randn(batch_size, embedding_dim)
        
        # Test enabled version
        high_freq_1, low_freq_1, freq_domain_1 = decoupler_enabled(test_input)
        assert high_freq_1.shape == (batch_size, embedding_dim)
        assert low_freq_1.shape == (batch_size, embedding_dim)
        assert freq_domain_1 is not None
        
        # Test disabled version
        high_freq_2, low_freq_2, freq_domain_2 = decoupler_disabled(test_input)
        assert high_freq_2.shape == (batch_size, embedding_dim)
        assert low_freq_2.shape == (batch_size, embedding_dim)
        assert freq_domain_2 is None
        
        # With disabled frequency domain, high and low freq should be similar
        assert torch.allclose(high_freq_2, low_freq_2, atol=1e-6)
        
        print("✓ Frequency Domain Decoupler tests passed")
        return True
        
    except Exception as e:
        print(f"✗ Frequency Domain Decoupler error: {e}")
        return False


def test_repetition_tracker():
    """Test repetition trackers"""
    print("\nTesting Repetition Trackers...")
    
    try:
        embedding_dim = 64
        batch_size = 4
        device = torch.device('cpu')
        
        # Test DTW tracker with DTW enabled
        dtw_tracker_enabled = DTWRepetitionTrackerAblated(
            embedding_dim, enable_dtw=True
        )
        
        # Test DTW tracker with DTW disabled
        dtw_tracker_disabled = DTWRepetitionTrackerAblated(
            embedding_dim, enable_dtw=False
        )
        
        # Test simple tracker
        simple_tracker = SimpleRepetitionTracker(embedding_dim)
        
        # Test all trackers
        for name, tracker in [
            ("DTW Enabled", dtw_tracker_enabled),
            ("DTW Disabled", dtw_tracker_disabled), 
            ("Simple", simple_tracker)
        ]:
            if hasattr(tracker, 'forward'):
                scores = tracker(batch_size, device)
            else:
                scores = tracker.forward(batch_size, device)
            
            assert scores.shape == (batch_size, 1)
            assert torch.all(scores >= 0) and torch.all(scores <= 1)
        
        print("✓ Repetition Tracker tests passed")
        return True
        
    except Exception as e:
        print(f"✗ Repetition Tracker error: {e}")
        return False


def test_unified_recommender():
    """Test unified recommender with different fusion settings"""
    print("\nTesting Unified Recommender...")
    
    try:
        embedding_dim = 64
        batch_size = 4
        
        # Create test inputs
        high_freq = torch.randn(batch_size, embedding_dim)
        low_freq = torch.randn(batch_size, embedding_dim)
        rep_scores = torch.rand(batch_size, 1)
        
        # Test with fusion enabled, dynamic weight
        recommender_fusion = UnifiedRecommenderAblated(
            embedding_dim, enable_fusion=True, fixed_weight=False
        )
        
        # Test with fusion disabled
        recommender_no_fusion = UnifiedRecommenderAblated(
            embedding_dim, enable_fusion=False, fixed_weight=False
        )
        
        # Test with fixed weight
        recommender_fixed = UnifiedRecommenderAblated(
            embedding_dim, enable_fusion=True, fixed_weight=True
        )
        
        # Test all variants
        for name, recommender in [
            ("Fusion Enabled", recommender_fusion),
            ("Fusion Disabled", recommender_no_fusion),
            ("Fixed Weight", recommender_fixed)
        ]:
            output = recommender(high_freq, low_freq, rep_scores)
            assert output.shape == (batch_size, embedding_dim)
        
        print("✓ Unified Recommender tests passed")
        return True
        
    except Exception as e:
        print(f"✗ Unified Recommender error: {e}")
        return False


def test_frf_ablation_model():
    """Test FRF ablation model initialization"""
    print("\nTesting FRF Ablation Model...")
    
    try:
        config = MockConfig()
        dataset = MockDataset()
        
        # Mock additional attributes needed by FRFAblation
        MockSequentialRecommender.n_items = 1000
        MockSequentialRecommender.max_seq_length = 50
        MockSequentialRecommender.ITEM_SEQ = 'item_seq'
        MockSequentialRecommender.ITEM_SEQ_LEN = 'item_seq_len'
        MockSequentialRecommender.POS_ITEM_ID = 'pos_item_id'
        MockSequentialRecommender.NEG_ITEM_ID = 'neg_item_id'
        MockSequentialRecommender.ITEM_ID = 'item_id'
        
        variants = [
            'FRF-Full',
            'FRF-w/o-Freq',
            'FRF-w/o-DTW',
            'FRF-w/o-Fusion',
            'FRF-w/o-HighFreq',
            'FRF-Fixed-Weight'
        ]
        
        for variant in variants:
            try:
                model = FRFAblation(config, dataset, variant=variant)
                
                # Check that the variant configuration is correct
                if variant == 'FRF-w/o-Freq':
                    assert not model.enable_freq_domain
                elif variant == 'FRF-w/o-DTW':
                    assert not model.enable_dtw
                elif variant == 'FRF-w/o-Fusion':
                    assert not model.enable_fusion
                elif variant == 'FRF-w/o-HighFreq':
                    assert not model.enable_high_freq_enhancement
                elif variant == 'FRF-Fixed-Weight':
                    assert model.fixed_weight
                
                print(f"✓ {variant} model initialized successfully")
                
            except Exception as e:
                print(f"✗ {variant} model initialization error: {e}")
                return False
        
        print("✓ All FRF Ablation Model variants tested")
        return True
        
    except Exception as e:
        print(f"✗ FRF Ablation Model error: {e}")
        return False


def test_forward_pass():
    """Test forward pass with mock data"""
    print("\nTesting Forward Pass...")
    
    try:
        config = MockConfig()
        dataset = MockDataset()
        
        # Mock additional attributes
        MockSequentialRecommender.n_items = 100
        MockSequentialRecommender.max_seq_length = 20
        MockSequentialRecommender.ITEM_SEQ = 'item_seq'
        MockSequentialRecommender.ITEM_SEQ_LEN = 'item_seq_len'
        MockSequentialRecommender.POS_ITEM_ID = 'pos_item_id'
        MockSequentialRecommender.NEG_ITEM_ID = 'neg_item_id'
        MockSequentialRecommender.ITEM_ID = 'item_id'
        
        # Create a simple LSTM encoder instead of Transformer
        class MockLSTMEncoder(nn.LSTM):
            def __init__(self):
                super().__init__(64, 64, 1, batch_first=True)
        
        model = FRFAblation(config, dataset, variant='FRF-Full')
        model.encoder = MockLSTMEncoder()  # Replace with simple LSTM
        
        # Create mock input data
        batch_size = 2
        seq_len = 10
        
        item_seq = torch.randint(1, 100, (batch_size, seq_len))
        item_seq_len = torch.tensor([seq_len, seq_len-2])
        
        # Test forward pass
        output = model.forward(item_seq, item_seq_len)
        
        assert output.shape == (batch_size, config['embedding_size'])
        assert not torch.isnan(output).any()
        
        print("✓ Forward pass test passed")
        return True
        
    except Exception as e:
        print(f"✗ Forward pass error: {e}")
        return False


def run_all_tests():
    """Run all tests"""
    print("=" * 60)
    print("FRF Ablation Study Component Tests")
    print("=" * 60)
    print(f"Timestamp: {datetime.now()}")
    print(f"PyTorch Version: {torch.__version__}")
    print("=" * 60)
    
    tests = [
        test_imports,
        test_frequency_decoupler,
        test_repetition_tracker,
        test_unified_recommender,
        test_frf_ablation_model,
        test_forward_pass
    ]
    
    passed = 0
    total = len(tests)
    
    for test in tests:
        try:
            if test():
                passed += 1
        except Exception as e:
            print(f"✗ Test {test.__name__} failed with exception: {e}")
    
    print("\n" + "=" * 60)
    print(f"Test Results: {passed}/{total} tests passed")
    print("=" * 60)
    
    if passed == total:
        print("🎉 All tests passed! Ablation study implementation is working correctly.")
        return True
    else:
        print("❌ Some tests failed. Please check the implementation.")
        return False


if __name__ == '__main__':
    success = run_all_tests()
    sys.exit(0 if success else 1)