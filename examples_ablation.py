#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2025/1/31
# @Author  : Generated for FRF Ablation Study
# @Description: Example usage of ablation study

"""
Example Usage of FRF Ablation Study
===================================

This script demonstrates how to use the FRF ablation study components
in different scenarios.
"""

import os
import sys

def example_single_variant():
    """Example: Run a single ablation variant"""
    print("Example 1: Running Single Ablation Variant")
    print("-" * 50)
    
    print("Command to run FRF without frequency domain:")
    print("python run_ablation.py --variant FRF-w/o-Freq --dataset ml-1m")
    print()
    
    print("Expected behavior:")
    print("- Frequency domain decoupling will be disabled")
    print("- Model will use standard sequence representation")
    print("- Should show performance degradation compared to full model")
    print()


def example_comprehensive_study():
    """Example: Run comprehensive ablation study"""
    print("Example 2: Comprehensive Ablation Study")
    print("-" * 50)
    
    print("Command to run all ablation variants:")
    print("python run_ablation.py --variant all --dataset ml-1m")
    print()
    
    print("This will run the following variants:")
    variants = [
        ("FRF-Full", "Complete model (baseline)"),
        ("FRF-w/o-Freq", "No frequency domain decoupling"),
        ("FRF-w/o-DTW", "No DTW repetition tracking"),
        ("FRF-w/o-Fusion", "No unified fusion mechanism"),
        ("FRF-w/o-HighFreq", "No high frequency enhancement"),
        ("FRF-Fixed-Weight", "Fixed repetition weights")
    ]
    
    for variant, description in variants:
        print(f"  - {variant}: {description}")
    print()


def example_custom_config():
    """Example: Using custom configuration"""
    print("Example 3: Custom Configuration")
    print("-" * 50)
    
    print("Create a custom config file (custom_ablation.yaml):")
    print("```yaml")
    print("# Custom ablation configuration")
    print("embedding_size: 256      # Larger embeddings")
    print("epochs: 50               # More training epochs")
    print("learning_rate: 0.0005    # Lower learning rate")
    print("cutoff_ratio: 0.7        # Different frequency split")
    print("high_freq_scale: 2.0     # Stronger enhancement")
    print("```")
    print()
    
    print("Run with custom config:")
    print("python run_ablation.py --variant all --config custom_ablation.yaml")
    print()


def example_analysis_workflow():
    """Example: Complete analysis workflow"""
    print("Example 4: Complete Analysis Workflow")
    print("-" * 50)
    
    print("Step 1: Run ablation study")
    print("python run_ablation.py --variant all --dataset ml-1m --output-dir results_ml1m")
    print()
    
    print("Step 2: Check results")
    print("ls results_ml1m/")
    print("# Should show:")
    print("# - ablation_report_ml-1m.txt")
    print("# - ablation_table_ml-1m.tex") 
    print("# - Individual result files for each variant")
    print()
    
    print("Step 3: View comprehensive report")
    print("cat results_ml1m/ablation_report_ml-1m.txt")
    print()
    
    print("Step 4: Include LaTeX table in paper")
    print("# Copy content from ablation_table_ml-1m.tex to your LaTeX document")
    print()


def example_direct_usage():
    """Example: Direct Python usage"""
    print("Example 5: Direct Python Usage")
    print("-" * 50)
    
    print("```python")
    print("from ablation_study import FRFAblation, run_ablation_study")
    print()
    print("# Create configuration")
    print("config = {")
    print("    'embedding_size': 128,")
    print("    'hidden_size': 128,")
    print("    'loss_type': 'BPR',")
    print("    # ... other parameters")
    print("}")
    print()
    print("# Initialize model with specific variant")
    print("model = FRFAblation(config, dataset, variant='FRF-w/o-DTW')")
    print()
    print("# Or run complete ablation study")
    print("result = run_ablation_study('FRF-w/o-Freq', 'ml-1m', 'config.yaml')")
    print("```")
    print()


def example_comparison_analysis():
    """Example: Analyzing component contributions"""
    print("Example 6: Component Contribution Analysis")
    print("-" * 50)
    
    print("After running the ablation study, analyze component contributions:")
    print()
    
    expected_results = [
        ("FRF-Full", "0.1234", "0.0987", "Baseline performance"),
        ("FRF-w/o-Freq", "0.1123", "0.0876", "~10% drop - freq domain important"),
        ("FRF-w/o-DTW", "0.1198", "0.0945", "~3% drop - DTW moderately important"),
        ("FRF-w/o-Fusion", "0.1156", "0.0923", "~6% drop - fusion mechanism valuable"),
        ("FRF-w/o-HighFreq", "0.1189", "0.0934", "~4% drop - enhancement helps"),
        ("FRF-Fixed-Weight", "0.1201", "0.0951", "~3% drop - dynamic weights useful")
    ]
    
    print("Expected analysis (example numbers):")
    print("Variant               | Recall@10 | NDCG@10 | Analysis")
    print("-" * 70)
    for variant, recall, ndcg, analysis in expected_results:
        print(f"{variant:<20} | {recall:<9} | {ndcg:<7} | {analysis}")
    print()
    
    print("Key insights:")
    print("1. Frequency domain decoupling has the largest impact")
    print("2. DTW tracking provides moderate improvements")
    print("3. Fusion mechanism is important for adaptive weighting")
    print("4. High frequency enhancement provides consistent gains")
    print("5. Dynamic weights outperform fixed weights")
    print()


def example_troubleshooting():
    """Example: Common troubleshooting scenarios"""
    print("Example 7: Troubleshooting")
    print("-" * 50)
    
    print("Check dependencies:")
    print("python run_ablation.py --check-deps")
    print()
    
    print("List available datasets:")
    print("python run_ablation.py --list-datasets")
    print()
    
    print("Common issues and solutions:")
    issues = [
        ("CUDA out of memory", "Reduce batch_size in config file"),
        ("RecBole not found", "pip install recbole"),
        ("Dataset not found", "Ensure dataset is in RecBole format"),
        ("Config file error", "Check YAML syntax and required fields"),
        ("Import errors", "Check Python environment and dependencies")
    ]
    
    for issue, solution in issues:
        print(f"  Problem: {issue}")
        print(f"  Solution: {solution}")
        print()


def main():
    """Main function showing all examples"""
    print("=" * 70)
    print("FRF Ablation Study - Usage Examples")
    print("=" * 70)
    print()
    
    examples = [
        example_single_variant,
        example_comprehensive_study,
        example_custom_config,
        example_analysis_workflow,
        example_direct_usage,
        example_comparison_analysis,
        example_troubleshooting
    ]
    
    for i, example_func in enumerate(examples, 1):
        example_func()
        if i < len(examples):
            input("Press Enter to continue to next example...")
            print("\n" + "=" * 70 + "\n")


if __name__ == '__main__':
    main()