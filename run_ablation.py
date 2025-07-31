#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2025/1/31
# @Author  : Generated for FRF Ablation Study
# @Description: Run script for FRF ablation study

"""
FRF Ablation Study Runner
================================================

This script runs comprehensive ablation studies for the FRF model
to evaluate the effectiveness of different components.

Ablation Variants:
- FRF-Full: Complete model with all components
- FRF-w/o-Freq: Remove frequency domain decoupling
- FRF-w/o-DTW: Remove DTW repetition tracker  
- FRF-w/o-Fusion: Remove unified recommender fusion
- FRF-w/o-HighFreq: Remove high frequency enhancement
- FRF-Fixed-Weight: Use fixed repetition weight

Usage Examples:
    # Run single variant
    python run_ablation.py --variant FRF-w/o-Freq --dataset ml-1m
    
    # Run all variants  
    python run_ablation.py --variant all --dataset ml-100k
    
    # Run with custom config
    python run_ablation.py --variant all --config custom_config.yaml
"""

import argparse
import os
import sys
import time
import json
from datetime import datetime
import subprocess

# Add current directory to path for importing ablation_study
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

try:
    from ablation_study import run_ablation_study, run_all_ablation_variants, generate_latex_table
    ABLATION_AVAILABLE = True
except ImportError as e:
    print(f"Warning: Could not import ablation_study module: {e}")
    ABLATION_AVAILABLE = False


def check_dependencies():
    """Check if required dependencies are available"""
    dependencies = ['torch', 'recbole', 'numpy']
    missing = []
    
    for dep in dependencies:
        try:
            __import__(dep)
        except ImportError:
            missing.append(dep)
    
    if missing:
        print(f"Missing dependencies: {missing}")
        print("Please install them using:")
        print(f"pip install {' '.join(missing)}")
        return False
    
    return True


def get_available_datasets():
    """Get list of available datasets"""
    # Common RecBole datasets for sequential recommendation
    return ['ml-1m', 'ml-100k', 'amazon-books', 'yelp', 'gowalla']


def validate_config_file(config_file):
    """Validate configuration file exists and is readable"""
    if not os.path.exists(config_file):
        print(f"Error: Configuration file {config_file} not found")
        return False
    
    try:
        import yaml
        with open(config_file, 'r') as f:
            yaml.safe_load(f)
        return True
    except Exception as e:
        print(f"Error: Invalid YAML configuration file {config_file}: {e}")
        return False


def run_single_variant(variant, dataset, config_file, output_dir):
    """Run single ablation variant"""
    print(f"\n{'='*60}")
    print(f"Running Ablation Study: {variant}")
    print(f"Dataset: {dataset}")
    print(f"Config: {config_file}")
    print(f"{'='*60}")
    
    start_time = time.time()
    
    try:
        if ABLATION_AVAILABLE:
            result = run_ablation_study(variant, dataset, config_file)
        else:
            # Fallback to subprocess if imports fail
            cmd = [
                'python', 'ablation_study.py',
                '--variant', variant,
                '--dataset', dataset,
                '--config', config_file
            ]
            
            print(f"Running command: {' '.join(cmd)}")
            subprocess.run(cmd, check=True)
            result = "Completed via subprocess"
        
        end_time = time.time()
        duration = end_time - start_time
        
        print(f"\n{variant} completed in {duration:.2f} seconds")
        
        # Save individual result
        result_file = os.path.join(output_dir, f"{variant}_{dataset}_result.txt")
        with open(result_file, 'w') as f:
            f.write(f"Ablation Study Result\n")
            f.write(f"Variant: {variant}\n")
            f.write(f"Dataset: {dataset}\n")
            f.write(f"Duration: {duration:.2f} seconds\n")
            f.write(f"Timestamp: {datetime.now()}\n")
            f.write(f"Result: {result}\n")
        
        return result
        
    except Exception as e:
        print(f"Error running {variant}: {str(e)}")
        return None


def run_comprehensive_ablation(dataset, config_file, output_dir):
    """Run comprehensive ablation study with all variants"""
    variants = [
        'FRF-Full',
        'FRF-w/o-Freq', 
        'FRF-w/o-DTW',
        'FRF-w/o-Fusion',
        'FRF-w/o-HighFreq',
        'FRF-Fixed-Weight'
    ]
    
    print(f"\n{'='*80}")
    print("COMPREHENSIVE FRF ABLATION STUDY")
    print(f"{'='*80}")
    print(f"Dataset: {dataset}")
    print(f"Config: {config_file}")
    print(f"Output Directory: {output_dir}")
    print(f"Variants: {len(variants)}")
    print(f"{'='*80}")
    
    results = {}
    total_start_time = time.time()
    
    for i, variant in enumerate(variants, 1):
        print(f"\n[{i}/{len(variants)}] Running {variant}...")
        
        result = run_single_variant(variant, dataset, config_file, output_dir)
        results[variant] = result
        
        # Save intermediate results
        intermediate_file = os.path.join(output_dir, f"intermediate_results_{dataset}.json")
        with open(intermediate_file, 'w') as f:
            json.dump(results, f, indent=2, default=str)
        
        print(f"Intermediate results saved to: {intermediate_file}")
        
        # Add delay between experiments to prevent resource conflicts
        if i < len(variants):
            print("Waiting 5 seconds before next experiment...")
            time.sleep(5)
    
    total_duration = time.time() - total_start_time
    
    # Generate comprehensive report
    report_file = os.path.join(output_dir, f"ablation_report_{dataset}.txt")
    latex_file = os.path.join(output_dir, f"ablation_table_{dataset}.tex")
    
    with open(report_file, 'w') as f:
        f.write("FRF Model Ablation Study - Comprehensive Report\n")
        f.write("=" * 60 + "\n")
        f.write(f"Dataset: {dataset}\n")
        f.write(f"Configuration: {config_file}\n")
        f.write(f"Total Duration: {total_duration:.2f} seconds\n")
        f.write(f"Timestamp: {datetime.now()}\n")
        f.write("=" * 60 + "\n\n")
        
        f.write("ABLATION VARIANTS TESTED:\n")
        f.write("-" * 30 + "\n")
        f.write("1. FRF-Full: Complete model with all components\n")
        f.write("2. FRF-w/o-Freq: Remove frequency domain decoupling\n")
        f.write("3. FRF-w/o-DTW: Remove DTW repetition tracker\n")
        f.write("4. FRF-w/o-Fusion: Remove unified recommender fusion\n")
        f.write("5. FRF-w/o-HighFreq: Remove high frequency enhancement\n")
        f.write("6. FRF-Fixed-Weight: Use fixed repetition weight\n\n")
        
        f.write("RESULTS SUMMARY:\n")
        f.write("-" * 20 + "\n")
        for variant, result in results.items():
            f.write(f"\n{variant}:\n")
            if result:
                f.write(f"  Status: Success\n")
                f.write(f"  Result: {result}\n")
            else:
                f.write(f"  Status: Failed\n")
        
        f.write(f"\n\nDetailed results are available in individual result files.\n")
    
    # Generate LaTeX table if ablation module is available
    if ABLATION_AVAILABLE and any(results.values()):
        try:
            latex_table = generate_latex_table(results)
            with open(latex_file, 'w') as f:
                f.write("% FRF Ablation Study LaTeX Table\n")
                f.write(f"% Generated on: {datetime.now()}\n")
                f.write(f"% Dataset: {dataset}\n\n")
                f.write(latex_table)
            
            print(f"\nLaTeX table saved to: {latex_file}")
        except Exception as e:
            print(f"Warning: Could not generate LaTeX table: {e}")
    
    print(f"\n{'='*80}")
    print("ABLATION STUDY COMPLETED")
    print(f"{'='*80}")
    print(f"Total Duration: {total_duration:.2f} seconds")
    print(f"Report: {report_file}")
    print(f"Results: {output_dir}")
    
    return results


def main():
    """Main function"""
    parser = argparse.ArgumentParser(
        description='FRF Model Ablation Study Runner',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    parser.add_argument('--variant', '-v', type=str, default='all',
                       choices=['all', 'FRF-Full', 'FRF-w/o-Freq', 'FRF-w/o-DTW', 
                               'FRF-w/o-Fusion', 'FRF-w/o-HighFreq', 'FRF-Fixed-Weight'],
                       help='Ablation variant to run (default: all)')
    
    parser.add_argument('--dataset', '-d', type=str, default='ml-1m',
                       choices=get_available_datasets(),
                       help='Dataset to use (default: ml-1m)')
    
    parser.add_argument('--config', '-c', type=str, default='ablation_config.yaml',
                       help='Configuration file (default: ablation_config.yaml)')
    
    parser.add_argument('--output-dir', '-o', type=str, default='ablation_results',
                       help='Output directory for results (default: ablation_results)')
    
    parser.add_argument('--check-deps', action='store_true',
                       help='Check dependencies and exit')
    
    parser.add_argument('--list-datasets', action='store_true',
                       help='List available datasets and exit')
    
    args = parser.parse_args()
    
    # Handle special commands
    if args.check_deps:
        if check_dependencies():
            print("All dependencies are available!")
        else:
            sys.exit(1)
        return
    
    if args.list_datasets:
        print("Available datasets:")
        for ds in get_available_datasets():
            print(f"  - {ds}")
        return
    
    # Validate inputs
    if not validate_config_file(args.config):
        sys.exit(1)
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    print("FRF Model Ablation Study Runner")
    print("=" * 50)
    print(f"Variant: {args.variant}")
    print(f"Dataset: {args.dataset}")
    print(f"Config: {args.config}")
    print(f"Output: {args.output_dir}")
    print("=" * 50)
    
    # Check dependencies
    if not check_dependencies():
        print("\nWarning: Some dependencies are missing.")
        print("The script will attempt to run but may fail.")
        response = input("Continue anyway? (y/N): ")
        if response.lower() != 'y':
            sys.exit(1)
    
    try:
        if args.variant == 'all':
            # Run comprehensive ablation study
            results = run_comprehensive_ablation(args.dataset, args.config, args.output_dir)
        else:
            # Run single variant
            result = run_single_variant(args.variant, args.dataset, args.config, args.output_dir)
            results = {args.variant: result}
        
        print(f"\nAblation study completed successfully!")
        print(f"Results saved in: {args.output_dir}")
        
    except KeyboardInterrupt:
        print("\n\nAblation study interrupted by user.")
        sys.exit(1)
    except Exception as e:
        print(f"\nError during ablation study: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()