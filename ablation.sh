#!/bin/bash
# -*- coding: utf-8 -*-
# @Time    : 2025/1/31
# @Author  : Generated for FRF Ablation Study
# @Description: Convenience script for running ablation studies

# FRF Ablation Study Convenience Script
# =====================================
# This script provides easy commands for running ablation studies

set -e  # Exit on any error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Default values
DATASET="ml-1m"
CONFIG="ablation_config.yaml"
OUTPUT_DIR="ablation_results"

# Helper function to print colored output
print_status() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

print_header() {
    echo -e "${BLUE}================================${NC}"
    echo -e "${BLUE}$1${NC}"
    echo -e "${BLUE}================================${NC}"
}

# Function to check dependencies
check_deps() {
    print_header "Checking Dependencies"
    python run_ablation.py --check-deps
}

# Function to list available datasets
list_datasets() {
    print_header "Available Datasets"
    python run_ablation.py --list-datasets
}

# Function to run single variant
run_single() {
    local variant=$1
    local dataset=${2:-$DATASET}
    local config=${3:-$CONFIG}
    
    print_header "Running Single Variant: $variant"
    print_status "Dataset: $dataset"
    print_status "Config: $config"
    
    python run_ablation.py --variant "$variant" --dataset "$dataset" --config "$config" --output-dir "$OUTPUT_DIR"
}

# Function to run all variants
run_all() {
    local dataset=${1:-$DATASET}
    local config=${2:-$CONFIG}
    
    print_header "Running Complete Ablation Study"
    print_status "Dataset: $dataset"
    print_status "Config: $config"
    print_status "This will take a while..."
    
    python run_ablation.py --variant all --dataset "$dataset" --config "$config" --output-dir "$OUTPUT_DIR"
}

# Function to run quick test (smaller config)
run_quick() {
    local dataset=${1:-"ml-100k"}  # Use smaller dataset for quick test
    
    print_header "Running Quick Ablation Test"
    print_status "Dataset: $dataset (using smaller dataset for speed)"
    
    # Create quick config
    cat > quick_config.yaml << EOF
# Quick test configuration
embedding_size: 64
hidden_size: 64
n_layers: 1
n_heads: 2
epochs: 5
train_batch_size: 1024
eval_batch_size: 1024
loss_type: 'BPR'
cutoff_ratio: 0.65
rep_weight: 0.75
high_freq_scale: 1.8
window_size: 3
field_separator: "\t"
USER_ID_FIELD: user_id
ITEM_ID_FIELD: item_id
TIME_FIELD: timestamp
load_col:
  inter: [user_id, item_id, timestamp]
MAX_ITEM_LIST_LENGTH: 50
ITEM_LIST_LENGTH_FIELD: item_length
train_neg_sample_args:
  distribution: uniform
  sample_num: 1
eval_setting: TO_LS,full
metrics: ["Recall", "MRR", "NDCG"]
topk: [10]
valid_metric: MRR@10
learner: adam
learning_rate: 0.001
stopping_step: 5
EOF

    python run_ablation.py --variant all --dataset "$dataset" --config quick_config.yaml --output-dir "quick_results"
    
    print_status "Quick test completed. Results in quick_results/"
}

# Function to show results
show_results() {
    local result_dir=${1:-$OUTPUT_DIR}
    
    print_header "Ablation Study Results"
    
    if [ ! -d "$result_dir" ]; then
        print_error "Results directory '$result_dir' not found!"
        return 1
    fi
    
    print_status "Results directory: $result_dir"
    echo ""
    
    # Show comprehensive report if available
    if [ -f "$result_dir/ablation_report_${DATASET}.txt" ]; then
        print_status "Comprehensive Report:"
        echo "=================================="
        cat "$result_dir/ablation_report_${DATASET}.txt"
        echo ""
    fi
    
    # Show LaTeX table if available
    if [ -f "$result_dir/ablation_table_${DATASET}.tex" ]; then
        print_status "LaTeX Table:"
        echo "============"
        cat "$result_dir/ablation_table_${DATASET}.tex"
        echo ""
    fi
    
    # List all result files
    print_status "All result files:"
    ls -la "$result_dir/"
}

# Function to clean results
clean() {
    print_header "Cleaning Results"
    
    if [ -d "$OUTPUT_DIR" ]; then
        print_warning "Removing $OUTPUT_DIR directory..."
        rm -rf "$OUTPUT_DIR"
        print_status "Results cleaned."
    else
        print_status "No results directory to clean."
    fi
    
    # Clean temporary files
    if [ -f "quick_config.yaml" ]; then
        rm quick_config.yaml
        print_status "Temporary config files cleaned."
    fi
}

# Function to validate setup
validate() {
    print_header "Validating Setup"
    
    # Check if required files exist
    required_files=("ablation_study.py" "run_ablation.py" "ablation_config.yaml")
    
    for file in "${required_files[@]}"; do
        if [ -f "$file" ]; then
            print_status "✓ $file found"
        else
            print_error "✗ $file missing"
            return 1
        fi
    done
    
    # Check Python syntax
    print_status "Checking Python syntax..."
    python -m py_compile ablation_study.py
    python -m py_compile run_ablation.py
    print_status "✓ Python syntax OK"
    
    # Check dependencies
    check_deps
    
    print_status "Setup validation completed."
}

# Function to show usage
show_help() {
    echo "FRF Ablation Study - Convenience Script"
    echo "======================================="
    echo ""
    echo "Usage: $0 [COMMAND] [OPTIONS]"
    echo ""
    echo "Commands:"
    echo "  help                  Show this help message"
    echo "  check                 Check dependencies"
    echo "  validate              Validate complete setup"
    echo "  list                  List available datasets"
    echo "  clean                 Clean result files"
    echo ""
    echo "Running ablation studies:"
    echo "  quick [dataset]       Run quick test with small config"
    echo "  all [dataset] [config] Run all ablation variants"
    echo "  single VARIANT [dataset] [config] Run single variant"
    echo ""
    echo "Viewing results:"
    echo "  results [dir]         Show ablation study results"
    echo ""
    echo "Available variants for 'single' command:"
    echo "  FRF-Full              Complete model (baseline)"
    echo "  FRF-w/o-Freq          Without frequency domain decoupling"
    echo "  FRF-w/o-DTW           Without DTW repetition tracker"
    echo "  FRF-w/o-Fusion        Without unified fusion mechanism"
    echo "  FRF-w/o-HighFreq      Without high frequency enhancement"
    echo "  FRF-Fixed-Weight      With fixed repetition weights"
    echo ""
    echo "Examples:"
    echo "  $0 quick                          # Quick test with ml-100k"
    echo "  $0 all ml-1m                      # Full study on MovieLens-1M"
    echo "  $0 single FRF-w/o-Freq ml-1m     # Test frequency ablation"
    echo "  $0 results                        # Show results"
    echo ""
    echo "Default values:"
    echo "  Dataset: $DATASET"
    echo "  Config:  $CONFIG"
    echo "  Output:  $OUTPUT_DIR"
}

# Main script logic
case "$1" in
    "help"|"--help"|"-h"|"")
        show_help
        ;;
    "check")
        check_deps
        ;;
    "validate")
        validate
        ;;
    "list")
        list_datasets
        ;;
    "clean")
        clean
        ;;
    "quick")
        run_quick "$2"
        ;;
    "all")
        run_all "$2" "$3"
        ;;
    "single")
        if [ -z "$2" ]; then
            print_error "Variant required for 'single' command"
            show_help
            exit 1
        fi
        run_single "$2" "$3" "$4"
        ;;
    "results")
        show_results "$2"
        ;;
    *)
        print_error "Unknown command: $1"
        show_help
        exit 1
        ;;
esac