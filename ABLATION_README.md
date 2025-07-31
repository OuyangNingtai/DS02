# FRF Model Ablation Study

This directory contains the implementation of comprehensive ablation studies for the FRF (Frequency Repetitive Fusion) model. The ablation study evaluates the effectiveness of different components in the model architecture.

## Overview

The FRF model consists of several key components:
1. **Frequency Domain Decoupler**: Separates user interests into high-frequency (repetitive) and low-frequency (exploratory) components
2. **DTW Repetition Tracker**: Uses Dynamic Time Warping to track repetitive interest patterns
3. **Unified Recommender**: Fusion mechanism that dynamically combines high-freq and low-freq representations
4. **High Frequency Enhancement**: Amplifies repetitive interest signals
5. **Dynamic Weight Mechanism**: Adaptively adjusts repetition vs. exploration weights

## Ablation Variants

The ablation study includes the following variants:

| Variant | Description |
|---------|-------------|
| **FRF-Full** | Complete model with all components (baseline) |
| **FRF-w/o-Freq** | Remove frequency domain decoupling, use standard sequence representation |
| **FRF-w/o-DTW** | Remove DTW repetition tracker, use simple repetition detection |
| **FRF-w/o-Fusion** | Remove unified recommender fusion, use simple concatenation |
| **FRF-w/o-HighFreq** | Remove high frequency enhancement (set high_freq_scale=1.0) |
| **FRF-Fixed-Weight** | Use fixed repetition weight instead of dynamic weight |

## Files

- `ablation_study.py` - Main ablation study implementation with all model variants
- `run_ablation.py` - Runner script for conducting ablation experiments
- `ablation_config.yaml` - Configuration file optimized for ablation studies
- `ABLATION_README.md` - This documentation file

## Quick Start

### Prerequisites

Make sure you have the required dependencies installed:

```bash
pip install torch recbole numpy pyyaml
```

### Running Single Variant

To run a specific ablation variant:

```bash
python run_ablation.py --variant FRF-w/o-Freq --dataset ml-1m
```

### Running All Variants

To run the complete ablation study with all variants:

```bash
python run_ablation.py --variant all --dataset ml-1m
```

### Custom Configuration

To use a custom configuration file:

```bash
python run_ablation.py --variant all --dataset ml-100k --config custom_config.yaml
```

## Detailed Usage

### Command Line Options

```bash
python run_ablation.py [OPTIONS]

Options:
  --variant, -v     Ablation variant to run (choices: all, FRF-Full, FRF-w/o-Freq, 
                    FRF-w/o-DTW, FRF-w/o-Fusion, FRF-w/o-HighFreq, FRF-Fixed-Weight)
  --dataset, -d     Dataset to use (choices: ml-1m, ml-100k, amazon-books, yelp, gowalla)
  --config, -c      Configuration file (default: ablation_config.yaml)
  --output-dir, -o  Output directory for results (default: ablation_results)
  --check-deps      Check dependencies and exit
  --list-datasets   List available datasets and exit
```

### Direct Python Usage

You can also import and use the ablation study module directly:

```python
from ablation_study import FRFAblation, run_ablation_study

# Run single variant
result = run_ablation_study('FRF-w/o-Freq', 'ml-1m', 'ablation_config.yaml')

# Use the model directly
config = {...}  # Your configuration
dataset = ...   # Your dataset
model = FRFAblation(config, dataset, variant='FRF-w/o-DTW')
```

## Configuration

The `ablation_config.yaml` file contains optimized parameters for ablation studies:

```yaml
# Core model parameters
embedding_size: 128
hidden_size: 128
n_layers: 2
n_heads: 2

# FRF-specific parameters
cutoff_ratio: 0.65              # High-low frequency split ratio
rep_weight: 0.75                # Base repetition weight
high_freq_scale: 1.8            # High frequency enhancement factor
window_size: 7                  # DTW window size

# Training parameters
epochs: 30
train_batch_size: 2048
learning_rate: 0.001
```

Key parameters for ablation:
- `cutoff_ratio`: Controls the frequency domain split point
- `high_freq_scale`: Enhancement factor for high frequency components
- `rep_weight`: Base weight for repetitive interests
- `window_size`: Historical window for DTW analysis

## Output

The ablation study generates several output files:

### Result Files

- `ablation_report_[dataset].txt` - Comprehensive report with all results
- `ablation_table_[dataset].tex` - LaTeX table ready for paper inclusion
- `[variant]_[dataset]_result.txt` - Individual result files for each variant
- `intermediate_results_[dataset].json` - Intermediate results in JSON format

### Example Output Structure

```
ablation_results/
├── ablation_report_ml-1m.txt
├── ablation_table_ml-1m.tex
├── FRF-Full_ml-1m_result.txt
├── FRF-w/o-Freq_ml-1m_result.txt
├── FRF-w/o-DTW_ml-1m_result.txt
├── FRF-w/o-Fusion_ml-1m_result.txt
├── FRF-w/o-HighFreq_ml-1m_result.txt
├── FRF-Fixed-Weight_ml-1m_result.txt
└── intermediate_results_ml-1m.json
```

### LaTeX Table

The generated LaTeX table is ready for inclusion in academic papers:

```latex
\begin{table}[h]
\centering
\caption{Ablation Study Results on FRF Model Components}
\label{tab:ablation_study}
\begin{tabular}{l|ccccc}
\hline
\textbf{Model Variant} & \textbf{Recall@10} & \textbf{NDCG@10} & \textbf{MRR@10} & \textbf{Hit@10} & \textbf{Precision@10} \\
\hline
FRF-Full & 0.1234 & 0.0987 & 0.0765 & 0.2345 & 0.0543 \\
FRF-w/o-Freq & 0.1123 & 0.0876 & 0.0654 & 0.2134 & 0.0432 \\
... 
\hline
\end{tabular}
\end{table}
```

## Performance Analysis

### Expected Results

Based on the model design, we expect:

1. **FRF-Full** should achieve the best overall performance
2. **FRF-w/o-Freq** should show significant performance drop, demonstrating the importance of frequency domain analysis
3. **FRF-w/o-DTW** should have moderate impact on repetition detection accuracy
4. **FRF-w/o-Fusion** should reduce the model's ability to adaptively balance repetition vs. exploration
5. **FRF-w/o-HighFreq** should show reduced performance on repetitive recommendation tasks
6. **FRF-Fixed-Weight** should be less adaptive to different user behaviors

### Metrics

The ablation study evaluates models using standard RecBole metrics:
- **Recall@10/20**: Proportion of relevant items in top-k recommendations
- **NDCG@10/20**: Normalized Discounted Cumulative Gain
- **MRR@10**: Mean Reciprocal Rank
- **Hit@10/20**: Hit rate for top-k recommendations
- **Precision@10/20**: Precision for top-k recommendations

## Datasets

The ablation study supports multiple datasets:

- **MovieLens-1M** (`ml-1m`): 1 million movie ratings
- **MovieLens-100K** (`ml-100k`): 100,000 movie ratings (faster experiments)
- **Amazon Books** (`amazon-books`): Amazon book purchase data
- **Yelp** (`yelp`): Yelp business review data
- **Gowalla** (`gowalla`): Location check-in data

For consistent comparison with existing literature, we recommend using **ml-1m** and **ml-100k**.

## Troubleshooting

### Common Issues

1. **RecBole Import Error**: Make sure RecBole is properly installed
   ```bash
   pip install recbole
   ```

2. **CUDA Memory Error**: Reduce batch size in configuration
   ```yaml
   train_batch_size: 1024
   eval_batch_size: 1024
   ```

3. **Configuration Error**: Validate YAML syntax and required fields

4. **Dataset Not Found**: Ensure dataset is available in RecBole format

### Checking Dependencies

```bash
python run_ablation.py --check-deps
```

### Listing Available Datasets

```bash
python run_ablation.py --list-datasets
```

## Citation

If you use this ablation study implementation in your research, please cite:

```
@article{frf2025,
  title={FreqRepFusion: A Frequency Domain Approach for Recommendation Systems},
  author={...},
  journal={...},
  year={2025}
}
```

## Contributing

To add new ablation variants:

1. Modify the `FRFAblation` class in `ablation_study.py`
2. Add the variant to `configure_ablation_variant()` method
3. Update the choices in `run_ablation.py`
4. Add documentation for the new variant

## License

This code is released under the same license as the main FRF project.