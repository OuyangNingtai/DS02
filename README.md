## FreqRepFusion: A Frequency Domain Approach for Recommendation Systems

This repository contains the implementation of FreqRepFusion (FRF), a novel frequency domain approach for sequential recommendation systems that analyzes user behavior patterns in the frequency domain to better capture repetitive interests.

### RecBole Framework

Built on the RecBole framework for recommendation system research.

## Features

- **Frequency Domain Analysis**: Decomposes user interests into high-frequency (repetitive) and low-frequency (exploratory) components
- **DTW Repetition Tracking**: Uses Dynamic Time Warping to track repetitive interest patterns
- **Unified Fusion Mechanism**: Adaptively combines repetitive and exploratory interests
- **Comprehensive Ablation Study**: Evaluates the effectiveness of each component

## Model Variants

### Core Models
- `FRF.py` - Main FreqRepFusion model with enhanced frequency domain processing
- `frf.py` - Simplified FRF implementation
- `FRF02.py` - FRF variant with DTW optimization
- `freqdualdynamic.py` - Frequency dual dynamic model

### Ablation Study
- `ablation_study.py` - Comprehensive ablation study implementation
- `run_ablation.py` - Runner script for ablation experiments
- `ablation.sh` - Convenience script for easy experiment execution

## Quick Start

### Basic Usage

```bash
# Run FRF model
python run_fdd\ \(1\).py --model FRF --dataset ml-1m --config test.yaml

# Run ablation study - single variant
python run_ablation.py --variant FRF-w/o-Freq --dataset ml-1m

# Run complete ablation study
python run_ablation.py --variant all --dataset ml-1m
```

### Using Convenience Script

```bash
# Make script executable
chmod +x ablation.sh

# Quick test
./ablation.sh quick

# Full ablation study
./ablation.sh all ml-1m

# Single variant test
./ablation.sh single FRF-w/o-DTW ml-1m

# View results
./ablation.sh results
```

## Ablation Study

The repository includes a comprehensive ablation study to evaluate component effectiveness:

| Variant | Description |
|---------|-------------|
| **FRF-Full** | Complete model (baseline) |
| **FRF-w/o-Freq** | Remove frequency domain decoupling |
| **FRF-w/o-DTW** | Remove DTW repetition tracker |
| **FRF-w/o-Fusion** | Remove unified recommender fusion |
| **FRF-w/o-HighFreq** | Remove high frequency enhancement |
| **FRF-Fixed-Weight** | Use fixed repetition weight |

See `ABLATION_README.md` for detailed documentation.

## Files Structure

```
├── FRF.py                    # Main FRF model implementation
├── frf.py                    # Simplified FRF model
├── FRF02.py                  # FRF with DTW optimization
├── freqdualdynamic.py        # Frequency dual dynamic model
├── ablation_study.py         # Ablation study implementation
├── run_ablation.py           # Ablation runner script
├── ablation.sh               # Convenience script
├── test.yaml                 # Model configuration
├── ablation_config.yaml      # Ablation study configuration
├── ABLATION_README.md        # Ablation study documentation
└── examples_ablation.py      # Usage examples
```

## Requirements

```bash
pip install torch recbole numpy pyyaml
```

## Citation

If you use this code in your research, please cite our paper:

```bibtex
@article{frf2025,
  title={FreqRepFusion: A Frequency Domain Approach for Recommendation Systems},
  author={...},
  journal={...},
  year={2025}
}
```
