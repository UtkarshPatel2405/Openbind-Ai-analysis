# OpenBind AI vs. Physics-Based FEP Benchmarking on EV-A71 2A Protease

This repository contains the benchmarking pipeline and analysis code comparing machine learning co-folding models (Boltz-1, Boltz-2, AF3, OpenFold3-p2, RosettaFold3, Protenix) against traditional docking and Relative Binding Free Energy (RBFE) calculations using the **OpenBind EV-A71 2A protease dataset**.

---

## Objective
To provide a mathematically honest, data-driven comparison of machine learning co-folding models against physics-based simulations on a single dense target containing:
* 925 crystallographic binding events.
* 601 experimental binding affinity ($K_D$) measurements.

---

## Features
* **Bootstrap Confidence Intervals**: 2000 resamples to compute 95% CIs for all correlation metrics.
* **Ensemble Analysis**: Integrates ML and physics scores to capture complementary binding signals.
* **Protocol Sensitivity**: Evaluates all 6 tested Rowan FEP parameters (charge models, starting configurations).

---

## Directory Structure
* `src/setup_workspace.py`: Automates cloning of the OpenBind and Rowan repositories.
* `src/rigorous_analysis.py`: Computes Pearson $r$ and Spearman $\rho$ correlations with bootstrap 95% CIs, p-values, and ensembles.
* `src/generate_figures.py`: Generates the 4-panel publication-quality benchmark figure.
* `results/rigorous_analysis.json`: Exported correlation metrics and bootstrap bounds.
* `results/rigorous_benchmark.png`: The 4-panel benchmark figure.

---

## Installation & Setup
1. Clone this repository:
   ```bash
   git clone https://github.com/yourusername/Openbind-Ai-analysis.git
   cd Openbind-Ai-analysis
   ```
2. Set up the workspace (clones source data):
   ```bash
   python src/setup_workspace.py
   ```
3. Install dependencies:
   ```bash
   pip install pandas scipy matplotlib numpy rdkit
   ```

---

## Running the Pipeline
To run the correlation analysis:
```bash
python src/rigorous_analysis.py
```

To generate the benchmark figures:
```bash
python src/generate_figures.py
```
