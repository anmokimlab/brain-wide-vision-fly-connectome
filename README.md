# brain-wide-vision-fly-connectome

Source code for the manuscript on the connectome analysis of brain-wide vision in *Drosophila*.

This work is available as a preprint on bioRxiv: https://www.biorxiv.org/content/10.64898/2026.02.02.700492v1

## Overview
This repository contains the analysis pipeline used to characterize the long-range projections linking the **optic lobe (OL)** and **central brain (CB)** using the FAFB-FlyWire (v783) dataset, with a cross-dataset replication on the male-CNS (`male-cns:v0.9`) connectome.

Based on synaptic polarity and connectivity, projection neurons are classified into feedforward projection (FFP), feedback projection (FBP), bidirectional projection (BDP), and bilateral projection (BLP) streams, and each stream is then analyzed in turn.

All analysis code lives in [`src/`](src/). The pipeline is numbered end-to-end (`s01`–`s21` data-processing scripts, `fig_*` figure scripts); see [`src/README.md`](src/README.md) for the full step-by-step reproduction guide and the per-figure data dependencies.

## Key Analysis Modules
* **Synaptic Polarity (postPI / prePI):** Per-neuron input/output synapse distributions across the right OL, left OL, and central brain, used to classify neurons into FFP / FBP / BDP / BLP streams.
* **Convergenc/Divergence & Centrality:** Whole-brain graph centrality (PageRank, betweenness) and per-type fan-in/out statistics comparing optic, feedforward, and central neurons.
* **FFP Modular Clustering:** Bipartite Leiden community detection on the FFP→central-brain output matrix, grouping clusters into functional superclusters.
* **FBP Loops:** Input superclass composition, disynaptic optic-lobe input, innervation-depth (layer) profiling against fitted reference layers, and recurrent-loop detection.
* **Segregation & Reciprocity:** Skeleton-based segregation index and synapse-weighted input/output partner reciprocity (weighted Jaccard).
* **Visual Field Mapping:** Receptive fields (RF) and projective fields (PF) of bilateral neurons mapped into visual-space (θ, φ) coordinates via fitted Mi1/Tm3/T4a reference columns.
* **Network Simulation:** A rate-model simulation of the MeMe / Sm / Dm2 / MeTu cell-type network (intact vs. ablated).
* **Cross-dataset Replication (MCNS):** The core analyses repeated on the male-CNS connectome — see [`src/MCNS/`](src/MCNS/).

## Prerequisites
The analysis pipeline is a hybrid of **MATLAB** and **Python**.

### Python Dependencies
* `fafbseg`, `navis` — FlyWire mesh / skeleton access
* `scikit-network` — bipartite Leiden clustering
* `neuprint-python` — male-CNS (MCNS) data export
* `pandas`, `numpy`, `scipy`

### MATLAB Requirements
* Statistics and Machine Learning Toolbox

## Data Access
* **FAFB (main analysis):** Neuron classifications, connections, and synapse coordinates are downloaded from the [FlyWire Codex](https://codex.flywire.ai/?dataset=fafb) (v783) and placed in `src/Codex_Data/`. Neuron and neuropil meshes/skeletons are pulled directly from FlyWire via `fafbseg` / `navis`.
* **MCNS (replication):** The male-CNS (`male-cns:v0.9`) tables are exported from [neuprint](https://neuprint.janelia.org) with `neuprint-python` (an auth token is required) by the `src/MCNS/Data_Processing/s01`–`s04` scripts.
* **Reference data:** The Mi1 column (θ/φ) table from Garner et al. 2024, *Nature* is provided in `src/Reference_Data/`.

Raw Codex downloads and neuprint exports are not committed to the repository; the input folders (`src/Codex_Data/`, `src/MCNS/MCNS_Data/`) are kept as placeholders and populated by following the steps above.

## Usage
The pipeline is numbered and meant to be run in order. In broad strokes:

1. **Get the data** — download the FlyWire Codex files into `src/Codex_Data/` (and, for the replication, run the MCNS export scripts).
2. **Classify the streams** — run `src/Data_Processing/s02_compute_postPI_prePI_all_neurons.m`, then `src/Figures/fig_1D_E_postPI_prePI_right_neurons.m`, which computes the postPI / prePI indices and writes the FFP / FBP / BDP root-id lists that the downstream figures depend on.
3. **Run the per-figure modules** — each main and supplementary figure has a `Data_Processing/s*` step (where needed) followed by a `Figures/fig_*` script, e.g. `s08_FFP_ADJ_matrix.m` → `s09_FFP_leiden_clustering.ipynb` → `fig_3C_D_E_F_G_FFP_clustering_plot.m` for the FFP clustering.

For the complete, figure-by-figure instructions, prerequisites, and the exact inputs/outputs of every script, see **[`src/README.md`](src/README.md)** (and `src/MCNS/` for the male-CNS replication).

## Citation
If you use this code or our findings, please cite:
> Kim, S.Y. and Kim, A.J. "Connectome analysis reveals brain-wide processing of visual features in *Drosophila*." bioRxiv (2026). https://www.biorxiv.org/content/10.64898/2026.02.02.700492v1
