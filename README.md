# brain-wide-vision-fly-connectome

Source code for the manuscript on the connectome analysis for the brain-wide vision in Drosophila.

## Overview
This repository contains the analysis pipeline used to characterize the long-range projections linking the **Optic Lobe (OL)** and **Central Brain (CB)** using the FAFB-FlyWire (v783) dataset. 

The code facilitates the classification of projection neurons into Feedforward (FF), Feedback (FB), and Bilateral (BL) streams based on synaptic polarity and connectivity logic.



## Key Analysis Modules
* **Synaptic Polarity:** Calculation of Directionality ($D$) and Spatial Bias ($B$) indices from synaptic coordinate distributions.
* **FF Modular Clustering:** Modularity-based community detection (Leiden) to identify functional superclusters in the CB.
* **Recurrent Loop Mapping:** Analysis of layer-specific reciprocity between FB outputs and disynaptic visual inputs.
* **Visual Field Inference:** Estimation of Receptive Fields (RF) and Projective Fields (PF) from synaptic densities.

## Prerequisites
The pipeline is a hybrid of **MATLAB** and **Python**.

### Python Dependencies
* `fafbseg`
* `navis`
* `scikit-network`
* `pandas`, `numpy`, `scipy`

### MATLAB Requirements
* Statistics and Machine Learning Toolbox
* 

## Data Access
Neural data and synapse coordinates are accessed via the [FlyWire Codex](https://codex.flywire.ai/). Users must provide their own `FlyWire API token` for data retrieval.

## Usage
1. **Classification:** Run `polarity_classification.py` to categorize neurons into FF, FB, and BD groups.
2. **Modularity:** Use `ff_clustering.py` to generate the bipartite graph and community assignments.
3. **Visual Mapping:** Run `bilateral_rf_pf.m` to map synapse coordinates to retinotopic space.

## Citation
If you use this code or our findings, please cite:
> Kim, S.Y. and Kim, A.J. "Connectome analysis reveals brain-wide processing of visual features in Drosophila." bioRxiv (2026): 2026-02.
