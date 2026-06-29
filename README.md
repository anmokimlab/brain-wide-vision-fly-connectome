# brain-wide-vision-fly-connectome

Source code for the manuscript on the connectome analysis for the brain-wide vision in Drosophila.

## Overview
This repository contains the analysis pipeline used to characterize the long-range projections linking the **optic lobe (OL)** and **central brain (CB)** using the FAFB-FlyWire (v783) dataset. 

The code facilitates the classification of projection neurons into feedforward projection (FFP), feedback projection (FBP), and bilateral projection (BLP) streams based on synaptic polarity and connectivity.



## Key Analysis Modules
* **Synaptic Polarity:** Calculation of directionality ($D$) and spatial bias ($B$) indices from synaptic coordinate distributions.
* **FF Modular Clustering:** Modularity-based community detection (Leiden) to identify functional superclusters in the central brain.
* **Recurrent Loop Mapping:** Analysis of layer-specific reciprocity between FBP outputs and disynaptic visual inputs.
* **Visual Field Inference:** Estimation of receptive fields (RF) and projective fields (PF) from synaptic densities.

## Prerequisites
The analysis pipeline is a hybrid of **MATLAB** and **Python**.

### Python Dependencies
* `fafbseg`
* `navis`
* `scikit-network`
* `pandas`, `numpy`, `scipy`

### MATLAB Requirements
* Statistics and Machine Learning Toolbox
  

## Data Access
Neural data and synapse coordinates are accessed via the [FlyWire Codex](https://codex.flywire.ai/). 

## Usage
1. **Classification:** Run `polarity_classification.py` to categorize neurons into FF, FB, and BD groups.
2. **Modularity:** Use `ff_clustering.py` to generate the bipartite graph and community assignments.
3. **Visual Mapping:** Run `bilateral_rf_pf.m` to map synapse coordinates to retinotopic space.

## Citation
If you use this code or our findings, please cite:
> Kim, S.Y. and Kim, A.J. "Connectome analysis reveals brain-wide processing of visual features in Drosophila." bioRxiv (2026): 2026-02.
