# Connectome Visual Projection Analysis (FAFB)

Code to reproduce the figures from the FAFB connectome (FlyWire Codex) analysis.

## 0. Get the Codex data

Download **all** data files from the FlyWire Codex and place them in the
`Codex_Data` folder:

- Source: **https://codex.flywire.ai/?dataset=fafb**
- Destination: `Codex_Data\`

Some later steps also pull neuron / neuropil meshes and skeletons straight from FlyWire
(the `Data_Processing\s01`, `s03`, `s13`, `s15` notebooks) — just `pip install fafbseg navis`
first and run them; they download what they need.

## Figure 1A — neuron rendering

Two steps: first download the meshes (Python), then render the figure (MATLAB).

1. **Download meshes** — run `Data_Processing\s01_fetch_neuropil_neuron.ipynb`.
   It fetches the meshes from FlyWire and writes them under `Processed_Data\`:
   - whole-brain background mesh → `Processed_Data\fig1a_brain_mesh\`
   - the optic-lobe neuropil region meshes (Me / Lo / LoP / AMe, L and R), plus the
     combined whole-optic-lobe meshes (Optic_R / Optic_L = Me + Lo + LoP + AMe)
     → `Processed_Data\optic_lobe_neuropil_mesh\` (used by `s11_FBP_output_layer.m` and
     `fig_5I_J_LC9_interconnection.m`)
   - the Figure 1A neuron meshes (LT51, VCH, LC14) → `Processed_Data\fig1a_neuron_mesh\`

2. **Render the figure** — run `Figures\fig_1A_neurons.m`.
   It loads the meshes from `Processed_Data\` and `consolidated_cell_types.csv`
   from `Codex_Data\`, colors the neurons by type, and renders the figure.

## Figure 1B — LC9 input/output synapses

1. **Download mesh** — the same notebook `Data_Processing\s01_fetch_neuropil_neuron.ipynb`
   also downloads the Figure 1B LC9 neuron mesh → `Processed_Data\fig1b_neuron_mesh\`.

2. **Render the figure** — run `Figures\fig_1B_LC9.m`.
   It loads the LC9 mesh from `Processed_Data\fig1b_neuron_mesh\` and the synapse
   coordinates (`fafb_v783_princeton_synapse_table.csv`) from `Codex_Data\`, then
   plots the neuron with its pre-/post-synaptic sites and their distribution.

## Figure 1D–E — postPI / prePI of right-hemisphere neurons

1. **Build the NPI table** — run `Data_Processing\s02_compute_postPI_prePI_all_neurons.m`.
   It reads `connections_no_threshold.csv`, `classification.csv`, and
   `consolidated_cell_types.csv` from `Codex_Data\`, computes the per-neuron
   synapse distributions, and saves the table `FAFBNPIs` to
   `Processed_Data\FAFB_NPI_thr0.mat`.

2. **Render the figures** — run `Figures\fig_1D_E_postPI_prePI_right_neurons.m`.
   It loads `Processed_Data\FAFB_NPI_thr0.mat`, computes the postPI / prePI
   indices, and produces three figures:
   - Figure 1: per-type mean postPI vs prePI scatter
   - Figure 2: PostPI + PrePI histogram
   - Figure 3: PostPI − PrePI histogram

   It also classifies the right-hemisphere types into feedforward (FFP),
   feedback (FBP), and bidirectional (BDP, further split into optic / central /
   real) groups, and saves these tables to `Processed_Data\right_neurons_thr0.mat`.
   In addition, the root_ids of the FFP, FBP, BDP (real-bidirectional), and Others
   (other bidirectional = optic + central BD) neurons are written (no header, root_id
   only) to `Processed_Data\right_FFP_root_ids.csv`, `Processed_Data\right_FBP_root_ids.csv`,
   `Processed_Data\right_BDP_root_ids.csv`, and `Processed_Data\right_Others_root_ids.csv`.

## Figure 1H — single neuron at a specific view angle

1. **Download meshes** — the same notebook `Data_Processing\s01_fetch_neuropil_neuron.ipynb`
   also downloads the Figure 1H neuron meshes (LPLC2, cLP01, LT52, LT58, PLP215)
   → `Processed_Data\fig1h_neuron_mesh\`.

2. **Render the figure** — run `Figures\fig_1H_neurons_specific_angle.m`.
   Set `Neuron_root_id` to the neuron you want to plot. It loads that mesh from
   `Processed_Data\fig1h_neuron_mesh\` and the synapse coordinates
   (`fafb_v783_princeton_synapse_table.csv`) from `Codex_Data\`, then renders the
   neuron with its synapses at a fixed view angle, the screen-Y synapse
   distribution, and the projected direction axes.

## Figure 2A — FFP projection neurons

Prerequisite: `Figures\fig_1D_E_postPI_prePI_right_neurons.m` must have been run
first, since it produces the FFP / FBP root_id lists in `Processed_Data\`.

1. **Download meshes** — run `Data_Processing\s03_FFP_FBP_mesh.ipynb`.
   It reads `Processed_Data\right_FFP_root_ids.csv` and
   `Processed_Data\right_FBP_root_ids.csv` and downloads each neuron's mesh into
   `Processed_Data\fig2a_neuron_mesh_FFP\` (FFP, used here) and
   `Processed_Data\fig4a_neuron_mesh_FBP\` (FBP, used by Figure 4A).

2. **Render the figure** — run `Figures\fig_2A_FFP_neurons.m`.
   It draws the background mesh from `Processed_Data\fig1a_brain_mesh\`, then plots
   every FFP neuron in `Processed_Data\fig2a_neuron_mesh_FFP\` colored by type (types
   looked up in `consolidated_cell_types.csv` from `Codex_Data\`), and renders the
   figure.

## Figure 2B — FFP optic-lobe → central-brain connectivity

Prerequisite: `Figures\fig_1D_E_postPI_prePI_right_neurons.m` must have been run
first, since it produces `Processed_Data\right_neurons_thr0.mat` (`RightFFP_NPIs`).

Run `Figures\fig_2B_FFP_neuropil_connection.m`. It loads
`Processed_Data\right_neurons_thr0.mat` and the Codex tables (`neurons.csv`,
`connections_no_threshold.csv`, `consolidated_cell_types.csv`) from `Codex_Data\`,
then produces:

- **Figure 1** — a bubble chart of FFP connections from each right optic-lobe
  region (Me R / Lo R / LoP R) to the top central-brain regions. Each bubble's
  size is the number of FFP neurons for that pair; an extra row and column are
  added on purpose so their bubbles encode the total neuron count per neuropil.
- **Figure 2** — a bar chart of the FFP neuron type distribution (counts entered
  manually). This type-distribution bar chart is the one shown in panel Figure 2A.

It also tabulates the per-neuropil neurotransmitter (NT) composition into the
`NT_<neuropil>` variables (`NT_Me`, `NT_Lo`, `NT_LoP`, `NT_PLP`, ...). Those
counts were copied into Excel to draw the NT pie charts, saved at
`Processed_Data\FFP_neuropil_NT.xlsx`.

## Figure 2C–E — input/output statistics of Optic / FFP / Central neurons

Prerequisite: `Figures\fig_1D_E_postPI_prePI_right_neurons.m` must have been run
first (it produces `Processed_Data\right_neurons_thr0.mat`).

1. **Build the graph centrality table** — run
   `Data_Processing\s04_compute_graph_centrality.m`. It builds the whole-brain
   connectivity graph from `connections_no_threshold.csv` (`Codex_Data\`) and
   computes PageRank and (un)weighted betweenness, saving them to
   `Processed_Data\allgraph_thr0.mat`.

2. **Build the fan-in/out table** — run
   `Data_Processing\s05_compute_FFP_optic_central_fan_in_out.m`. For the right FFP,
   optic, and central neurons it summarizes input/output partners (partner-neuron
   count, partner-type count, synapse count) grouped by partner type, attaches the
   centrality values, aggregates per cell type, and saves
   `Processed_Data\FFP_optic_central_fan_in_out.mat`.

3. **Render the figures** — run `Figures\fig_2C_D_E_plot_FFP_optic_central.m`.
   It loads `Processed_Data\FFP_optic_central_fan_in_out.mat` and draws box plots
   comparing the Optic / Feedforward / Central groups (with Wilcoxon rank-sum
   tests printed to the console): in-neuron number, in-neuron type number,
   in-neuron ratio, out-neuron number, out-neuron type number, out-neuron ratio,
   unweighted betweenness, and PageRank.

## Figure 2F — axo-axonic same-type input of FFP types (optic-lobe proximity-probability reference)

Prerequisite: `Figures\fig_1D_E_postPI_prePI_right_neurons.m` must have been run
first (it produces `Processed_Data\right_neurons_thr0.mat`, `RightFFP_NPIs`). The
optic-lobe proximity step also needs the combined Optic_R mesh from
`Data_Processing\s01_fetch_neuropil_neuron.ipynb`
(`Processed_Data\optic_lobe_neuropil_mesh\`) and the `intriangulation` helper
(`Helper_Function\`) on the MATLAB path.

1. **Compute the data** — run `Data_Processing\s06_FFP_axo_axonic_compute.m`. From
   `connections_no_threshold.csv` and `fafb_v783_princeton_synapse_table.csv`
   (`Codex_Data\`) it computes, for every FFP cell type, the axo-axonic input from
   the same type (the fraction of each neuron's total input synapses arising from
   central-brain connections between neurons of the same type, averaged per type),
   and flags high same-type-input outlier types with the upper Tukey fence
   (Q3 + 1.5 × IQR). For the outlier types only it then computes the optic-lobe
   proximity probability P. Each neuron's OL **input** synapses (the neuron is post)
   are reduced to a single centroid and the OL distance between two neurons is the
   Euclidean distance between their centroids; the same-type CB connection strength W is the
   directional synapse count (presynaptic → partner). P is the per-neuron pooled
   (stratified) probability that a connected partner lies closer in the OL than an
   unconnected one, computed at two thresholds (W > 0, W ≥ 5). Results are saved to
   `Processed_Data\FFP_axo_axonic_outlier_data.mat` (plus `FFP_axo_axonic_perType.csv`
   and `FFP_axo_axonic_outlier_proximity.csv`).

2. **Render the figures** — run `Figures\fig_2F_G_FFP_axo_axonic_plot.m`. It loads
   the saved `.mat` and draws five figures:
   - **Figure 1** (paper panel **2F**) — boxplot of per-type axo-axonic same-type
     input for all FFP types with the upper Tukey fence; this also **defines** the
     high same-type-input outlier types.
   - **Figures 2–5** — the optic-lobe proximity *probability* P readout (per-type P
     bars; per-type P vs axo-axonic input at W > 0 and at W ≥ 5; per-type P boxplot
     against the 0.5 chance line). These are kept only as a **reference** analysis
     and are **not** paper panels. The paper's optic-lobe panels **2G**, **2H**, and
     **S1D** are instead produced by the RF-distance script, which computes the
     **RF distance** (the optic-lobe synapse-centroid distance, µm) — see the next
     section (`Data_Processing\s07_FFP_axo_axonic_RF.m` /
     `Figures\fig_2G_H_FFP_axo_axonic_RF_plot.m`).

## Figure 2F — MeTu3c same-type axo-axonic matrix (example)

The example same-type interconnection matrix shown for **MeTu3c**: among the MeTu3c
neurons of the right hemisphere it visualizes the central-brain (CB) same-type
synapse matrix, ordered by hierarchical clustering of the neurons' connectivity
profiles.

Prerequisite: `Figures\fig_1D_E_postPI_prePI_right_neurons.m` must have been run
first (it produces `Processed_Data\right_neurons_thr0.mat`, `RightFFP_NPIs`).

Run `Figures\fig_2F_MeTu3c_axoaxonic_matrix.m`. From `connections_no_threshold.csv`
(`Codex_Data\`) it selects the MeTu3c neurons, builds the directional same-type
synapse matrix `M(i,j)` (i = pre, j = post) split into an optic-lobe part
(`inter_Optic`) and a central-brain part (`inter_Central`, neuropil that is neither
optic lobe nor UNASGD), and computes a single leaf order by hierarchical clustering
(`pdist` → `linkage('average')` → `optimalleaforder`) of the combined connectivity
profile `[interConnection interConnection']` (`interConnection = inter_Optic +
inter_Central`). It produces:

- **Figure 1** (paper panel **2F**) — the central-brain same-type interconnection
  matrix `inter_Central` reordered by the clustered leaf order, drawn with a reversed
  gray colormap (0 = white, large = black), `clim([0 10])`. The title also reports
  the mean central out-degree (average number of distinct same-type partners a
  neuron sends to).

This is an example-type visualization with no saved data output; change `target_type`
at the top of the script to render the matrix for a different FFP type.

## Figure 2G–H (+ S1D) — RF distance vs same-type connection strength

The same hypothesis as Figure 2F, read out as a **distance** rather than a
proximity *probability*: it computes the **RF distance** — the Euclidean distance
between two neurons' optic-lobe (OL) **input**-synapse centroids (input = neuron is
post) (µm) — as a function of
central-brain (CB) same-type connection strength W. Computed for the high
same-type-input (axo-axonic) outlier types only.

Prerequisite: `Figures\fig_1D_E_postPI_prePI_right_neurons.m` must have been run
first (it produces `Processed_Data\right_neurons_thr0.mat`, `RightFFP_NPIs`). The
proximity step also needs the combined Optic_R mesh from
`Data_Processing\s01_fetch_neuropil_neuron.ipynb`
(`Processed_Data\optic_lobe_neuropil_mesh\`) and the `intriangulation` helper
(`Helper_Function\`) on the MATLAB path.

1. **Compute the data** — run `Data_Processing\s07_FFP_axo_axonic_RF.m`. From
   `connections_no_threshold.csv` and `fafb_v783_princeton_synapse_table.csv`
   (`Codex_Data\`) it flags the high same-type-input outlier types with the upper
   Tukey fence (Q3 + 1.5 × IQR), then for those types computes per directional
   neuron pair the RF distance (OL input-synapse centroid distance, µm) and the
   same-type CB synapse count W.
   Pairs are binned by W (W = 0, 0 < W < 5, W ≥ 5); the per-type value is the mean
   over neurons of each neuron's bin-mean distance. A connected (W > 0) vs
   unconnected (W = 0) split is also saved. Results go to
   `Processed_Data\FFP_axo_axonic_RF_proximity.csv` (plus
   `FFP_axo_axonic_RF_outlier_types.csv`, the fence-flagged types).

2. **Render the figures** — run `Figures\fig_2G_H_FFP_axo_axonic_RF_plot.m`. It loads the
   CSV and draws three figures:
   - **Figure 1** (paper panel **2G**) — per-type 2D scatter of connected (W > 0) RF
     distance (x) vs unconnected (W = 0) RF distance (y), with the y = x reference
     line; points above the line have connected partners closer than unconnected ones.
   - **Figure 2** (paper panel **2H**) — per-type mean RF distance by W bin (boxplot),
     with a Kruskal-Wallis omnibus test and adjacent-bin one-sided ranksum stars
     (Bonferroni-corrected).
   - **Figure 3** (paper panel **S1D**) — connected (W > 0) vs unconnected (W = 0)
     grouped bar chart per type, sorted by D_unconnected / D_connected descending;
     error bar = std. A one-sided signed-rank test of the ratio against 1 is printed.

## Figure 3A — convergence of FFP cell types onto central-brain neurons

Prerequisite: `Figures\fig_1D_E_postPI_prePI_right_neurons.m` must have been run
first, since it produces `Processed_Data\right_neurons_thr0.mat` (`RightFFP_NPIs`).

Run `Figures\fig_3A_FFP_per_CB_neurons.m`. It loads
`Processed_Data\right_neurons_thr0.mat` and the Codex tables
(`connections_princeton.csv`, `consolidated_cell_types.csv`) from `Codex_Data\`.
It keeps the FFP → central-brain connections, maps each presynaptic FFP neuron to
its consolidated cell type, and counts how many distinct FFP cell types converge
onto each central-brain (CB) postsynaptic neuron. It then produces:

- **Figure 1** — histogram of the number of distinct FFP cell types per CB neuron.
- **Figure 2** — bar chart of the fraction of CB neurons receiving 1, 2–10, or
  >10 FFP cell types.

## Figure 3C–G — Leiden clustering of FFP output projections

Prerequisite: `Figures\fig_1D_E_postPI_prePI_right_neurons.m` must have been run
first, since it produces `Processed_Data\right_neurons_thr0.mat` (`RightFFP_NPIs`).

This figure is built in three steps: build the adjacency matrix (MATLAB), run the
Leiden clustering (Python), then render the figures (MATLAB).

1. **Build the FFP output adjacency matrix** — run
   `Data_Processing\s08_FFP_ADJ_matrix.m`. It loads
   `Processed_Data\right_neurons_thr0.mat` (`RightFFP_NPIs`) and
   `connections_princeton.csv` / `classification.csv` from `Codex_Data\`, keeps the
   FFP → central-brain output connections (optic-lobe neuropils dropped, VPN targets
   dropped), and builds the FFP (rows) × post-neuron (cols) synapse matrix
   `WantMatrix_Out`. It writes to `Processed_Data\`:
   - `right_FFP_table.csv` — FFP neuron table (matrix rows)
   - `right_FFP_output_matrix_CB_no_VPN_thr0.csv` / `.mat` — the adjacency matrix
   - `post_right_FFP_CB_no_VPN_thr0.csv` — post-neuron root_ids (matrix columns)
   - `post_neurons_FFP_opticlobe_central_no_VPN_thr0.mat` — same post-neuron ids (`Post_Want`)

2. **Run the Leiden clustering** — run
   `Data_Processing\s09_FFP_leiden_clustering.ipynb` (requires the `scikit-network`
   package). It reads the three CSVs from step 1 and runs bipartite Leiden
   clustering on the adjacency matrix, writing to `Processed_Data\`:
   - `leiden_right_FFP_output_CB_no_VPN_thr0_resol3_100000.csv` — FFP-row cluster assignments
   - `leiden_post_right_FFP_CB_no_VPN_thr0_100000.csv` — post-column cluster assignments
   - `leiden_right_FFP_intercluster_connectivity_CB_no_VPN_thr0_resol3_100000.csv` — inter-cluster connectivity
   - `dendrogram_raw_CB_no_VPN_thr0_resol3_100000.csv` — raw-similarity dendrogram (linkage matrix)

3. **Render the figures** — run `Figures\fig_3C_D_E_F_G_FFP_clustering_plot.m`. It
   loads the s08/s09 outputs above plus the Codex tables
   (`connections_no_threshold.csv`, `consolidated_cell_types.csv`,
   `classification.csv`, `fafb_v783_princeton_synapse_table.csv`) from `Codex_Data\`.
   The 32 Leiden clusters are grouped into 7 hand-defined **superclusters**
   (`SuperClusterOrder`). It produces:
   - **Figure 1** — FFP-cluster dendrogram, reordered by leaf order.
   - **Figure 2** — inter-cluster connectivity heatmap, reordered by leaf order.
   - **Figure 3** — cluster-level projection similarity (weighted Jaccard), within- vs
     between-cluster horizontal box plot + Wilcoxon rank-sum test.
   - **Figure 4** — same, at the supercluster level.
   - **Figure 5** — cluster-level output-synapse distance, within- vs between-cluster
     box plot + Wilcoxon.
   - **Figure 6** — same, at the supercluster level.
   - **Figure 7** — per-type clustering consistency vs type size (scatter) with the
     consistency CDF on a shared X-axis.

   The per-supercluster output neuropil targets are computed in the variable
   `SuperCluster_neuropil_targets` (a 7×1 cell array; for each supercluster, column 1
   is the neuropil name and column 2 is the total output synapse count, sorted
   descending). It is computed but not plotted here. The script also tabulates, per
   supercluster, the post-synaptic neuron-type counts, neuropil synapse counts, and
   per-type synapse counts into the struct `SuperCluster_targets`, saved to
   `Processed_Data\FFP_supercluster_targets.mat` for Figure S1E–J.

## Figure 4A — FBP projection neurons

Prerequisite: `Figures\fig_1D_E_postPI_prePI_right_neurons.m` must have been run
first, since it produces the FFP / FBP root_id lists in `Processed_Data\`.

1. **Download meshes** — run `Data_Processing\s03_FFP_FBP_mesh.ipynb` (the same
   notebook used for Figure 2A). It downloads the FBP neuron meshes into
   `Processed_Data\fig4a_neuron_mesh_FBP\`.

2. **Render the figure** — run `Figures\fig_4A_FBP_neurons.m`.
   It draws the background mesh from `Processed_Data\fig1a_brain_mesh\`, then plots
   every FBP neuron in `Processed_Data\fig4a_neuron_mesh_FBP\` colored by type (types
   looked up in `consolidated_cell_types.csv` from `Codex_Data\`), and renders
   the figure.

## Figure 4B — FBP optic-lobe ↔ central-brain connectivity

Prerequisite: `Figures\fig_1D_E_postPI_prePI_right_neurons.m` must have been run
first, since it produces `Processed_Data\right_neurons_thr0.mat` (`RightFBP_NPIs`).

Run `Figures\fig_4B_FBP_neuropil_connection.m`. It loads
`Processed_Data\right_neurons_thr0.mat` and the Codex tables (`neurons.csv`,
`connections_no_threshold.csv`, `consolidated_cell_types.csv`) from `Codex_Data\`,
then produces:

- **Figure 1** — a bubble chart of FBP connections linking each right optic-lobe
  region (Me R / Lo R / LoP R) to the top central-brain regions. Each bubble's size
  is the number of FBP neurons for that pair; an extra row and column are added on
  purpose so their bubbles encode the total neuron count per neuropil.
- **Figure 2** — a bar chart of the FBP neuron type distribution (counts entered
  manually). This type-distribution bar chart is the one shown in panel Figure 4A.

It also tabulates the per-neuropil neurotransmitter (NT) composition into the
`NT_<neuropil>` variables (`NT_Me`, `NT_Lo`, `NT_LoP`, `NT_PLP`, ...). Those counts
were copied into Excel to draw the NT pie charts, saved at
`Processed_Data\FBP_neuropil_NT.xlsx`.

## FBP output-neuropil classification

Run `Data_Processing\s10_FBP_out_opticlobes.m`. It loads
`Processed_Data\right_neurons_thr0.mat` (`RightFBP_by_type`) and the Codex tables
(`consolidated_cell_types.csv`, `connections_no_threshold.csv`) from `Codex_Data\`,
and, for every FBP cell type, computes the fraction (%) of its output synapses
landing in each optic-lobe region, collapsing left/right into AME / ME / LO / LOP
(stored in the variable `RightFBP_OutNeuropils`). These values were written to
`Processed_Data\FBP_output_neuropils.xlsx`, which was then used to classify each FBP
neuron type by the optic-lobe region it targets.

(The local function `seeConnection_by_region`, which summarizes each type's
input/output partners and neuropils, is defined as a local function at the bottom of
the script.)

## Figure 4C (+ S2A) — superclass composition of FBP inputs

Prerequisite: `Figures\fig_1D_E_postPI_prePI_right_neurons.m` must have been run
first, since it produces `Processed_Data\right_neurons_thr0.mat` (`RightFBP_NPIs`).

Run `Figures\fig_4C_FBP_input_superclass.m`. It loads
`Processed_Data\right_neurons_thr0.mat` and the Codex tables
(`connections_no_threshold.csv`, `classification.csv`, `consolidated_cell_types.csv`)
from `Codex_Data\`. For every FBP cell type it aggregates the **central-brain input
neurons** (optic-lobe inputs excluded) by FlyWire `super_class`. It also builds a
"central-reassigned" version (`input_superclass2`) in which each type's `central`
input fraction is redistributed to the superclasses of those central neurons' own
top-5 inputs. FBP types are grouped into **Me / Lo / Lop / Multi** by their target
optic lobe (hardcoded index lists from the FBP output-neuropil classification, s10).
It produces:

- **Figure 1** (paper panel **S2A**) — per-type input superclass composition
  (one stacked bar per FBP type).
- **Figure 2** (paper panel **4C**) — per-group (Me / Lo / Lop / Multi) averaged
  input superclass composition.
- **Figure 3** (paper panel **S2A**) — per-type composition after the central
  reassignment.
- **Figure 4** (paper panel **4C**) — per-group composition after the central
  reassignment.

Two further stacked bars (drawn by the local function `plotStackedSelectedBar`) show
only the selected superclasses {ascending, central, descending, visual_centrifugal,
visual_projection} plus an "Others" bucket (a per-type and a central-reassigned version).
The output-side analysis is not used and has been removed.

## Figure 4D (+ S2B) — disynaptic optic-lobe input to FBP neurons

Prerequisite: `Figures\fig_1D_E_postPI_prePI_right_neurons.m` must have been run
first, since it produces `Processed_Data\right_neurons_thr0.mat` (`RightFBP_NPIs`).

Run `Figures\fig_4D_FBP_dysynaptic_optclobe.m`. It loads
`Processed_Data\right_neurons_thr0.mat` and the Codex tables
(`connections_no_threshold.csv`, `consolidated_cell_types.csv`) from `Codex_Data\`.
For every FBP cell type it estimates the **disynaptic optic-lobe input** along the
path `optic lobe -> upstream (central, non-optic) neuron -> FBP neuron`: it finds
each FBP type's upstream non-optic partners, then measures how much optic-lobe input
those upstream neurons receive, weighted by their synapse contribution, broken down
into six regions (Me R, Me L, Lo R, Lo L, LoP R, LoP L). FBP types are grouped into
**Me / Lo / Lop / Multi** by their target optic lobe (hardcoded index lists from the
FBP output-neuropil classification, s10). It produces:

- **Figure 1** (paper panel **4D**) — `imagesc` of the mean ± std disynaptic
  optic-lobe input fraction per neuropil (rows) × FBP group (columns: CB-Me / CB-LO /
  CB-LOP / CB-Multi), with the mean ± std annotated in each cell.
- **Figure 2** (paper panel **S2B**) — per Me-targeting FBP type, the disynaptic Me
  input split into contralateral (Me L) vs ipsilateral (Me R).
- **Figure 3** (paper panel **S2B**) — the same for Lo-targeting types (Lo L vs Lo R).
- **Figure 4** (paper panel **S2B**) — the same for Lop-targeting types (LoP L vs LoP R).

The upstream-partner helper `seeConnection_root_id_NoOptic` is defined as a local
function at the bottom of the script.

## Figure 4E (+ S2C, S2D) — FBP output and disynaptic-input depth layers

Prerequisite: `Figures\fig_1D_E_postPI_prePI_right_neurons.m` must have been run first
(it produces `Processed_Data\right_neurons_thr0.mat`, `RightFBP_NPIs`). Then build the
two synapse-layer `.mat` files by running, in order:

1. **`Data_Processing\s11_FBP_output_layer.m`** — for each optic lobe (Me / Lo / LoP)
   it fits a **reference layer** (PCA + parabolic surface; Me = `Dm6` dendrite,
   Lo = `LT1a` dendrite, LoP = `T5a` axon) and measures the depth of the FBP group's
   **own output (axon) synapses** relative to that layer, for the right (ipsilateral)
   and left (contralateral) hemispheres. Its reference-layer PCA figures are the
   **output-side S2C** panels. It saves the six matrices `FBP_Synapse_{ME,LO,LOP}_{R,L}`
   to `Processed_Data\FBP_output_synapse_layer.mat`.

2. **`Data_Processing\s12_FBP_upstream_layer.m`** — the upstream counterpart: same
   reference layers, but it profiles the in-neuropil synapses of each FBP type's
   **upstream (central, non-optic) partners** (`seeConnection_root_id_NoOptic`),
   weighted by each partner's synapse contribution. Its reference-layer PCA figures are
   the **upstream-side S2C** panels. It saves the six matrices
   `FBP_Upstream_Synapse_{ME,LO,LOP}_{R,L}` to
   `Processed_Data\FBP_upstream_synapse_layer.mat`.

Both scripts load the optic-lobe neuropil meshes from
`Processed_Data\optic_lobe_neuropil_mesh\` and need the `intriangulation` helper
(`Helper_Function\`) on the MATLAB path. They group FBP types into Me / Lo / Lop by
their target optic lobe (hardcoded index lists from the FBP output-neuropil
classification, s10).

Finally, run `Figures\fig_4E_FBP_layer.m`. It loads both `.mat` files and, for each
optic lobe, renders one RGB image whose columns are per-FBP-type triplets — upstream
input ipsi (blue) · FBP output (orange) · upstream input contra (green) — stacked over
the innervation-depth axis. The three optic-lobe images (ME / LO / LOP) are the panels
**4E** and **S2D**.

## Figure 4F — overall depth-profile histograms

Prerequisite: the two synapse-layer `.mat` files from `s11_FBP_output_layer.m` and
`s12_FBP_upstream_layer.m` (see above) must already exist in `Processed_Data\`.

Run `Figures\fig_4F_FBP_layer_overall_histogram.m`. It loads
`Processed_Data\FBP_output_synapse_layer.mat` and
`Processed_Data\FBP_upstream_synapse_layer.mat`, replaces missing bins with 0
(`fillmissing`), and for each optic lobe (ME / LO / LOP) collapses all FBP types into a
single mean ± 1 SD profile along the innervation-depth axis, overlaying three traces:
FBP output (orange), upstream input ipsilateral (blue), and upstream input
contralateral (green). The three plots are panel **4F**. The mean ± SD band is drawn by the local
function `shaded_band` at the bottom of the script.

## Figure 4G (+ S2E) — input/output depth overlap

Prerequisite: the two synapse-layer `.mat` files from `s11_FBP_output_layer.m` and
`s12_FBP_upstream_layer.m` (see above) must already exist in `Processed_Data\`.

Run `Figures\fig_4G_FBP_layer_overlap.m`. It loads
`Processed_Data\FBP_output_synapse_layer.mat`,
`Processed_Data\FBP_upstream_synapse_layer.mat`, and
`Processed_Data\right_neurons_thr0.mat`. For each FBP type it Gaussian-smooths the
output depth profile (`X`) and the upstream-input depth profiles (`Y_L`, `Y_R`, each
scaled by its share of the total upstream synapses), then computes the **overlap
coefficient** between them as the per-depth `sum(min(...))`, separately for the
contralateral (L) and ipsilateral (R) input. It produces four figures:

- **Figures 1, 2, 3** (paper panel **S2E**) — per-neuron overlap bar charts for ME / LO /
  LOP (green = contralateral input, blue = ipsilateral input).
- **Figure 4** (paper panel **4G**) — a notched box plot of the overlap distribution
  across all neurons, grouped by region and side (Me L/R, Lo L/R, LoP L/R).

The overall mean ± std of the overlap (across all neurons, and per side) is also printed
to the console.

## Figure 4H — loop-forming FBP neuron percentage

Prerequisite: `Figures\fig_1D_E_postPI_prePI_right_neurons.m` must have been run first (it
produces `Processed_Data\right_neurons_thr0.mat`, `RightFBP_NPIs`).

Run `Figures\fig_4H_FBP_loop.m`. It loads `connections_no_threshold.csv` and
`consolidated_cell_types.csv` from `Codex_Data\` and `RightFBP_NPIs` from
`Processed_Data\right_neurons_thr0.mat`, and builds the binary neuron-level optic-lobe (OL)
and central-brain (CB) adjacency. For every right-hemisphere FBP (feedback) neuron `a` it
takes its **disynaptic-visual-input** partners — the CB input neurons `b` (b → a in CB) that
themselves receive OL input — and tests whether `a` closes a feedback loop back onto any of
them through its OL output:

- **direct**   : `a -> b` (one OL output hop) ⇒ the cycle `a -> b -> a`
- **indirect** : `a -> x -> b` (two OL hops, x ≠ b) ⇒ the cycle `a -> x -> b -> a`

A focal neuron "forms a loop" (direct / indirect) if it reaches ≥ 1 such partner. FBP types
are grouped into Me / Lo / Lop / Multi by their target optic lobe (hardcoded index lists from
the FBP output-neuropil classification, s10). It produces:

- **Figure 1** (paper panel **4H**) — per FB neuropil subgroup (Me / Lo / Lop / Multi), the
  **percentage of that subgroup's focal FBP neurons that form a direct / indirect feedback
  loop** (each focal neuron counted once; the per-subgroup percentages are also printed to the
  console).

The `subgroupLoopCounts` helper is a local function at the bottom of the script.

## Figure 5A — example neuron morphologies along their synapse axis

1. **Download meshes** — the same notebook `Data_Processing\s01_fetch_neuropil_neuron.ipynb`
   also downloads the Figure 5A neuron meshes (LC9, LT43, LT52)
   → `Processed_Data\fig5a_neuron_mesh\`.

2. **Render the figures** — run `Figures\fig_5A_BDP_neurons.m`.
   It loads the synapse coordinates (`fafb_v783_princeton_synapse_table.csv`) from
   `Codex_Data\` and, for each of the three neurons, aligns the mesh and its
   pre-/post-synaptic sites to the synapse-cloud PCA axes, then renders the neuron
   (gray) with its input (post) and output (pre) synapses at a per-neuron view angle
   (`view([az el])` stored alongside each neuron). The three neurons are drawn in
   **Figures 1, 2, 3** (LC9 / LT43 / LT52).

## Figure 5B–D — bias direction and PostPI / PrePI of bidirectional neurons

Prerequisite: `Figures\fig_1D_E_postPI_prePI_right_neurons.m` must have been run first,
since it produces `Processed_Data\right_neurons_thr0.mat` (the `RightFFP_*`, `RightFBP_*`,
and bidirectional `RightBDP_real_*` tables).

Run `Figures\fig_5B_C_D_BDP_postPI_prePI.m`. It loads
`Processed_Data\right_neurons_thr0.mat` and compares the feedforward (FFP), feedback
(FBP), and real-bidirectional (BDP) groups, highlighting three example BDP neuron types
(LC9, LT43, LT52). Set `SeeType` to 1 to aggregate per cell type or 0 to use per-neuron
values (the LC9 / LT43 / LT52 examples are always per neuron). It produces:

- **Figure 1** (paper panel **5B**) — bias-direction scatter of each group's mean
  `Post+Pre` (B) vs `Post-Pre` (D) with SD error bars, plus the three example types.
- **Figure 2** (paper panel **5C**) — PostPI box plot across LC9 / LT43 / LT52 / BDP /
  FFP / FBP.
- **Figure 3** (paper panel **5D**) — the same for PrePI.

The BDP / FFP / FBP pairwise Wilcoxon rank-sum tests (Bonferroni-corrected) and
probabilistic-dominance values for PostPI and PrePI are printed to the console. The
`ternary` helper is a local function at the bottom of the script.

## Figure 5E — axon/dendrite segregation index

Prerequisite: `Figures\fig_1D_E_postPI_prePI_right_neurons.m` must have been run first,
since it writes the per-group root_id lists (`right_FFP_root_ids.csv`,
`right_FBP_root_ids.csv`, `right_BDP_root_ids.csv`, `right_Others_root_ids.csv`) to
`Processed_Data\`.

1. **Compute the segregation index** — run `Data_Processing\s13_segregation_index.ipynb`.
   For each group (FFP / FBP / BDP / Others) it reads the root_id list, skeletonizes each
   neuron from FlyWire, splits it into axon/dendrite, and computes its segregation index.
   It writes `Processed_Data\right_<group>_segregation_index.csv` (columns: `root_id`,
   `segregation_index`) for each group. (Figure 5E uses FFP / FBP / BDP; the Others group
   is used by Figure 5F.)

2. **Render the figure** — run `Figures\fig_5E_segregation_index_plot.m`. It loads
   `Processed_Data\right_neurons_thr0.mat` (the `RightFFP_NPIs`, `RightFBP_NPIs`,
   `RightBDP_real_NPIs` tables), attaches the segregation index from the three CSVs by
   `root_id`, and averages it per cell type. It produces:

   - **Figure 1** (paper panel **5E**) — segregation-index box plot across LC9 / LT43 /
     LT52 / BDP / FFP / FBP.

   The BDP / FFP / FBP pairwise Wilcoxon rank-sum tests (Bonferroni-corrected) and
   probabilistic-dominance values are printed to the console. The
   `attach_segregation_index`, `summarize_by_type`, and `ternary` helpers are local
   functions at the bottom of the script.

## Figure 5F–G — input/output reciprocity vs segregation index

Prerequisites: `Figures\fig_1D_E_postPI_prePI_right_neurons.m` (produces
`Processed_Data\right_neurons_thr0.mat`) and `Data_Processing\s13_segregation_index.ipynb`
(produces the `right_<group>_segregation_index.csv` files) must have been run first.

1. **Build the reciprocity table** — run `Data_Processing\s14_BDP_reciprocal.m`. For every
   right-hemisphere FFP, FBP, real-bidirectional (BDP), and Others (other bidirectional =
   optic + central BD) neuron it measures how much the input-partner set overlaps with the
   output-partner set, using the synapse-weighted Jaccard index (Σmin / Σmax over the union
   of partner neuron ids). It loads `connections_no_threshold.csv` from `Codex_Data\` and
   `right_neurons_thr0.mat`, averages the per-neuron values per cell type, keeps the three
   example BDP types (LC9 / LT43 / LT52) per neuron in `Want`, and saves `Type_FFP`,
   `Type_FBP`, `Type_BDP`, `Type_Others`, `Want` to `Processed_Data\BDP_reciprocity.mat`.
   The `reciprocity_for_neurons`, `weighted_jaccard`, and `summarize_reciprocity` helpers
   are local functions at the bottom of the script.

2. **Render the figures** — run `Figures\fig_5F_G_segregation_reciprocal.m`. It loads
   `Processed_Data\BDP_reciprocity.mat` and attaches the per-type mean segregation index
   from the `right_<group>_segregation_index.csv` files (by `root_id`). It produces:

   - **Figure 1** (paper panel **5F**) — per-type segregation index (x) vs reciprocity (y)
     scatter for FFP / FBP / BDP / Others, with a linear fit + 95% prediction interval (fit
     statistics annotated, fit over all four groups) and the three example BDP types
     highlighted.
   - **Figure 2** (paper panel **5G**) — reciprocity box plot across LC9 / LT43 / LT52 /
     BDP / FFP / FBP.

   The BDP / FFP / FBP pairwise Wilcoxon rank-sum tests (Bonferroni-corrected) and
   probabilistic-dominance values are printed to the console. The `attach_mean_segregation`
   and `ternary` helpers are local functions at the bottom of the script.

## Figure 5I–J — LC9–LC9 interconnection vs proximity

Prerequisites: `Figures\fig_1D_E_postPI_prePI_right_neurons.m` must have been run first
(it produces `Processed_Data\right_neurons_thr0.mat`, used here for the right-hemisphere
LC9 neurons), and `Data_Processing\s01_fetch_neuropil_neuron.ipynb` must have produced the
combined right-optic-lobe mesh `Processed_Data\optic_lobe_neuropil_mesh\Optic_R_vertices.csv`
/ `Optic_R_faces.csv`. The `intriangulation` helper (`Helper_Function\`) must be on the
MATLAB path.

Run `Figures\fig_5I_J_LC9_interconnection.m`. It loads `connections_no_threshold.csv` and
`fafb_v783_princeton_synapse_table.csv` from `Codex_Data\` and the LC9 root_ids from
`right_neurons_thr0.mat`. It builds the LC9→LC9 central-brain interconnection matrix
(synapse counts), and for each LC9 neuron averages its input synapses (the neuron is post) inside the
right optic lobe (kept via `intriangulation` against the Optic_R mesh) into a single OL
centroid. Each pair's distance is then the Euclidean distance between the two neurons'
OL synapse centroids, and connection weight is related to proximity over all neuron pairs.
It produces:

- **Figure 1** (paper panel **5I**) — box plot of inter-centroid distance for pairs that
  are strongly connected (W ≥ 5), connected (W > 0), or unconnected (W = 0), with the
  connected-vs-unconnected Wilcoxon rank-sum p-value in the title.
- **Figure 2** (paper panel **5J**) — connection probability vs distance (grouped bars for
  W > 0 and W ≥ 5 on the left axis) overlaid with the per-bin pair count on the right axis.

## Figure 6B — bilateral (BLP) PostPI / PrePI of bilateral neurons

Prerequisite: `Data_Processing\s02_compute_postPI_prePI_all_neurons.m` must have been run
first (it produces `Processed_Data\FAFB_NPI_thr0.mat`).

Run `Figures\fig_6B_postPI_prePI_BLP_neurons.m`. It loads
`Processed_Data\FAFB_NPI_thr0.mat` and computes, per neuron, the **BLP** (bilateral) PostPI
/ PrePI (right vs left optic-lobe lateralization of input / output synapses). It selects
BLP neurons in three steps — (1) keep neurons with ≥ 5 synapses bridging the right and left
optic lobes, (2) drop neurons whose input or output is mostly central (> 70%), (3) drop
cell types kept in fewer than 20% of their neurons (and, where a type has a 'right' neuron,
its 'left' neurons are removed first) — then aggregates per cell type. It produces:

- **Figure 1** (paper panel **6B**) — per-type mean BLP PostPI vs PrePI scatter, colored by
  superclass.

It then classifies the BLP types into right (R), left (L), and bidirectional (BD) groups by
`Mean_BLP_PostPI − Mean_BLP_PrePI` and saves `BLP_R_by_type`, `BLP_L_by_type`,
`BLP_BD_type`, `BLP_R_NPIs`, `BLP_L_NPIs`, `BLP_BD_NPIs` to
`Processed_Data\BLP_neurons_thr0.mat`. In addition, the root_ids of the BLP_R neurons are
written (no header, root_id only) to `Processed_Data\BLP_R_root_ids.csv`.

## Figure 6C — BLP_R projection neurons

Prerequisite: `Figures\fig_6B_postPI_prePI_BLP_neurons.m` must have been run first, since it
produces the BLP_R root_id list `Processed_Data\BLP_R_root_ids.csv`.

1. **Download meshes** — run `Data_Processing\s15_BLP_mesh.ipynb`. It reads
   `Processed_Data\BLP_R_root_ids.csv` and downloads each neuron's mesh into
   `Processed_Data\fig6c_neuron_mesh_BLP`.

2. **Render the figure** — run `Figures\fig_6C_BLP_neurons.m`. It draws the background mesh
   from `Processed_Data\fig1a_brain_mesh\`, then plots every BLP_R neuron in
   `Processed_Data\fig6c_neuron_mesh_BLP` colored by type (types looked up in
   `consolidated_cell_types.csv` from `Codex_Data\`), and exports the image (paper panel
   **6C**).

## Figure 6D — BLP_R optic-lobe ↔ optic-lobe connectivity

Prerequisite: `Figures\fig_6B_postPI_prePI_BLP_neurons.m` must have been run first, since it
produces `Processed_Data\BLP_neurons_thr0.mat` (`BLP_R_NPIs`).

Run `Figures\fig_6D_BLP_neuropil_connection.m`. It loads
`Processed_Data\BLP_neurons_thr0.mat` and the Codex tables (`neurons.csv`,
`connections_no_threshold.csv`, `consolidated_cell_types.csv`) from `Codex_Data\`, then
produces:

- **Figure 1** (paper panel **6D**) — a bubble chart of BLP_R connections linking each right
  optic-lobe region (Me R / Lo R / LoP R) to each left optic-lobe region (Me L / Lo L /
  LoP L). Each bubble's size is the number of neurons for that pair; an extra row and column
  are added on purpose so their bubbles encode the total neuron count per neuropil.
- **Figure 2** (paper panel **6C**) — a bar chart of the BLP_R neuron type distribution
  (counts entered manually). This type-distribution bar chart is the one shown in panel
  Figure 6C.

It also tabulates the per-neuropil neurotransmitter (NT) composition into the
`NT_<neuropil>` variables (`NT_Me_R`, `NT_Lo_R`, `NT_LoP_R`, `NT_Me_L`, `NT_Lo_L`,
`NT_LoP_L`). Those counts were copied into Excel to draw the NT pie charts, saved at
`Processed_Data\BLP_neuropil_NT.xlsx`.

## Figure 6E–G (+ S4) — BLP_R projective / receptive fields in visual space

This figure maps each BLP_R neuron's input synapses (receptive field, RF) and output
synapses (projective field, PF) into visual-space (theta, phi) coordinates. It is built in
three steps: build the reference columns, fit the reference-layer surfaces, then render the
per-neuron-type fields.

Prerequisites: `Figures\fig_6B_postPI_prePI_BLP_neurons.m` must have been run first (it
produces `Processed_Data\BLP_neurons_thr0.mat`, `BLP_R_NPIs` / `BLP_R_by_type`), and `s01`
must have produced the optic-lobe neuropil meshes in
`Processed_Data\optic_lobe_neuropil_mesh\`. The `intriangulation` helper
(`Helper_Function\`) must be on the MATLAB path. One further input you must provide:

- **Mi1 column reference** — `Reference_Data\Mi1_columns_DustinGarner.csv`, the Mi1
  column (theta/phi) table from Garner et al. 2024, *Nature*
  (https://www.nature.com/articles/s41586-024-07967-z).

The neuron SWC skeletons (`Codex_Data\sk_lod1_783_healed\`) are **not** needed and do not have
to be downloaded: they were only used to derive skeleton voxels and the output-synapse SWC
ancestors, which `fig_6E_G_BLP_RFs_PFs.m` never reads, so that skeleton code has been removed.

1. **Build the reference columns** — run `Data_Processing\s16_Mi1_Tm3_T4a_matching.m`. It
   loads the Codex tables (`consolidated_cell_types.csv`, `connections_no_threshold.csv`,
   `fafb_v783_princeton_synapse_table.csv`) from `Codex_Data\`, the optic-lobe neuropil
   meshes, and the Mi1 reference CSV. It builds `Mi1_Columns` (anchor columns
   with denoised medulla output synapses), then `Tm3_Columns` and `T4a_Columns` (each
   inheriting the synapse-weighted (theta, phi) of its upstream Mi1 columns, with denoised
   lobula / lobula-plate output synapses). All three are saved to
   `Processed_Data\Mi1_Tm3_T4a_columns.mat` (`-v7.3`).

2. **Fit the reference-layer surfaces** — run `Data_Processing\s17_making_layers.m`. For
   each reference cell type (Me = Dm6 dendrite, Lo = LT1a dendrite, LoP = T5a axon) it
   collects the reference synapses inside the target neuropil (right and left), denoises
   them, runs PCA, and fits a parabolic surface (poly22 for Me/Lo, poly33 for LoP). The
   per-hemisphere mean, PCA coefficients, and surface fit are stored in `PCA_Basis` and
   saved to `Processed_Data\neuropil_PCA_basis.mat`. (It also shows the per-config
   reference-layer fit, right and left.)

3. **Render the fields** — run `Figures\fig_6E_G_BLP_RFs_PFs.m`. It loads
   `Processed_Data\Mi1_Tm3_T4a_columns.mat`, `Processed_Data\neuropil_PCA_basis.mat`, and
   `Processed_Data\BLP_neurons_thr0.mat`, plus the Codex tables and optic-lobe meshes. For
   every BLP_R cell type it maps each input/output synapse to its nearest reference column's
   (theta, phi) (projecting onto the PCA + parabolic surface via `closest_point_poly22`),
   accumulates synapse-weighted 2D histograms over (phi, theta), and Gaussian-smooths them.
   It produces, **per cell type**, a per-neuron IN (orange = RF) + OUT (blue = PF) heatmap
   overlay and a type-level overlay of per-neuron contours (IN solid / OUT dashed). These
   per-cell-type figures are the panels used in
   **6E**, **6G**, and **S4**. The per-type result tables (`all_BLP_by_type`, including the
   RF/PF heatmaps) are saved to `Processed_Data\BLP_RFs_PFs.mat` (`-v7.3`). The
   `map_pos_to_weighted_angles`, `closest_point_poly22`, `distinguishable_colors`, and the
   heatmap/smoothing helpers are local functions at the bottom of the script.

## Figure 6F / 6H / 6J — BLP_R RF/PF size and center offset

Prerequisite: `Figures\fig_6E_G_BLP_RFs_PFs.m` must have been run first (it produces
`Processed_Data\BLP_RFs_PFs.mat`, `all_BLP_by_type` with the RF/PF heatmaps).

1. **Compute the metrics** — run `Data_Processing\s18_BLP_RFs_PFs_area_center.m`. From each
   neuron's normalized RF (input) and PF (output) heatmap it measures, at a contour
   threshold, the field size (deg²) and the weighted center (phi circular, theta linear),
   and the RF→PF center offset (total distance plus signed phi/theta components), including
   the phi-mirrored RF center. The metric columns are appended to `all_BLP_by_type` and
   saved (with a small `all_BLP_by_type_meta`) to
   `Processed_Data\BLP_RFs_PFs_area_center.mat` (`-v7.3`).

2. **Render the figures** — run `Figures\fig_6F_H_J_BLP_area_center_plot.m`. It loads
   `Processed_Data\BLP_RFs_PFs_area_center.mat` and produces (in creation order):

   - **Figure 1** — per-type grouped bar of RF vs PF area (not a paper panel).
   - **Figure 2** (paper panel **6F**) — per-type RF→PF center offset, phi (dx)
     vs theta (dy).
   - **Figure 3** (paper panel **6H**) — the same for the phi-mirrored RF center.
   - **Figure 4** (paper panel **6J**) — per-type mean RF area vs PF area scatter with a 1:1
     line; Pearson/Spearman correlations (individual neurons and per-type means) are printed
     to the console.

   A per-type mean RF/PF area table is also printed to the console.

# Supplementary figures

Some supplementary panels are produced by the main-figure scripts above; the rest have
their own scripts. They are listed below in panel order.

**Supplementary panels already produced by main-figure scripts** (which main figure makes them):

- **S1D** — `fig_2G_H_FFP_axo_axonic_RF_plot.m` (Figure 2G–H), connected vs unconnected RF distance bars
- **S2A** — `fig_4C_FBP_input_superclass.m` (Figure 4C)
- **S2B** — `fig_4D_FBP_dysynaptic_optclobe.m` (Figure 4D)
- **S2C / S2D** — `fig_4E_FBP_layer.m` (Figure 4E)
- **S2E** — `fig_4G_FBP_layer_overlap.m` (Figure 4G)
- **S4** — `fig_6E_G_BLP_RFs_PFs.m` (Figure 6E–G)

**Supplementary figures with their own scripts:** S1A–B (`fig_S1A_B_postPI_prePI_left_neurons.m`),
S1C (`fig_S1C_right_left_symmetry.m`), S1E–J (`fig_S1E_J_superclusters.m`), S3A–B
(`fig_S3A_B_BDP_reciprocity_OL_CB.m`), S3C–E (`fig_S3C_D_E_LC9_LT43_LT52_piechart.m`),
S5A–B (`fig_S5A_B_MeMe_connection_graph.m`), S5D–H (`fig_S5D_E_F_G_H_MeMe_simulation.m`).
Details below, in panel order.

## Figure S1A–B — postPI / prePI of left-hemisphere neurons

Prerequisite: `Data_Processing\s02_compute_postPI_prePI_all_neurons.m` must have been run
first (it produces `Processed_Data\FAFB_NPI_thr0.mat`).

Run `Figures\fig_S1A_B_postPI_prePI_left_neurons.m`. It loads `Processed_Data\FAFB_NPI_thr0.mat`
and computes, per neuron, the left-hemisphere PostPI / PrePI (left optic lobe vs central
brain). It selects left optic-lobe neurons in three steps — (1) keep neurons with ≥ 5
synapses bridging the left optic lobe and central brain, (2) drop neurons that are mostly in
the right optic lobe (> 70%), (3) drop cell types kept in fewer than 20% of their neurons —
then aggregates per cell type. It produces:

- **Figure 1** (paper panel **S1A**) — per-type mean left PostPI vs PrePI scatter, colored by
  superclass.
- **Figure 2** (paper panel **S1B**) — left PostPI + PrePI histogram for the visual_projection /
  visual_centrifugal / optic / central superclasses.
- **Figure 3** (paper panel **S1B**) — the same for left PostPI − PrePI.

## Figure S1C — right vs left symmetry

Prerequisite: `Data_Processing\s02_compute_postPI_prePI_all_neurons.m` must have been run
first (it produces `Processed_Data\FAFB_NPI_thr0.mat`).

Run `Figures\fig_S1C_right_left_symmetry.m`. It loads `Processed_Data\FAFB_NPI_thr0.mat`,
computes the right and left PostPI / PrePI, and selects the right and left optic-lobe neuron
sets using the same filters as `fig_1D_E_postPI_prePI_right_neurons.m` and
`fig_S1A_B_postPI_prePI_left_neurons.m` respectively (≥ 5 bridging synapses; drop neurons
dominated > 70% by the opposite optic lobe; drop cell types kept in < 20% of their neurons).
For each cell type present on both sides it measures the distance between its mean right
(PostPI, PrePI) and mean left (PostPI, PrePI) positions. It produces:

- **Figure 1** (paper panel **S1C**) — bar chart of the mean right-vs-left distance per
  superclass (VPN / VCN / Optic / Central), with error bars and jittered per-type points.

## Figure S1E–J — FFP supercluster output composition

Prerequisite: `Figures\fig_3C_D_E_F_G_FFP_clustering_plot.m` must have been run first (it
produces `Processed_Data\FFP_supercluster_targets.mat`, the per-supercluster post-synaptic
targets struct `SuperCluster_targets`).

Run `Figures\fig_S1E_J_superclusters.m`. It loads
`Processed_Data\FFP_supercluster_targets.mat` and, for each of the 7 FFP superclusters,
draws pie charts (small slices pooled into "Others"; the `collapse_small` and `draw_pie`
helpers are local functions):

- post-synaptic neuron-type distribution by neuron count (threshold 5%),
- post-synaptic neuropil distribution by synapse count (threshold 5%),
- post-synaptic neuron-type distribution by synapse count (threshold 1%).

## Figure S2 — FBP input / layer analyses (produced by the Figure 4 scripts)

These panels are generated by the main Figure 4 scripts; see the linked sections above.

- **S2A** — per-type FBP input superclass composition (raw and central-reassigned), from
  `Figures\fig_4C_FBP_input_superclass.m` (see Figure 4C).
- **S2B** — per-type disynaptic optic-lobe input split contra vs ipsi (Me / Lo / LoP), from
  `Figures\fig_4D_FBP_dysynaptic_optclobe.m` (see Figure 4D).
- **S2C** — reference-layer PCA + parabola fits per optic lobe, from
  `Data_Processing\s11_FBP_output_layer.m` and `s12_FBP_upstream_layer.m` (see Figure 4E).
- **S2D** — per-optic-lobe FBP output / disynaptic-input depth RGB images, from
  `Figures\fig_4E_FBP_layer.m` (see Figure 4E).
- **S2E** — per-neuron input/output depth overlap bar charts, from
  `Figures\fig_4G_FBP_layer_overlap.m` (see Figure 4G).

## Figure S3 — BDP reciprocity within the optic lobe vs central brain

Prerequisite: `Figures\fig_1D_E_postPI_prePI_right_neurons.m` must have been run first (it
produces `Processed_Data\right_neurons_thr0.mat`, `RightBDP_real_type` / `RightBDP_real_NPIs`).

Run `Figures\fig_S3A_B_BDP_reciprocity_OL_CB.m`. It loads
`Processed_Data\right_neurons_thr0.mat` and the Codex tables (`connections_no_threshold.csv`,
`consolidated_cell_types.csv`) from `Codex_Data\`. For every right-side real-bidirectional
(BDP) neuron it computes the weighted Jaccard index between its input-partner and
output-partner neurons **within the same compartment** — once restricted to the optic
lobe (`WJI_OL_total`) and once to the central brain (`WJI_CB_total`), with partners weighted
by synapse count — and averages the per-neuron values per cell type (the same neuron-level
convention as `s14_BDP_reciprocal.m`). It produces:

- **Figure 1** (paper panel **S3A**) — box plot of `WJI_OL_total` vs `WJI_CB_total` with the
  mean overlaid and a Wilcoxon signed-rank p-value.
- **Figure 2** (paper panel **S3B**) — scatter of `WJI_OL_total` vs `WJI_CB_total` across BDP
  types, with the top types labeled.

The paired (signed-rank) and unpaired (rank-sum) statistics are also printed to the console.
The helpers `reciprocity_within_npl`, `weighted_jaccard`, and `labelTopN_byScore` are local
functions at the bottom of the script.

## Figure S3C–E — LC9 / LT43 / LT52 input/output composition pie charts

Prerequisite: the Codex tables `connections_no_threshold.csv`, `consolidated_cell_types.csv`,
and `classification.csv` must be in `Codex_Data\`.

Run `Figures\fig_S3C_D_E_LC9_LT43_LT52_piechart.m`. It first builds each example neuron's
input/output partner tables directly (via the local function `seeConnection_by_region`, so no
precomputed `.mat` is needed). Then, for each of the three example bidirectional neuron types
— **LC9 (panel S3C)**, **LT43 (panel S3D)**, **LT52 (panel S3E)** — it draws 8 pie charts
(partner-type pies are thresholded by a per-neuron synapse count and colored by each partner
type's most common superclass; the `draw_overview_pie`, `draw_type_pie`, and
`seeConnection_by_region` helpers are local functions). The 8 figures per neuron are:

1. **Input overview** — total input synapses split optic lobe vs central brain.
2. **Output overview** — total output synapses split optic lobe vs central brain.
3. **Inputs total** — input partner-type composition, all neuropils.
4. **Outputs total** — output partner-type composition, all neuropils.
5. **Inputs from central** — input partner types from the central brain only.
6. **Inputs from optic** — input partner types from the optic lobes only.
7. **Outputs to central** — output partner types in the central brain only.
8. **Outputs to optic** — output partner types in the optic lobes only.

## Figure S4 — BLP_R RF/PF fields per cell type (produced by the Figure 6E–G script)

The per-cell-type RF/PF heatmap and contour figures from
`Figures\fig_6E_G_BLP_RFs_PFs.m` (see Figure 6E–G) are the panels used in **S4** (as well as
6E / 6G).

## Figure S5 — MeMe network connectivity and rate-model simulation

The Figure S5 panels are a small rate-model simulation of the MeMe / Sm / Dm2 /
MeTu cell-type network and the connectivity it is built from. The pipeline is
three data-generation steps plus two figure scripts.

**Note on data:** the precomputed `.mat` inputs (connection graph, per-type RFs,
stimuli; some > 100 MB) are **not** shipped — regenerate them with the s-scripts
below. The RF step reuses the optic-lobe reference data from the main Figure 6
pipeline.

1. **Generate the rolling stimuli** — run `Data_Processing\s19_MeMe_roll_stimuli.m`.
   It writes the two "sky roll" stimuli (upper hemisphere bright, constant angular
   speed, 1 s pre + 2 s roll + 1 s post @ 50 fps) to `Processed_Data\`:
   `Roll_60deg_cw.mat` / `.mp4` and `Roll_60deg_ccw.mat` / `.mp4` (`pat`:
   arenaheight × arenawidth × nFrames). Edit `total_roll_deg` to change the angle.

2. **Build the connection graph** — run `Data_Processing\s20_MeMe_connection.m`.
   It reads `connections_no_threshold.csv` and `consolidated_cell_types.csv` from
   `Codex_Data\`, removes minority-NT (connectome-error) cells and sparse type-pair
   edges, and saves the neuron-level connectivity used by the simulation —
   `pair_table` (signed synapses per (pre,post) pair), `neuron_table`, and the
   type-pair `edge_table` — to `Processed_Data\MeMe_AllSimNeurons_Graph.mat`. No
   figure is drawn (the connection-graph figure is the L/R-split one below).

3. **Build the per-type input RFs** — run `Data_Processing\s21_MeMe_generate_RFs.m`.
   For each simulation type (Dm2, MeMe_e01, MeMe_e02, Sm07, MeTu1, MeTu3c) it maps
   every neuron's medulla input synapses to visual-space (theta, phi) via the
   nearest Mi1 reference column (projected onto the PCA + parabolic medulla layer),
   accumulates a synapse-weighted, Gaussian-smoothed RF, and writes
   `Processed_Data\<type>_RFs.mat` (`Hin_norm` per neuron). It reuses
   `Processed_Data\Mi1_Tm3_T4a_columns.mat` (`Mi1_Columns`, from
   `s16_Mi1_Tm3_T4a_matching.m`) and `Processed_Data\neuropil_PCA_basis.mat`
   (`PCA_Basis`, from `s17_making_layers.m`), the medulla meshes
   `Processed_Data\optic_lobe_neuropil_mesh\Me_{R,L}_{vertices,faces}.csv` (from
   `s01`), and `fafb_v783_princeton_synapse_table.csv` /
   `consolidated_cell_types.csv` / `classification.csv` from `Codex_Data\`. The
   `intriangulation` helper (`Helper_Function\`) must be on the MATLAB path.

Then the two figure scripts:

**Figure S5A–B — L/R-split connection graph** — run
`Figures\fig_S5A_B_MeMe_connection_graph.m`. From the Codex tables
(`connections_no_threshold.csv`, `consolidated_cell_types.csv`,
`classification.csv`) it splits each cell type into left/right nodes and builds the
type-pair connectivity (edge weight = mean synapse count over connected (pre,post)
pairs; sign from the dominant neurotransmitter, GABA/GLUT = inhibitory). It produces:

- **Figure 1** (paper panel **S5A**) — the L/R-split connection digraph (orange =
  excitatory, blue = inhibitory; edge width ∝ |weight|; nodes draggable, edge
  midpoints clickable for details).
- **Figure 2** (paper panel **S5B**) — the same edge weights as a signed
  source(pre) × target(post) connection matrix.

**Figure S5D–H — rate-model simulation** — run
`Figures\fig_S5D_E_F_G_H_MeMe_simulation.m`. It loads the three Processed_Data
inputs above and integrates `tau dr/dt = -r + f(baseline + drive(t) + W r)`
(`f = tanh`, `r ≥ 0`) for the intact network ("Normal") and one with the
MeMe_e01 ↔ MeMe_e02 interconnections cut ("Ablated"). Run it **once per stimulus**
via the `stim_name` setting:

- **Figure 1** (paper panel **S5D**) — RF-quadrant subgroup color key (stim-independent).
- **Figures 2 & 3** (paper panels **S5E** for `stim_name = 'Roll_60deg_cw'`,
  **S5G** for `'Roll_60deg_ccw'`) — per-type mean response by RF subgroup, Normal
  (Figure 2) and Ablated (Figure 3).
- **Figure 4** (paper panel **S5F** for cw, **S5H** for ccw) — per-type Left vs Right
  mean final activity (line slope = L-R asymmetry), Normal (blue) vs Ablated (red).
  Left/right and upper/lower are split at each type's RF center (median φ / median θ);
  the annotation `L-R: dn → da (p%)` gives the normal/ablated left-right difference
  and the percent change in `|L-R|`. Dm2 excluded.
- **Figure 5** (Supplementary **Video S1** for cw, **Video S2** for ccw) — per-type
  spatial response movie over time: Normal vs Ablated vs their difference
  (Ablated − Normal), one row per type. Saved as `MeMe_sim_video__<stim>.mp4` in the
  project root; set `make_video = false` to skip.

So the cw run gives S5D / S5E / S5F / Video S1 and the ccw run gives S5G / S5H /
Video S2. The `build_activation`, `rf_quadrant_color` and `vid_col_vals` helpers are
local functions at the bottom of the script.

## Figure S6 — MCNS (male-CNS) cross-dataset replication

The Figure S6 panels are the male-CNS (MCNS) replication of the main FAFB analyses. The full
MCNS pipeline and the figures it produces are documented in the **MCNS (male-CNS) dataset**
section below.

## Folder structure

```
Connectome_Visual_Projection_Analysis/
├── Codex_Data/         # raw FlyWire Codex downloads (you provide these)
├── Reference_Data/     # external reference tables (you provide these)
│   └── Mi1_columns_DustinGarner.csv   # Mi1 columns (Garner et al. 2024, Nature)
├── Data_Processing/    # scripts/notebooks that generate data
│   ├── s01_fetch_neuropil_neuron.ipynb     # downloads optic-lobe neuropil + Fig 1A/1B/1H/5A meshes
│   ├── s02_compute_postPI_prePI_all_neurons.m  # per-neuron synapse / NPI table
│   ├── s03_FFP_FBP_mesh.ipynb               # downloads FFP / FBP neuron meshes
│   ├── s04_compute_graph_centrality.m       # graph centrality (PageRank, betweenness)
│   ├── s05_compute_FFP_optic_central_fan_in_out.m  # per-type fan-in/out + centrality
│   ├── s06_FFP_axo_axonic_compute.m          # Fig 2F: per-type axo-axonic same-type input + OL proximity-probability reference
│   ├── s07_FFP_axo_axonic_RF.m               # Fig 2G-H/S1D: outlier RF distance (OL input-synapse centroid dist, µm) binned by same-type CB strength W
│   ├── s08_FFP_ADJ_matrix.m                  # FFP -> central-brain output adjacency matrix
│   ├── s09_FFP_leiden_clustering.ipynb       # bipartite Leiden clustering of the matrix
│   ├── s10_FBP_out_opticlobes.m              # per-FBP-type output neuropil fractions
│   ├── s11_FBP_output_layer.m               # S2C: FBP output-synapse depth per optic lobe
│   ├── s12_FBP_upstream_layer.m             # FBP upstream-synapse depth per optic lobe
│   ├── s13_segregation_index.ipynb            # per-neuron axon/dendrite segregation index
│   ├── s14_BDP_reciprocal.m                    # FFP/FBP/BDP/Others input-output partner reciprocity (WJI)
│   ├── s15_BLP_mesh.ipynb                      # downloads BLP_R neuron meshes
│   ├── s16_Mi1_Tm3_T4a_matching.m              # Mi1/Tm3/T4a reference columns (theta/phi)
│   ├── s17_making_layers.m                     # reference-layer PCA bases (neuropil_PCA_basis.mat)
│   ├── s18_BLP_RFs_PFs_area_center.m           # RF/PF area & center metrics
│   ├── s19_MeMe_roll_stimuli.m                 # Fig S5: rolling visual stimuli (Roll_60deg_cw/ccw)
│   ├── s20_MeMe_connection.m                   # Fig S5: neuron-level connection graph for the simulation
│   └── s21_MeMe_generate_RFs.m                 # Fig S5: per-type input RFs (<type>_RFs.mat)
├── Processed_Data/     # generated data (meshes, .mat) — created by the scripts
└── Figures/            # figure-generating scripts
    ├── fig_1A_neurons.m
    ├── fig_1B_LC9.m
    ├── fig_1D_E_postPI_prePI_right_neurons.m
    ├── fig_1H_neurons_specific_angle.m
    ├── fig_2A_FFP_neurons.m
    ├── fig_2B_FFP_neuropil_connection.m
    ├── fig_2C_D_E_plot_FFP_optic_central.m
    ├── fig_2F_G_FFP_axo_axonic_plot.m
    ├── fig_2F_MeTu3c_axoaxonic_matrix.m       # Fig 2F: MeTu3c same-type central axo-axonic matrix (clustered)
    ├── fig_2G_H_FFP_axo_axonic_RF_plot.m
    ├── fig_3A_FFP_per_CB_neurons.m
    ├── fig_3C_D_E_F_G_FFP_clustering_plot.m
    ├── fig_4A_FBP_neurons.m
    ├── fig_4B_FBP_neuropil_connection.m
    ├── fig_4C_FBP_input_superclass.m
    ├── fig_4D_FBP_dysynaptic_optclobe.m
    ├── fig_4E_FBP_layer.m
    ├── fig_4F_FBP_layer_overall_histogram.m
    ├── fig_4G_FBP_layer_overlap.m
    ├── fig_4H_FBP_loop.m
    ├── fig_5A_BDP_neurons.m
    ├── fig_5B_C_D_BDP_postPI_prePI.m
    ├── fig_5E_segregation_index_plot.m
    ├── fig_5F_G_segregation_reciprocal.m
    ├── fig_5I_J_LC9_interconnection.m
    ├── fig_6B_postPI_prePI_BLP_neurons.m
    ├── fig_6C_BLP_neurons.m
    ├── fig_6D_BLP_neuropil_connection.m
    ├── fig_6E_G_BLP_RFs_PFs.m
    ├── fig_6F_H_J_BLP_area_center_plot.m
    ├── fig_S1A_B_postPI_prePI_left_neurons.m
    ├── fig_S1C_right_left_symmetry.m
    ├── fig_S1E_J_superclusters.m
    ├── fig_S3A_B_BDP_reciprocity_OL_CB.m
    ├── fig_S3C_D_E_LC9_LT43_LT52_piechart.m
    ├── fig_S5A_B_MeMe_connection_graph.m       # S5A/B: L/R-split connection graph + matrix
    └── fig_S5D_E_F_G_H_MeMe_simulation.m       # S5D-H + Video S1/S2: MeMe network rate-model simulation
```

# MCNS (male-CNS) dataset

The analyses above were built on the FAFB (FlyWire Codex) dataset. The `MCNS` subfolder
mirrors this project's structure (`MCNS_Data` / `Data_Processing` / `Processed_Data` /
`Figures`) to repeat the same kind of analysis on the male-CNS (`male-cns:v0.9`) connectome
from neuprint (https://neuprint.janelia.org).

## 0. Export the male-CNS data

Unlike FAFB (where you download the Codex CSVs by hand), the MCNS raw tables are exported
from neuprint with Python. Requirements:

- **neuprint-python** — `pip install neuprint-python` (https://github.com/connectome-neuprint/neuprint-python).
  `s04_export_synapse_coordinates.py` also uses neuprint's internal query helpers
  (`neuprint.queries.synapses._fetch_synapse_connections`, `neuprint.queries.neuroncriteria.NeuronCriteria`),
  so a recent version is required.
- a **neuprint auth token** — get it from https://neuprint.janelia.org (Account menu). Either set
  the `NEUPRINT_APPLICATION_CREDENTIALS` environment variable, or save the token to
  `MCNS/Data_Processing/.neuprint_token`.
- `tqdm` and `ujson` (`pip install tqdm ujson`).

Run the four export scripts (each writes one CSV into `MCNS/MCNS_Data/`, dataset
`male-cns:v0.9` by default; pass `--server` / `--dataset` / `--output` to override):

1. `Data_Processing/s01_export_connections.py` → `male-cns-v0.9-all-connections.csv`
   (`pre_root_id, post_root_id, neuropil, syn_count, nt_type`; primary ROIs + a "NotPrimary" row).
2. `Data_Processing/s02_export_classification.py` → `male-cns-v0.9-classification.csv`
   (`root_id, super_class, side`; super_class normalized to ascending / descending / central / optic).
3. `Data_Processing/s03_export_primary_types.py` → `male-cns-v0.9-primary-types.csv`
   (`root_id, primary_type, flywireType`).
4. `Data_Processing/s04_export_synapse_coordinates.py` → `male-cns-v0.9-synapse-coordinates.csv`
   (`pre_x, pre_y, pre_z, pre_root_id, post_root_id, neuropil`; one row per synapse). This is large
   and resumable — progress is tracked in `<output>.progress.json` and partial rows in `<output>.part`.

These four CSVs are the MCNS analogue of the FAFB Codex tables and feed the downstream
MCNS analysis. (A fifth neuprint export, `s10_export_roi_mesh_csv.py`, provides the optic-lobe
ROI meshes used by the Figure 4 layer analysis; it is described in the Figure 4 section below.)

## MCNS — per-neuron NPI table

Run `MCNS/Data_Processing/s05_compute_postPI_prePI_all_neurons.m` (the MCNS analogue of the
FAFB `s02_compute_postPI_prePI_all_neurons.m`). It reads the three exported tables
(`male-cns-v0.9-primary-types.csv`, `male-cns-v0.9-all-connections.csv`,
`male-cns-v0.9-classification.csv`) from `MCNS/MCNS_Data/`, groups neuropils into right optic
lobe / left optic lobe / central brain (central = everything but the optic lobes and
`NotPrimary`), computes each neuron's In/Out synapse counts per group, and saves the table
`MCNSNPIs` to `MCNS/Processed_Data/MCNS_NPI_thr0.mat`.

## MCNS — S6A: postPI / prePI of right-hemisphere neurons

Run `MCNS/Figures/fig_1_postPI_prePI_right_neurons.m` (the MCNS analogue of the FAFB
`fig_1D_E_postPI_prePI_right_neurons.m`). It loads `MCNS/Processed_Data/MCNS_NPI_thr0.mat`,
computes the right/left postPI / prePI, selects right optic-lobe neurons (≥ 5 bridging
synapses; drop neurons > 70% in the left optic lobe; drop types kept in < 20% of their
neurons), aggregates per cell type, and produces three figures:

- Figure 1: per-type mean postPI vs prePI scatter.
- Figure 2: PostPI + PrePI histogram.
- Figure 3: PostPI − PrePI histogram.

It also classifies the right-hemisphere types into feedforward (FFP), feedback (FBP), and
bidirectional (BDP, further split into optic / central / real) groups, saves these tables to
`MCNS/Processed_Data/right_neurons_thr0.mat`, and writes the FFP / FBP / BDP / Others root_id
lists (no header, root_id only) to `MCNS/Processed_Data/right_*_root_ids.csv`.

## MCNS — S6B: fan-in/out of Optic / FFP / Central neurons

1. **Build the fan-in/out table** — run
   `MCNS/Data_Processing/s06_compute_FFP_optic_central_fan_in_out.m` (the MCNS analogue of the
   FAFB `s05_compute_FFP_optic_central_fan_in_out.m`, without the graph-centrality columns).
   For the right FFP, optic, and central neurons it summarizes each neuron's input/output
   partners grouped by partner type (partner-neuron count, partner-type count, synapse count)
   and aggregates per cell type, saving `type_RightFFP_InOut`, `type_Optic_InOut`,
   `type_Central_InOut` to `MCNS/Processed_Data/FFP_optic_central_fan_in_out.mat`.

2. **Render the figures** — run `MCNS/Figures/fig_2_FFP_optic_central.m`. It loads
   `FFP_optic_central_fan_in_out.mat` and draws Optic / Feedforward / Central box plots
   (with Wilcoxon rank-sum tests printed to the console) for in-neuron number / type number /
   ratio and out-neuron number / type number / ratio (panel **S6B**). The `draw_metric_boxplot`,
   `pairwise_ranksum`, and `ternary` helpers are local functions.

## MCNS — S6C: Leiden clustering of FFP output projections

Prerequisite: `MCNS/Figures/fig_1_postPI_prePI_right_neurons.m` must have been run first,
since it produces `MCNS/Processed_Data/right_neurons_thr0.mat` (`RightFFP_NPIs`). This is the
MCNS analogue of the FAFB Figure 3C–G pipeline, built in three steps: build the adjacency
matrix (MATLAB), run the Leiden clustering (Python), then render the figures (MATLAB).

1. **Build the FFP output adjacency matrix** — run
   `MCNS/Data_Processing/s07_FFP_ADJ_matrix.m` (the MCNS analogue of the FAFB
   `s08_FFP_ADJ_matrix.m`). It loads `MCNS/Processed_Data/right_neurons_thr0.mat`
   (`RightFFP_NPIs`) and `male-cns-v0.9-all-connections.csv` / `male-cns-v0.9-classification.csv`
   from `MCNS/MCNS_Data/`, keeps the FFP → central-brain output connections (optic-lobe
   neuropils dropped, VPN targets dropped), and builds the FFP (rows) × post-neuron (cols)
   synapse matrix `WantMatrix_Out`. It writes to `MCNS/Processed_Data/`:
   - `right_FFP_table.csv` — FFP neuron table (matrix rows)
   - `right_FFP_output_matrix_CB_no_VPN_thr0.csv` / `.mat` — the adjacency matrix
   - `post_right_FFP_CB_no_VPN_thr0.csv` — post-neuron root_ids (matrix columns)
   - `post_neurons_FFP_opticlobe_central_no_VPN_thr0.mat` — same post-neuron ids (`Post_Want`)

2. **Run the Leiden clustering** — run
   `MCNS/Data_Processing/s08_FFP_leiden_clustering.ipynb` (the MCNS analogue of the FAFB
   `s09_FFP_leiden_clustering.ipynb`; requires the `scikit-network` package). It reads the
   three CSVs from step 1 and runs bipartite Leiden clustering (resolution 3.6) on the
   adjacency matrix. Row clusters that contain a single neuron (singletons) and the post-neuron
   columns connected only to those singletons are then removed, and the biadjacency is
   re-aggregated by cluster. It writes to `MCNS/Processed_Data/`:
   - `leiden_right_FFP_output_CB_no_VPN_thr0_resol3.6_100000_row_singletons_removed.csv` — FFP-row cluster assignments
   - `leiden_post_right_FFP_CB_no_VPN_thr0_100000_row_singletons_removed.csv` — post-column cluster assignments
   - `leiden_right_FFP_intercluster_connectivity_CB_no_VPN_thr0_resol3.6_100000_row_singletons_removed.csv` — inter-cluster connectivity
   - `dendrogram_raw_CB_no_VPN_thr0_resol3.6_100000_row_singletons_removed.csv` — raw-similarity dendrogram (linkage matrix)
   - `leiden_right_FFP_kept_rows_*` / `leiden_right_FFP_removed_rows_*` and
     `leiden_post_right_FFP_kept_cols_*` / `leiden_post_right_FFP_removed_cols_*` — kept/removed row/column bookkeeping

3. **Render the figures** — run `MCNS/Figures/fig_3_FFP_clustering_plot.m` (the MCNS analogue
   of the FAFB `fig_3C_D_E_F_G_FFP_clustering_plot.m`, trimmed to the two clustering-overview
   panels). It loads the s07 / s08 outputs above plus `male-cns-v0.9-all-connections.csv` and
   `male-cns-v0.9-primary-types.csv` from `MCNS/MCNS_Data/`, labels each FFP-row and
   post-column cluster by its dominant cell types (down to 5% of the cluster size), and
   produces:
   - **Figure 1** — FFP-cluster dendrogram, reordered by leaf order.
   - **Figure 2** — inter-cluster connectivity heatmap, reordered by leaf order.

   The 32 Leiden clusters are also grouped into 11 hand-defined **superclusters**
   (`SuperClusterOrder`); for each, the script tabulates the post-synaptic neuron-type
   counts, neuropil synapse counts, and per-type synapse counts into the struct
   `SuperCluster_targets`, saved to `MCNS/Processed_Data/FFP_supercluster_targets.mat`
   (computed but not plotted here). The `get_leaf_order_from_linkage` helper is a local
   function at the bottom of the script.

## MCNS — FBP output-neuropil classification

Prerequisite: `MCNS/Figures/fig_1_postPI_prePI_right_neurons.m` must have been run first,
since it produces `MCNS/Processed_Data/right_neurons_thr0.mat` (`RightFBP_by_type`).

Run `MCNS/Data_Processing/s09_FBP_out_opticlobes.m` (the MCNS analogue of the FAFB
`s10_FBP_out_opticlobes.m`). It loads `MCNS/Processed_Data/right_neurons_thr0.mat`
(`RightFBP_by_type`) and the MCNS tables (`male-cns-v0.9-primary-types.csv`,
`male-cns-v0.9-all-connections.csv`) from `MCNS/MCNS_Data/`, and, for every FBP cell type,
computes the fraction (%) of its output synapses landing in each optic-lobe region,
collapsing left/right into AME / ME / LO / LOP (stored in the variable
`RightFBP_OutNeuropils`). These values were saved by hand (separately, outside the script)
to `MCNS/Processed_Data/FBP_output_neuropils.xlsx`, which was then used to classify each FBP
neuron type by the optic-lobe region it targets. The resulting target grouping (Me / Lo /
Lop / Multi) is reused by the downstream MCNS S6D code.

(The local function `seeConnection_by_region`, which summarizes each type's input/output
partners and neuropils, is defined at the bottom of the script.)

## MCNS — S6D: disynaptic optic-lobe input to FBP neurons

Prerequisite: `MCNS/Figures/fig_1_postPI_prePI_right_neurons.m` must have been run first
(it produces `MCNS/Processed_Data/right_neurons_thr0.mat`, `RightFBP_NPIs`). The FBP target
grouping (Me / Lo / Lop / Multi) comes from the FBP output-neuropil classification
(`MCNS/Processed_Data/FBP_output_neuropils.xlsx`, see above); the hardcoded index lists are
embedded in the script.

Run `MCNS/Figures/fig_4A_FBP_dysynaptic_optclobe.m` (the MCNS analogue of the FAFB
`fig_4D_FBP_dysynaptic_optclobe.m`). It loads `MCNS/Processed_Data/right_neurons_thr0.mat`
and the MCNS tables (`male-cns-v0.9-primary-types.csv`, `male-cns-v0.9-all-connections.csv`)
from `MCNS/MCNS_Data/`. For every FBP cell type it estimates the **disynaptic optic-lobe
input** along the path `optic lobe -> upstream (central, non-optic) neuron -> FBP neuron`:
it finds each FBP type's upstream non-optic partners, then measures how much optic-lobe input
those upstream neurons receive, weighted by their synapse contribution, broken down into six
regions (Me R, Me L, Lo R, Lo L, LoP R, LoP L). FBP types are grouped into Me / Lo / Lop /
Multi by their target optic lobe. It produces:

- **Figure 1** — `imagesc` of the mean ± std disynaptic optic-lobe input fraction per neuropil
  (rows) × FBP group (columns: CB-Me / CB-LO / CB-LOP / CB-Multi), with the mean ± std
  annotated in each cell.
- **Figure 2** — per Me-targeting FBP type, the disynaptic Me input split into contralateral
  (Me L) vs ipsilateral (Me R).
- **Figure 3** — the same for Lo-targeting types (Lo L vs Lo R).
- **Figure 4** — the same for Lop-targeting types (LoP L vs LoP R).

The upstream-partner helper `seeConnection_root_id_NoOptic` is defined as a local function at
the bottom of the script.

## MCNS — S6D: FBP innervation-depth layers

This is the MCNS analogue of the FAFB Figure 4E–G depth-layer analysis, built from the optic-lobe
ROI meshes plus three Data_Processing steps (s10 / s11 / s12) and three figure scripts. Common
prerequisites: `MCNS/Figures/fig_1_postPI_prePI_right_neurons.m` must have been run first (it
produces `MCNS/Processed_Data/right_neurons_thr0.mat`, `RightFBP_NPIs`), and the FBP target
grouping (Me / Lo / Lop / Multi) comes from the FBP output-neuropil classification
(`MCNS/Processed_Data/FBP_output_neuropils.xlsx`, see above); the hardcoded index lists are
embedded in each script.

1. **Export the ROI meshes** — run `MCNS/Data_Processing/s10_export_roi_mesh_csv.py` (a neuprint
   export, the MCNS analogue of the FAFB `s01_fetch_neuropil_neuron.ipynb` mesh fetch; uses the
   same neuprint client / token as the other exports). For each ROI in
   `{ME(R), ME(L), LO(R), LO(L), LOP(R), LOP(L)}` it writes `vertices.csv` (`x, y, z`, in 8 nm
   voxels) and `faces.csv` (`v1, v2, v3`, zero-based) into
   `MCNS/Processed_Data/roi_mesh_csv/<ROI>/`.

2. **Build the FBP output-synapse depth profiles** — run
   `MCNS/Data_Processing/s11_FBP_output_layer.m` (the MCNS analogue of the FAFB
   `s11_FBP_output_layer.m`). It loads `MCNS/Processed_Data/right_neurons_thr0.mat`,
   `male-cns-v0.9-synapse-coordinates.csv` and `male-cns-v0.9-primary-types.csv` from
   `MCNS/MCNS_Data/`, and the ROI meshes from `MCNS/Processed_Data/roi_mesh_csv/`. For each optic
   lobe it fits a **reference layer** (PCA + parabolic surface; Me = `Dm6` dendrite,
   Lo = `LT1a` dendrite, LoP = `T5a` axon) and measures the depth of the FBP group's **own output
   (axon) synapses** relative to that layer, for the right (ipsilateral) and left (contralateral)
   hemispheres. Unlike the FAFB version (which tests mesh membership with `intriangulation`),
   synapses are assigned to a neuropil by the synapse table's `neuropil` column. It saves the six
   matrices `FBP_Synapse_{ME,LO,LOP}_{R,L}` to
   `MCNS/Processed_Data/FBP_output_synapse_layer.mat`. The `load_roi_mesh` helper is a local
   function at the bottom of the script.

3. **Build the FBP upstream-synapse depth profiles** — run
   `MCNS/Data_Processing/s12_FBP_upstream_layer.m` (the MCNS analogue of the FAFB
   `s12_FBP_upstream_layer.m`). Same reference layers, but it profiles the in-neuropil synapses of
   each FBP type's **upstream (central, non-optic) partners** (`seeConnection_root_id_NoOptic`),
   weighted by each partner's synapse contribution. It saves the six matrices
   `FBP_Upstream_Synapse_{ME,LO,LOP}_{R,L}` to
   `MCNS/Processed_Data/FBP_upstream_synapse_layer.mat`.

Three figure scripts then load the two `.mat` files from steps 2–3:

- `MCNS/Figures/fig_4S_FBP_layer_figure.m` (the MCNS analogue of the FAFB `fig_4E_FBP_layer.m`)
  also loads `right_neurons_thr0.mat` and, for each optic lobe (ME / LO / LOP), renders one RGB
  image whose columns are per-FBP-type triplets — upstream input ipsilateral (blue) · FBP output
  (orange) · upstream input contralateral (green) — stacked over the innervation-depth axis
  (Figures 1–3).
- `MCNS/Figures/fig_4B_FBP_layer_overall_histogram.m` (the MCNS analogue of the FAFB
  `fig_4F_FBP_layer_overall_histogram.m`) collapses all FBP types into a single mean ± 1 SD profile
  per optic lobe (ME / LO / LOP), overlaying FBP output (orange), upstream input ipsilateral
  (blue), and upstream input contralateral (green) (`shaded_band` is a local function).
- `MCNS/Figures/fig_4C_FBP_layer_overlap.m` (the MCNS analogue of the FAFB
  `fig_4G_FBP_layer_overlap.m`) also loads `right_neurons_thr0.mat`. For each FBP neuron it
  Gaussian-smooths the output depth profile and the upstream-input depth profiles (each scaled by
  its share of the total upstream synapses), then computes the **overlap coefficient** as the
  per-depth `sum(min(...))`, separately for the contralateral (L) and ipsilateral (R) input. It
  produces per-neuron overlap bar charts for ME / LO / LOP (Figures 1–3) and a box plot of the overlap distribution
  across all neurons grouped by region and side (Figure 4); the overall mean ± std is printed to
  the console.

(In addition, `MCNS/Figures/fig_4A_FBP_dysynaptic_optclobe.m`, above, covers the disynaptic
optic-lobe input branch of MCNS S6D; it works directly from the connection table and does not
use these depth-layer `.mat` files.)

## MCNS — S6E: segregation index vs input/output reciprocity

This is the MCNS analogue of the FAFB Figure 5E–G segregation / reciprocity analysis, built from
two Data_Processing steps (s13 / s14) and one figure script. Common prerequisite:
`MCNS/Figures/fig_1_postPI_prePI_right_neurons.m` must have been run first — it produces
`MCNS/Processed_Data/right_neurons_thr0.mat` and writes the per-group root_id lists
(`right_FFP_root_ids.csv`, `right_FBP_root_ids.csv`, `right_BDP_root_ids.csv`,
`right_Others_root_ids.csv`) to `MCNS/Processed_Data/`.

1. **Compute the segregation index** — run `MCNS/Data_Processing/s13_segregation_index.ipynb` (the
   MCNS analogue of the FAFB `s13_segregation_index.ipynb`; requires `neuprint-python`, `navis`, and
   a neuprint token). For each group (FFP / FBP / BDP / Others) it reads the root_id list, fetches
   each neuron's skeleton and synapses from neuprint (`male-cns:v0.9`), splits it into axon/dendrite,
   and computes its segregation index. It writes `MCNS/Processed_Data/right_<group>_segregation_index.csv`
   (columns: `root_id`, `segregation_index`) for each group. (Unlike the FAFB version, which
   skeletonizes from FlyWire, the skeleton/synapses are fetched from neuprint.)

2. **Compute input/output reciprocity** — run `MCNS/Data_Processing/s14_BDP_reciprocal.m` (the MCNS
   analogue of the FAFB `s14_BDP_reciprocal.m`). For every right-hemisphere FFP, FBP, real-bidirectional
   (BDP), and other bidirectional (Others = optic + central BD) neuron it measures how much the
   input-partner set overlaps with the output-partner set, using the synapse-weighted Jaccard index
   (Σmin / Σmax over the union of partner neuron ids). It loads `male-cns-v0.9-all-connections.csv`
   from `MCNS/MCNS_Data/` and `right_neurons_thr0.mat`, averages the per-neuron values per cell type,
   keeps the three example BDP types (LC9 / LT43 / LT52) per neuron in `Want`, and saves `Type_FFP`,
   `Type_FBP`, `Type_BDP`, `Type_Others`, `Want` to `MCNS/Processed_Data/BDP_reciprocity.mat`. Because
   the MCNS connection table is large, each group's input/output partner weights are pre-aggregated
   once with `groupsummary` into lookup Maps before the per-neuron Jaccard. The `reciprocity_for_neurons`,
   `aggregate_id_map`, `weighted_jaccard`, and `summarize_reciprocity` helpers are local functions.

3. **Render the figures** — run `MCNS/Figures/fig_5A_segregation_reciprocal.m` (the MCNS analogue of
   the FAFB `fig_5F_G_segregation_reciprocal.m`). It loads `BDP_reciprocity.mat` (step 2) and attaches
   the per-type mean segregation index from the `right_<group>_segregation_index.csv` files (step 1,
   by `root_id`). It produces:
   - **Figure 1** — per-type segregation index (x) vs reciprocity (y) scatter for FFP / FBP / BDP /
     Others, with a linear fit + 95% prediction interval (fit statistics annotated, fit over all four
     groups) and the three example BDP types highlighted.
   - **Figure 2** — reciprocity box plot across LC9 / LT43 / LT52 / BDP / FFP / FBP.

   The BDP / FFP / FBP pairwise Wilcoxon rank-sum tests (Bonferroni-corrected) and probabilistic-dominance
   values are printed to the console. The `attach_mean_segregation` and `ternary` helpers are local
   functions at the bottom of the script.

## MCNS — S6E: LC9–LC9 interconnection vs proximity

Prerequisite: `MCNS/Figures/fig_1_postPI_prePI_right_neurons.m` must have been run first (it produces
`MCNS/Processed_Data/right_neurons_thr0.mat`, `RightBDP_real_NPIs`, used here for the right-hemisphere
LC9 neurons).

Run `MCNS/Figures/fig_5B_LC9_interconnection.m` (the MCNS analogue of the FAFB
`fig_5I_J_LC9_interconnection.m`). It loads `male-cns-v0.9-all-connections.csv` and
`male-cns-v0.9-synapse-coordinates.csv` from `MCNS/MCNS_Data/` and the LC9 root_ids from
`right_neurons_thr0.mat`. The optic-lobe and central-brain neuropil groups are defined inline (the
optic lobes hardcoded, central = every other neuropil except `NotPrimary`, the same convention as
`s05_compute_postPI_prePI_all_neurons.m`). It builds the LC9 -> LC9 central-brain interconnection
matrix (synapse counts), and for each LC9 neuron averages its input synapses (the neuron is post) inside the right
optic lobe into a single OL centroid (synapses assigned to the optic lobe by the synapse table's
`neuropil` column, rather than by mesh `intriangulation` as in the FAFB version). Each pair's distance
is then the Euclidean distance between the two neurons' OL synapse centroids, and connection weight is
related to proximity over all neuron pairs. It produces:

- **Figure 1** — box plot of inter-centroid distance for pairs that are strongly connected (W ≥ 5),
  connected (W > 0), or unconnected (W = 0), with the connected-vs-unconnected Wilcoxon rank-sum
  p-value in the title.
- **Figure 2** — connection probability vs distance (grouped bars for W > 0 and W ≥ 5 on the left
  axis) overlaid with the per-bin pair count on the right axis.

Synapse coordinates are scaled from 8 nm voxels to metres, so distances and the plot axes are in
metres.

## MCNS — S6F: BLP_R projective / receptive fields in visual space

This is the MCNS analogue of the FAFB Figure 6E–G / 6F–J RF/PF analysis. It maps each BLP_R
neuron's input synapses (receptive field, RF) and output synapses (projective field, PF) into
visual-space (theta, phi) coordinates, built from four Data_Processing steps (s15–s18) and two
figure scripts. Unlike the FAFB version (which assigns synapses to a neuropil via mesh
`intriangulation` and uses SWC skeletons), the MCNS scripts assign synapses to a region by the
synapse table's `neuropil` column and do not use skeletons.

Common prerequisites (required inputs not produced by the s-scripts below):
- `MCNS/Processed_Data/BLP_neurons_thr0.mat` (`BLP_R_NPIs`, `BLP_R_by_type`) — the bilateral (BLP)
  neuron classification (MCNS analogue of the FAFB `fig_6B_postPI_prePI_BLP_neurons.m` output).
- `Reference_Data/Mi1_columns_DustinGarner.csv` — the FlyWire Mi1 column (theta/phi) table from
  Garner et al. 2024, *Nature* (shared with the FAFB `Reference_Data`).
- The optic-lobe ROI meshes in `MCNS/Processed_Data/roi_mesh_csv/` (from `s10_export_roi_mesh_csv.py`).

1. **Match neuprint Mi1 to FlyWire Mi1** — run
   `MCNS/Data_Processing/s15_FAFB_MCNS_Mi1_matching.ipynb` (requires `neuprint-python`, `fafbseg`,
   `navis`, and a neuprint token). For each neuprint Mi1 it fetches the skeleton, finds the FlyWire
   Mi1 with the smallest directed-mean-nearest skeleton distance, and writes
   `neuprint-mi1-to-flywire-v783-directed-mean-nearest-skeleton-matches.csv` — the match table read
   downstream by `s16` — to `MCNS/Processed_Data/` (its default `output_dir`). It writes a few
   bookkeeping files alongside it (per-side skeleton summaries, a skeleton-node `.npz`, and
   fetch-issue CSVs); only the matches CSV is needed downstream.

2. **Build the Mi1 / Tm3 / T4a reference columns** — run
   `MCNS/Data_Processing/s16_Mi1_Tm3_T4a_matching.m` (the MCNS analogue of the FAFB
   `s16_Mi1_Tm3_T4a_matching.m`). It maps each neuprint Mi1 to its FlyWire match (step 1) and
   inherits that column's (theta, phi) from the FlyWire Mi1 reference; `Tm3_Columns` / `T4a_Columns`
   inherit the synapse-weighted (theta, phi) of their upstream Mi1 columns. Each column also gets its
   denoised medulla / lobula / lobula-plate output synapses. Saves `Mi1_Columns`, `Tm3_Columns`,
   `T4a_Columns` to `MCNS/Processed_Data/Mi1_Tm3_T4a_columns.mat` (`-v7.3`).

3. **Fit the reference-layer surfaces** — run `MCNS/Data_Processing/s17_making_layers.m` (the MCNS
   analogue of the FAFB `s17_making_layers.m`). For each reference cell type (Me = Dm6 dendrite,
   Lo = LT1a dendrite, LoP = T5a axon) it collects the reference synapses inside the target neuropil
   (right and left), denoises them, runs PCA, and fits a parabolic surface (poly22 for Me/Lo, poly33
   for LoP). The per-hemisphere mean, PCA coefficients, and surface fit are stored in `PCA_Basis` and
   saved to `MCNS/Processed_Data/neuropil_PCA_basis.mat`. (It also shows the per-config reference-layer
   fit, right and left.)

4. **Render the fields** — run `MCNS/Figures/fig_6A_BLP_RFs_PFs.m` (the MCNS analogue of the FAFB
   `fig_6E_G_BLP_RFs_PFs.m`). It loads `Mi1_Tm3_T4a_columns.mat` (step 2), `neuropil_PCA_basis.mat`
   (step 3), `BLP_neurons_thr0.mat`, and the MCNS connection / synapse-coordinate tables. For every
   BLP_R cell type it maps each input/output synapse to its nearest reference column's (theta, phi)
   (projecting onto the PCA + parabolic surface via `closest_point_poly22`), accumulates
   synapse-weighted 2D histograms over (phi, theta), and Gaussian-smooths them. It produces, per cell
   type, a per-neuron IN (orange = RF) + OUT (blue = PF) heatmap overlay and a type-level overlay of
   per-neuron contours. The per-type result tables
   (`all_BLP_by_type`, including the RF/PF heatmaps) are saved to
   `MCNS/Processed_Data/BLP_RFs_PFs.mat` (`-v7.3`).

5. **Compute the area/center metrics** — run `MCNS/Data_Processing/s18_BLP_RFs_PFs_area_center.m`
   (the MCNS analogue of the FAFB `s18_BLP_RFs_PFs_area_center.m`). From each neuron's normalized RF
   (input) and PF (output) heatmap it measures, at a contour threshold, the field size (deg²) and the
   weighted center (phi circular, theta linear), and the RF→PF center offset (total distance plus
   signed phi/theta components), including the phi-mirrored RF center. The metric columns are appended
   to `all_BLP_by_type` and saved (with a small `all_BLP_by_type_meta`) to
   `MCNS/Processed_Data/BLP_RFs_PFs_area_center.mat` (`-v7.3`).

6. **Render the area/center plots** — run `MCNS/Figures/fig_6B_BLP_area_center_plot.m` (the MCNS
   analogue of the FAFB `fig_6F_H_J_BLP_area_center_plot.m`). It loads `BLP_RFs_PFs_area_center.mat`
   (step 5) and produces per-type grouped bars of RF vs PF area, RF→PF center offset and its
   phi-mirrored version, plus a per-type mean RF area vs PF area scatter with a 1:1 line and
   Pearson/Spearman correlations; a per-type mean RF/PF area table is printed to the console.

```
MCNS/
├── MCNS_Data/         # neuprint exports (generated by the scripts below)
│   ├── male-cns-v0.9-all-connections.csv
│   ├── male-cns-v0.9-classification.csv
│   ├── male-cns-v0.9-primary-types.csv
│   └── male-cns-v0.9-synapse-coordinates.csv
├── Data_Processing/   # data export / processing scripts
│   ├── s01_export_connections.py          # neuprint -> all-connections CSV
│   ├── s02_export_classification.py       # neuprint -> classification CSV
│   ├── s03_export_primary_types.py        # neuprint -> primary-types CSV
│   ├── s04_export_synapse_coordinates.py  # neuprint -> synapse-coordinates CSV
│   ├── s05_compute_postPI_prePI_all_neurons.m       # per-neuron synapse / NPI table
│   ├── s06_compute_FFP_optic_central_fan_in_out.m   # per-type fan-in/out (Optic/FFP/Central)
│   ├── s07_FFP_ADJ_matrix.m                         # FFP -> central-brain output adjacency matrix
│   ├── s08_FFP_leiden_clustering.ipynb              # bipartite Leiden clustering of the matrix
│   ├── s09_FBP_out_opticlobes.m                     # per-FBP-type output neuropil fractions
│   ├── s10_export_roi_mesh_csv.py                   # neuprint -> optic-lobe ROI meshes (CSV)
│   ├── s11_FBP_output_layer.m                       # FBP output-synapse depth per optic lobe
│   ├── s12_FBP_upstream_layer.m                     # FBP upstream-synapse depth per optic lobe
│   ├── s13_segregation_index.ipynb                  # per-neuron axon/dendrite segregation index
│   ├── s14_BDP_reciprocal.m                         # FFP/FBP/BDP/Others input-output reciprocity (WJI)
│   ├── s15_FAFB_MCNS_Mi1_matching.ipynb             # neuprint Mi1 -> FlyWire Mi1 skeleton matches
│   ├── s16_Mi1_Tm3_T4a_matching.m                  # Mi1/Tm3/T4a reference columns (theta/phi)
│   ├── s17_making_layers.m                          # reference-layer PCA bases (neuropil_PCA_basis.mat)
│   └── s18_BLP_RFs_PFs_area_center.m               # BLP RF/PF area & center metrics
├── Processed_Data/    # generated data (created by the analysis)
│   └── roi_mesh_csv/   # optic-lobe ROI meshes (exported by s10_export_roi_mesh_csv.py)
└── Figures/           # figure-generating scripts
    ├── fig_1_postPI_prePI_right_neurons.m  # postPI/prePI + FFP/FBP/BDP classification
    ├── fig_2_FFP_optic_central.m           # Optic/FFP/Central fan-in/out box plots
    ├── fig_3_FFP_clustering_plot.m         # FFP Leiden clustering dendrogram + connectivity
    ├── fig_4A_FBP_dysynaptic_optclobe.m    # FBP disynaptic optic-lobe input
    ├── fig_4B_FBP_layer_overall_histogram.m # FBP depth-profile mean ± SD histograms
    ├── fig_4C_FBP_layer_overlap.m          # FBP input/output depth overlap
    ├── fig_4S_FBP_layer_figure.m           # FBP output/upstream depth RGB images
    ├── fig_5A_segregation_reciprocal.m     # segregation index vs reciprocity + box plot
    ├── fig_5B_LC9_interconnection.m        # LC9–LC9 interconnection vs proximity
    ├── fig_6A_BLP_RFs_PFs.m                # BLP RF/PF heatmaps in visual space
    └── fig_6B_BLP_area_center_plot.m       # BLP RF/PF area & center plots
```
