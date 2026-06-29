%% 1. Load data and initialize
clear all; close all; clc

baseDir = fileparts(fileparts(mfilename('fullpath')));  % resolved from this script's location (portable; no absolute paths)

% Codex connectivity (force root_id columns to int64)
opt = detectImportOptions(fullfile(baseDir, 'Codex_Data', 'connections_no_threshold.csv'));
opt = setvartype(opt, 'pre_root_id', 'int64');
opt = setvartype(opt, 'post_root_id', 'int64');
FAFBConnections = readtable(fullfile(baseDir, 'Codex_Data', 'connections_no_threshold.csv'), opt);

% Neuropil region definitions. Central = everything but the optic lobes and UNASGD.
FAFBNeuropils = unique(FAFBConnections.neuropil);
FAFBNeuropil_OpticLobeRight = {'LA_R';'AME_R';'ME_R';'LO_R';'LOP_R'};
FAFBNeuropil_OpticLobeLeft  = {'LA_L';'AME_L';'ME_L';'LO_L';'LOP_L'};
FAFBNeuropil_Central = FAFBNeuropils(~ismember(FAFBNeuropils, FAFBNeuropil_OpticLobeRight));
FAFBNeuropil_Central = FAFBNeuropil_Central(~ismember(FAFBNeuropil_Central, FAFBNeuropil_OpticLobeLeft));
FAFBNeuropil_Central = FAFBNeuropil_Central(~ismember(FAFBNeuropil_Central, 'UNASGD'));

% Synapse coordinates (force root_id columns to int64)
opt = detectImportOptions(fullfile(baseDir, 'Codex_Data', 'fafb_v783_princeton_synapse_table.csv'));
opt = setvartype(opt, 'pre_root_id_720575940', 'int64');
opt = setvartype(opt, 'post_root_id_720575940', 'int64');
FAFB_synapse_coordinates = readtable(fullfile(baseDir, 'Codex_Data', 'fafb_v783_princeton_synapse_table.csv'), opt);
FAFB_synapse_coordinates.Properties.VariableNames(11) = "pre_root_id";
FAFB_synapse_coordinates.Properties.VariableNames(12) = "post_root_id";
FAFB_synapse_coordinates.pre_root_id  = FAFB_synapse_coordinates.pre_root_id  + int64(720575940000000000);
FAFB_synapse_coordinates.post_root_id = FAFB_synapse_coordinates.post_root_id + int64(720575940000000000);

% Right-hemisphere LC9 (real-bidirectional) neurons, from
% Figures/fig_1D_E_postPI_prePI_right_neurons.m
load(fullfile(baseDir, 'Processed_Data', 'right_neurons_thr0.mat'), 'RightBDP_real_NPIs')
LC9_R_root_ids = RightBDP_real_NPIs.root_id(strcmp(RightBDP_real_NPIs.type, 'LC9'));

%% 2. LC9 -> LC9 central interconnection matrix
% Map each LC9 root_id to a 1..N matrix index, keep only LC9->LC9 connections in
% the central brain, and accumulate synapse counts into an N x N matrix.
[~, loc_pre]  = ismember(FAFBConnections.pre_root_id,  LC9_R_root_ids);
[~, loc_post] = ismember(FAFBConnections.post_root_id, LC9_R_root_ids);
is_LC9_to_LC9 = (loc_pre > 0) & (loc_post > 0);

LC9_conns = FAFBConnections(is_LC9_to_LC9, :);
row_idx = loc_pre(is_LC9_to_LC9);
col_idx = loc_post(is_LC9_to_LC9);
is_central = ismember(LC9_conns.neuropil, FAFBNeuropil_Central);

N = numel(LC9_R_root_ids);
LC9_inter_Central = accumarray([row_idx(is_central), col_idx(is_central)], ...
    LC9_conns.syn_count(is_central), [N, N]);

%% 3. LC9 synapse centroid inside the right optic lobe
% For each LC9 neuron, take its input synapses (the neuron is post) falling inside
% the right optic lobe and average their coordinates into a single centroid point.
% Optic_R is the combined whole right-optic-lobe mesh built by
% Data_Processing/s01_fetch_neuropil_neuron.ipynb.
% Requires the intriangulation helper (Helper_Function/) on the MATLAB path.
meshDir = fullfile(baseDir, 'Processed_Data', 'optic_lobe_neuropil_mesh');
Optic_R_Faces    = readmatrix(fullfile(meshDir, 'Optic_R_faces.csv'));
Optic_R_Vertices = readmatrix(fullfile(meshDir, 'Optic_R_vertices.csv'));
Optic_R_Faces = Optic_R_Faces + 1;   % MATLAB indexing correction

LC9_centroid_Optic_R = nan(N, 3);   % one OL centroid per neuron (NaN if no OL synapses)
for i = 1:N
    Neuron_root_id = LC9_R_root_ids(i);

    in_xyz  = [FAFB_synapse_coordinates.post_x(FAFB_synapse_coordinates.post_root_id == Neuron_root_id), ...
               FAFB_synapse_coordinates.post_y(FAFB_synapse_coordinates.post_root_id == Neuron_root_id), ...
               FAFB_synapse_coordinates.post_z(FAFB_synapse_coordinates.post_root_id == Neuron_root_id)];
    syn_xyz = in_xyz;   % input synapses only (neuron is post)

    syn_in_OL = syn_xyz(intriangulation(Optic_R_Vertices, Optic_R_Faces, syn_xyz), :);
    if ~isempty(syn_in_OL)
        LC9_centroid_Optic_R(i, :) = mean(syn_in_OL, 1);   % centroid of OL synapses
    end
end

%% 4. Pairwise distance between optic-lobe synapse centroids
% Distance between two LC9 neurons = Euclidean distance between their OL centroids.
dist_matrix_avg = zeros(N);
for i = 1:N
    if any(isnan(LC9_centroid_Optic_R(i, :))), continue; end
    for j = i+1:N
        if any(isnan(LC9_centroid_Optic_R(j, :))), continue; end
        d = norm(LC9_centroid_Optic_R(i, :) - LC9_centroid_Optic_R(j, :));
        dist_matrix_avg(i, j) = d;
        dist_matrix_avg(j, i) = d;
    end
end

%% 5. Connection weight vs distance vectors (directed; both off-diagonal directions, no self-pairs)
% Keep the central interconnection weight directed (W(i,j) = synapses i -> j) instead
% of symmetrising it. Each unordered pair therefore contributes two entries, i->j and
% j->i, which share the same centroid distance but may carry different weights.
W_dir = LC9_inter_Central;
W_dir(logical(eye(N))) = 0;
dist_matrix_avg(logical(eye(N))) = 0;

mask = ~logical(eye(N));   % all off-diagonal entries (both directions)
W_vector = W_dir(mask);
D_vector = dist_matrix_avg(mask);

%% 6. Figure 1 (panel 5I): inter-synaptic distance, connected vs unconnected
connected_dists_5 = D_vector(W_vector >= 5);
connected_dists   = D_vector(W_vector > 0);
unconnected_dists = D_vector(W_vector == 0);

p_ranksum = ranksum(connected_dists, unconnected_dists);

figure(1); set(gcf,'Color','w')
group = [zeros(size(connected_dists_5)); ones(size(connected_dists)); 2*ones(size(unconnected_dists))];
boxplot([connected_dists_5; connected_dists; unconnected_dists], group, ...
    'Labels', {'Connected>=5','Connected','Unconnected'}, 'Symbol','');
ylabel('Inter-centroid Distance (nm)');
title(['Distance Comparison (p = ', num2str(p_ranksum, '%.4e'), ')']);
set(gca,'Box','off','TickDir','out')
ylim([0 160000])

%% 7. Figure 2 (panel 5J): connection probability vs distance, with pair frequency
binSize = 10000;   % 10 um (coordinates in nm)
maxD = max(D_vector);
bins = 0:binSize:maxD+binSize;

prob_0    = zeros(length(bins)-1, 1);   % P(connected, W > 0)  per distance bin
prob_5    = zeros(length(bins)-1, 1);   % P(strongly connected, W >= 5) per bin
num_pairs = zeros(length(bins)-1, 1);   % number of neuron pairs per bin
center_D  = zeros(length(bins)-1, 1);

for i = 1:length(bins)-1
    idx = (D_vector >= bins(i)) & (D_vector < bins(i+1));
    num_pairs(i) = sum(idx);
    if num_pairs(i) > 0
        prob_0(i) = sum(W_vector(idx) > 0)  / num_pairs(i);
        prob_5(i) = sum(W_vector(idx) >= 5) / num_pairs(i);
    end
    center_D(i) = (bins(i) + bins(i+1)) / 2;
end

figure(2); set(gcf,'Color','w','Name','Connectivity Probability and Frequency')

% Left axis: connection probability
yyaxis left
hBar = bar(center_D, [prob_0, prob_5], 'grouped');
hBar(1).FaceColor = [0.7 0.7 0.7]; hBar(1).EdgeColor = 'none';   % W > 0
hBar(2).FaceColor = [0.2 0.4 0.8]; hBar(2).EdgeColor = 'none';   % W >= 5
ylabel('Probability of Connection');
ylim([0 1]);

% Right axis: pair frequency (sample size)
yyaxis right
hLine = plot(center_D, num_pairs, '-s', 'Color', [0.8 0.2 0.2], 'LineWidth', 1.5, 'MarkerSize', 5);
ylabel('Number of Pairs (Frequency)');
ax = gca; ax.YColor = [0.8 0.2 0.2];
ylim([0 550]);

xlabel('Distance between neurons (nm)');
title('Connectivity Probability vs. Distance with Pair Frequency');
legend([hBar(1), hBar(2), hLine], ...
    {'Threshold: W > 0', 'Threshold: W \geq 5', 'Total Pair Count (N)'}, 'Box', 'off');
set(gca, 'Box','off', 'TickDir','out', 'FontSize',12, 'XTick',center_D, 'YTick',0:100:500);
