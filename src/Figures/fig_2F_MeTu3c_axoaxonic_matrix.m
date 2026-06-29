%% fig_2F_MeTu3c_axoaxonic_matrix
% Figure 2F: the same-type (axo-axonic) central-brain interconnection matrix for a
% single FFP cell type (MeTu3c). For the neurons of one type in the right hemisphere
% it builds the directional same-type synapse matrix split into optic-lobe and
% central-brain parts, clusters the neurons by their combined connectivity profile,
% and renders the central-brain matrix in that clustered leaf order.
%
%   - target neurons = RightFFP_NPIs with type == target_type (right hemisphere)
%   - M(i,j) = directional same-type synapse count, i = pre, j = post
%   - inter_Optic   : synapses in the right optic lobe
%   - inter_Central : synapses in the central brain (neither optic lobe nor UNASGD)
%   - leaf order is computed from the combined matrix [interConnection interConnection']
%     (interConnection = inter_Optic + inter_Central) so rows/cols share one ordering.
%
% Requires Processed_Data/right_neurons_thr0.mat (RightFFP_NPIs), produced by
% Figures/fig_1D_E_postPI_prePI_right_neurons.m, and Codex_Data/connections_no_threshold.csv.
clear all; close all; clc;

baseDir = fileparts(fileparts(mfilename('fullpath')));  % portable; no absolute paths

%% ===== Settings =====
target_type = 'MeTu3c';   % FFP type to display

%% ===== Connection table (no threshold) =====
opt = detectImportOptions(fullfile(baseDir, 'Codex_Data', 'connections_no_threshold.csv'));
opt = setvartype(opt,'pre_root_id','int64');
opt = setvartype(opt,'post_root_id','int64');
FAFBConnections = readtable(fullfile(baseDir, 'Codex_Data', 'connections_no_threshold.csv'),opt);

% Neuropil groups: central = neither optic lobe (L/R) nor UNASGD
FAFBNeuropils = unique(FAFBConnections.neuropil);
FAFBNeuropil_OpticLobeRight = {'LA_R';'AME_R';'ME_R';'LO_R';'LOP_R'};
FAFBNeuropil_OpticLobeLeft  = {'LA_L';'AME_L';'ME_L';'LO_L';'LOP_L'};
FAFBNeuropil_Central = FAFBNeuropils(~ismember(FAFBNeuropils,FAFBNeuropil_OpticLobeRight));
FAFBNeuropil_Central = FAFBNeuropil_Central(~ismember(FAFBNeuropil_Central,FAFBNeuropil_OpticLobeLeft));
FAFBNeuropil_Central = FAFBNeuropil_Central(~ismember(FAFBNeuropil_Central,'UNASGD'));

%% ===== Target neurons : RightFFP type == target_type =====
load(fullfile(baseDir, 'Processed_Data', 'right_neurons_thr0.mat'), 'RightFFP_NPIs')
RightFF_NPIs = RightFFP_NPIs;    % alias (this repo names the table RightFFP_NPIs)

target_root_ids = RightFF_NPIs.root_id(strcmp(RightFF_NPIs.type, target_type));
N = numel(target_root_ids);
fprintf('target type = %s | N = %d\n', target_type, N);
if N < 2
    error('type "%s" has fewer than 2 neurons (N=%d); cannot build interconnection matrix.', target_type, N);
end

%% ===== Same-type (self-to-self) connection matrix, split by neuropil =====
[~, loc_pre]  = ismember(FAFBConnections.pre_root_id,  target_root_ids);
[~, loc_post] = ismember(FAFBConnections.post_root_id, target_root_ids);

% keep only connections whose pre AND post are both in the target type
is_self_to_self = (loc_pre > 0) & (loc_post > 0);
self_conns = FAFBConnections(is_self_to_self, :);
row_idx = loc_pre(is_self_to_self);
col_idx = loc_post(is_self_to_self);
syn_counts = self_conns.syn_count;

is_optic   = ismember(self_conns.neuropil, FAFBNeuropil_OpticLobeRight);
is_central = ismember(self_conns.neuropil, FAFBNeuropil_Central);

inter_Optic     = accumarray([row_idx(is_optic),   col_idx(is_optic)],   syn_counts(is_optic),   [N, N]);
inter_Central   = accumarray([row_idx(is_central), col_idx(is_central)], syn_counts(is_central), [N, N]);
interConnection = inter_Optic + inter_Central;

% mean out-degree (avg number of distinct partners a neuron sends to), central only
meandeg = @(M) mean(sum( (M>0) & ~eye(N), 2));
deg_central = meandeg(inter_Central);

%% ===== Hierarchical-clustering leaf order (shared rows/cols) =====
X = [interConnection interConnection'];   % N x 2N : combine out- and in-profiles
Y = pdist(X, 'euclidean');
Z = linkage(Y, 'average');
order = optimalleaforder(Z, Y);

%% ===== Figure 1 : central-brain interconnection matrix (clustered order) =====
figure(1); set(gcf,'Color','w')
imagesc(inter_Central(order, order))
colormap(flipud(gray));   % reversed gray (0 = white, large = black)
set(gca,'Box','off','TickDir','out','YDir','normal');
axis('square');
clim([0 10])
title(sprintf('%s interConnection central | avg %.1f out-targets/neuron', target_type, deg_central));
