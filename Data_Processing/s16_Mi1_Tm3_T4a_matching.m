%% s16 — Build Mi1 / Tm3 / T4a reference columns (visual-space theta/phi)
% - Mi1_Columns: anchor columns (theta/phi per column come from the Mi1 reference CSV),
%   with denoised medulla output-synapse locations.
% - Tm3_Columns / T4a_Columns: each Tm3/T4a column inherits the (theta,phi) of its
%   upstream Mi1 columns (synapse-weighted), plus its denoised lobula / lobula-plate
%   output-synapse locations.
%
% Saves Mi1_Columns, Tm3_Columns, T4a_Columns to
% Processed_Data/Mi1_Tm3_T4a_columns.mat (-v7.3), read by
% Figures/fig_6E_G_BLP_RFs_PFs.m.
%
% Requires the intriangulation helper (Helper_Function/) on the MATLAB path.

clear all; close all; clc;

baseDir = fileparts(fileparts(mfilename('fullpath')));  % resolved from this script's location (portable; no absolute paths)

%% Load shared data
typeOpts = detectImportOptions(fullfile(baseDir, 'Codex_Data', 'consolidated_cell_types.csv'));
typeOpts = setvartype(typeOpts, 'root_id', 'int64');
FAFBConsolidated_type = readtable(fullfile(baseDir, 'Codex_Data', 'consolidated_cell_types.csv'), typeOpts);

connOpts = detectImportOptions(fullfile(baseDir, 'Codex_Data', 'connections_no_threshold.csv'));
connOpts = setvartype(connOpts, 'pre_root_id', 'int64');
connOpts = setvartype(connOpts, 'post_root_id', 'int64');
FAFBConnections = readtable(fullfile(baseDir, 'Codex_Data', 'connections_no_threshold.csv'), connOpts);

synOpts = detectImportOptions(fullfile(baseDir, 'Codex_Data', 'fafb_v783_princeton_synapse_table.csv'));
synOpts = setvartype(synOpts, 'pre_root_id_720575940', 'int64');
synOpts = setvartype(synOpts, 'post_root_id_720575940', 'int64');
FAFB_synapse_coordinates = readtable(fullfile(baseDir, 'Codex_Data', 'fafb_v783_princeton_synapse_table.csv'), synOpts);
FAFB_synapse_coordinates.Properties.VariableNames(11) = "pre_root_id";
FAFB_synapse_coordinates.Properties.VariableNames(12) = "post_root_id";
FAFB_synapse_coordinates.pre_root_id  = FAFB_synapse_coordinates.pre_root_id  + int64(720575940000000000);
FAFB_synapse_coordinates.post_root_id = FAFB_synapse_coordinates.post_root_id + int64(720575940000000000);

% Optic-lobe neuropil meshes (downloaded by s01_fetch_neuropil_neuron.ipynb)
meshDir = fullfile(baseDir, 'Processed_Data', 'optic_lobe_neuropil_mesh');
[Me_R_Faces,  Me_R_Vertices]  = load_mesh(meshDir, 'Me_R');
[Me_L_Faces,  Me_L_Vertices]  = load_mesh(meshDir, 'Me_L');
[Lo_R_Faces,  Lo_R_Vertices]  = load_mesh(meshDir, 'Lo_R');
[Lo_L_Faces,  Lo_L_Vertices]  = load_mesh(meshDir, 'Lo_L');
[LoP_R_Faces, LoP_R_Vertices] = load_mesh(meshDir, 'LoP_R');
[LoP_L_Faces, LoP_L_Vertices] = load_mesh(meshDir, 'LoP_L');

% Mi1 reference columns (theta/phi per column).
% Source: Garner et al. 2024, Nature (https://www.nature.com/articles/s41586-024-07967-z)
mi1Opts = detectImportOptions(fullfile(baseDir, 'Reference_Data', 'Mi1_columns_DustinGarner.csv'));
mi1Opts = setvartype(mi1Opts, 'mi1_root_id', 'int64');
Mi1_Columns = readtable(fullfile(baseDir, 'Reference_Data', 'Mi1_columns_DustinGarner.csv'), mi1Opts);
Mi1_Columns = rmmissing(Mi1_Columns);
Mi1_Columns.Properties.VariableNames(1) = "root_id";

%% ====================== Mi1 columns ======================
% Output synapses inside the medulla, then denoise
Mi1_Columns = attach_out_syn_in_region(Mi1_Columns, FAFB_synapse_coordinates, ...
    {Me_R_Vertices, Me_R_Faces; Me_L_Vertices, Me_L_Faces});
Mi1_Columns = denoise_out_syn(Mi1_Columns);

%% ====================== Tm3 columns ======================
Tm3_Columns = table(FAFBConsolidated_type.root_id(strcmp(FAFBConsolidated_type.primary_type,'Tm3')), ...
    'VariableNames', {'root_id'});

% Inherit (theta,phi) from upstream Mi1 columns (synapse-weighted)
Tm3_Columns = match_upstream_mi1(Tm3_Columns, FAFBConnections, Mi1_Columns);

% Output synapses inside the lobula, then denoise
Tm3_Columns = attach_out_syn_in_region(Tm3_Columns, FAFB_synapse_coordinates, ...
    {Lo_R_Vertices, Lo_R_Faces; Lo_L_Vertices, Lo_L_Faces});
Tm3_Columns = denoise_out_syn(Tm3_Columns);

%% ====================== T4a columns ======================
T4a_Columns = table(FAFBConsolidated_type.root_id(strcmp(FAFBConsolidated_type.primary_type,'T4a')), ...
    'VariableNames', {'root_id'});

T4a_Columns = match_upstream_mi1(T4a_Columns, FAFBConnections, Mi1_Columns);

% Output synapses inside the lobula plate, then denoise
T4a_Columns = attach_out_syn_in_region(T4a_Columns, FAFB_synapse_coordinates, ...
    {LoP_R_Vertices, LoP_R_Faces; LoP_L_Vertices, LoP_L_Faces});
T4a_Columns = denoise_out_syn(T4a_Columns);

%% Save
save(fullfile(baseDir, 'Processed_Data', 'Mi1_Tm3_T4a_columns.mat'), ...
    'Mi1_Columns', 'Tm3_Columns', 'T4a_Columns', '-v7.3');

%% ====================== Local functions ======================
function [F, V] = load_mesh(meshDir, prefix)
% Read a neuropil mesh; faces are 0-based, corrected to 1-based.
F = readmatrix(fullfile(meshDir, [prefix '_faces.csv'])) + 1;
V = readmatrix(fullfile(meshDir, [prefix '_vertices.csv']));
end

function T = match_upstream_mi1(T, FAFBConnections, Mi1_Columns)
% For each row, collect upstream Mi1 partners and inherit their (theta,phi),
% weighted by synapse count. Drops rows with no upstream Mi1.
zeroIdx = [];
for i = 1:size(T,1)
    upstreamConnections = FAFBConnections(FAFBConnections.post_root_id==T.root_id(i),:);
    upstreamConnectionsMi1 = upstreamConnections(ismember(upstreamConnections.pre_root_id, Mi1_Columns.root_id),:);
    upstreamMi1 = table(unique(upstreamConnectionsMi1.pre_root_id), 'VariableNames', {'root_id'});
    if isempty(upstreamMi1)
        zeroIdx = [zeroIdx i]; %#ok<AGROW>
        continue;
    end
    for j = 1:size(upstreamMi1,1)
        idx_connection = upstreamConnectionsMi1.pre_root_id==upstreamMi1.root_id(j);
        upstreamMi1.syn_count(j) = sum(upstreamConnectionsMi1.syn_count(idx_connection));

        idx_column = Mi1_Columns.root_id==upstreamMi1.root_id(j);
        upstreamMi1.theta(j)      = Mi1_Columns.theta(idx_column);
        upstreamMi1.phi(j)        = Mi1_Columns.phi(idx_column);
        upstreamMi1.hemisphere(j) = Mi1_Columns.hemisphere(idx_column);
    end
    T.upstreamMi1{i}  = upstreamMi1;
    T.hemisphere(i)   = upstreamMi1.hemisphere(1);
    T.theta{i}        = [upstreamMi1.theta upstreamMi1.syn_count/sum(upstreamMi1.syn_count)];
    T.phi{i}          = [upstreamMi1.phi   upstreamMi1.syn_count/sum(upstreamMi1.syn_count)];
end
T(zeroIdx,:) = [];
end

function T = attach_out_syn_in_region(T, FAFB_synapse_coordinates, regions)
% Keep each neuron's output (pre) synapses that fall inside the region meshes.
% Drops rows with no output synapses inside.
zeroIdx = [];
for i = 1:size(T,1)
    out_syn_idx = find(ismember(FAFB_synapse_coordinates.pre_root_id, T.root_id(i)));
    out_syn_loc = FAFB_synapse_coordinates(out_syn_idx, 1:3);
    if isempty(out_syn_loc)
        zeroIdx = [zeroIdx i]; %#ok<AGROW>
        continue;
    end
    arr = table2array(out_syn_loc);
    inside = false(size(arr,1),1);
    for r = 1:size(regions,1)
        inside = inside | intriangulation(regions{r,1}, regions{r,2}, arr);
    end
    out_syn_loc = out_syn_loc(inside, :);
    if isempty(out_syn_loc)
        zeroIdx = [zeroIdx i]; %#ok<AGROW>
        continue;
    end
    T.out_syn_loc{i} = out_syn_loc;
end
T(zeroIdx,:) = [];
end

function T = denoise_out_syn(T)
% Pool output synapses per hemisphere, denoise the point cloud (pcdenoise),
% and keep each neuron's synapses that survived. Drops rows with none left.
out_syn_loc_R = [];
out_syn_loc_L = [];
for i = 1:size(T,1)
    if isempty(T.out_syn_loc{i}), continue; end
    if strcmp(T.hemisphere{i},'R')
        out_syn_loc_R = [out_syn_loc_R; table2array(T.out_syn_loc{i})]; %#ok<AGROW>
    elseif strcmp(T.hemisphere{i},'L')
        out_syn_loc_L = [out_syn_loc_L; table2array(T.out_syn_loc{i})]; %#ok<AGROW>
    end
end

R_denoise = denoise_cloud(out_syn_loc_R);
L_denoise = denoise_cloud(out_syn_loc_L);

zeroIdx = [];
for i = 1:size(T,1)
    if isempty(T.out_syn_loc{i}), continue; end
    arr = table2array(T.out_syn_loc{i});
    idx = ismember(arr, R_denoise, 'rows') | ismember(arr, L_denoise, 'rows');
    if sum(idx)==0
        zeroIdx = [zeroIdx i]; %#ok<AGROW>
        continue;
    end
    T.out_syn_loc_denoise{i} = array2table(arr(idx,:), 'VariableNames', {'x','y','z'});
end
T(zeroIdx,:) = [];
end

function loc = denoise_cloud(pts)
% pcdenoise wrapper; returns [] for an empty input.
if isempty(pts)
    loc = [];
    return
end
pc = pcdenoise(pointCloud(pts), 'NumNeighbors', 20, 'Threshold', 1);
loc = pc.Location;
end
