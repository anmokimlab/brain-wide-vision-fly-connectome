%% s16 — Build Mi1 / Tm3 / T4a reference columns (visual-space theta/phi)  [MCNS]
% MCNS analogue of the FAFB Data_Processing/s16_Mi1_Tm3_T4a_matching.m, combining
% the MCNS Mi1_columns.m / Mi1_Tm3_Matching.m / Mi1_T4a_Matching.m scripts.
%
% - Mi1_Columns: anchor columns. Each neuprint Mi1 is matched to its FlyWire Mi1
%   (s15 directed-mean-nearest skeleton matches) and inherits that column's
%   (theta, phi) from the FlyWire Mi1 reference; medulla output synapses are kept
%   and denoised.
% - Tm3_Columns / T4a_Columns: each column inherits the (theta, phi) of its upstream
%   Mi1 columns (synapse-weighted), plus its denoised lobula / lobula-plate output
%   synapses.
%
% Saves Mi1_Columns, Tm3_Columns, T4a_Columns to
% Processed_Data/Mi1_Tm3_T4a_columns.mat (-v7.3).
%
% Note on voxels/skeletons: unlike the FAFB s15 (which loads SWC skeletons to assign
% synapses to a neuropil via mesh intriangulation and to trace SWC ancestors of the
% output synapses), the MCNS synapse table carries a `neuropil` column, so synapses
% are assigned to a region by that column and no skeleton voxels / ancestors are used.

clear all; close all; clc;

baseDir = fileparts(fileparts(mfilename('fullpath')));  % resolved from this script's location (portable; no absolute paths)

%% Load shared data
opt = detectImportOptions(fullfile(baseDir, 'MCNS_Data', 'male-cns-v0.9-primary-types.csv'));
opt = setvartype(opt,'root_id','int64');
MCNSConsolidatedTypes = readtable(fullfile(baseDir, 'MCNS_Data', 'male-cns-v0.9-primary-types.csv'),opt);

opt = detectImportOptions(fullfile(baseDir, 'MCNS_Data', 'male-cns-v0.9-all-connections.csv'));
opt = setvartype(opt,'pre_root_id','int64');
opt = setvartype(opt,'post_root_id','int64');
MCNSConnections = readtable(fullfile(baseDir, 'MCNS_Data', 'male-cns-v0.9-all-connections.csv'),opt);

% Synapse coordinates (8 nm voxels -> metres; columns 1:3 = pre_x/pre_y/pre_z)
opt = detectImportOptions(fullfile(baseDir, 'MCNS_Data', 'male-cns-v0.9-synapse-coordinates.csv'));
opt = setvartype(opt,'pre_root_id','int64');
opt = setvartype(opt,'post_root_id','int64');
MCNS_synapse_coordinates = readtable(fullfile(baseDir, 'MCNS_Data', 'male-cns-v0.9-synapse-coordinates.csv'),opt);
MCNS_synapse_coordinates.pre_x = MCNS_synapse_coordinates.pre_x*8e-9;
MCNS_synapse_coordinates.pre_y = MCNS_synapse_coordinates.pre_y*8e-9;
MCNS_synapse_coordinates.pre_z = MCNS_synapse_coordinates.pre_z*8e-9;

% FlyWire Mi1 reference columns (mi1_root_id, theta, phi), from Garner et al. 2024,
% Nature (https://www.nature.com/articles/s41586-024-07967-z). Shared with the FAFB
% Reference_Data at the project root.
opt = detectImportOptions(fullfile(baseDir, '..', 'Reference_Data', 'Mi1_columns_DustinGarner.csv'));
opt = setvartype(opt,'mi1_root_id','int64');
Mi1_flywire = readtable(fullfile(baseDir, '..', 'Reference_Data', 'Mi1_columns_DustinGarner.csv'),opt);
Mi1_flywire = rmmissing(Mi1_flywire);

% neuprint Mi1 -> FlyWire Mi1 matches, produced by
% Data_Processing/s15_FAFB_MCNS_Mi1_matching.ipynb
opt = detectImportOptions(fullfile(baseDir, 'Processed_Data', 'neuprint-mi1-to-flywire-v783-directed-mean-nearest-skeleton-matches.csv'));
opt = setvartype(opt,'neuprint_body_id','int64');
opt = setvartype(opt,'flywire_root_id','int64');
Mi1_Matching = readtable(fullfile(baseDir, 'Processed_Data', 'neuprint-mi1-to-flywire-v783-directed-mean-nearest-skeleton-matches.csv'),opt);

%% ====================== Mi1 columns ======================
% Each neuprint Mi1 inherits the (theta, phi) of its matched FlyWire Mi1 column.
Mi1_Columns = table(Mi1_Matching.neuprint_body_id, 'VariableNames', {'root_id'});
Mi1_Columns.hemisphere = Mi1_Matching.neuprint_side_raw;
Mi1_Columns.theta = nan(height(Mi1_Columns),1);
Mi1_Columns.phi   = nan(height(Mi1_Columns),1);
for i = 1:height(Mi1_Columns)
    idx = Mi1_flywire.mi1_root_id == Mi1_Matching.flywire_root_id(i);
    if any(idx)
        Mi1_Columns.theta(i) = Mi1_flywire.theta(idx);
        Mi1_Columns.phi(i)   = Mi1_flywire.phi(idx);
    end
end
Mi1_Columns = rmmissing(Mi1_Columns);

% Output synapses inside the medulla, then denoise
Mi1_Columns = attach_out_syn_in_region(Mi1_Columns, MCNS_synapse_coordinates, {'ME(R)','ME(L)'});
Mi1_Columns = denoise_out_syn(Mi1_Columns);

%% ====================== Tm3 columns ======================
Tm3_Columns = table(MCNSConsolidatedTypes.root_id(strcmp(MCNSConsolidatedTypes.primary_type,'Tm3')), ...
    'VariableNames', {'root_id'});

% Inherit (theta, phi) from upstream Mi1 columns (synapse-weighted)
Tm3_Columns = match_upstream_mi1(Tm3_Columns, MCNSConnections, Mi1_Columns);

% Output synapses inside the lobula, then denoise
Tm3_Columns = attach_out_syn_in_region(Tm3_Columns, MCNS_synapse_coordinates, {'LO(R)','LO(L)'});
Tm3_Columns = denoise_out_syn(Tm3_Columns);

%% ====================== T4a columns ======================
T4a_Columns = table(MCNSConsolidatedTypes.root_id(strcmp(MCNSConsolidatedTypes.primary_type,'T4a')), ...
    'VariableNames', {'root_id'});

T4a_Columns = match_upstream_mi1(T4a_Columns, MCNSConnections, Mi1_Columns);

% Output synapses inside the lobula plate, then denoise
T4a_Columns = attach_out_syn_in_region(T4a_Columns, MCNS_synapse_coordinates, {'LOP(R)','LOP(L)'});
T4a_Columns = denoise_out_syn(T4a_Columns);

%% Save
save(fullfile(baseDir, 'Processed_Data', 'Mi1_Tm3_T4a_columns.mat'), ...
    'Mi1_Columns', 'Tm3_Columns', 'T4a_Columns', '-v7.3');

%% ====================== Local functions ======================
function T = match_upstream_mi1(T, MCNSConnections, Mi1_Columns)
% For each row, collect upstream Mi1 partners and inherit their (theta, phi),
% weighted by synapse count. Drops rows with no upstream Mi1.
zeroIdx = [];
for i = 1:size(T,1)
    upstreamConnections = MCNSConnections(MCNSConnections.post_root_id==T.root_id(i),:);
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

function T = attach_out_syn_in_region(T, MCNS_synapse_coordinates, neuropilNames)
% Keep each neuron's output (pre) synapses whose `neuropil` is in neuropilNames.
% out_syn_loc is stored as an [x y z] numeric array. Drops rows with none.
zeroIdx = [];
for i = 1:size(T,1)
    out_syn_idx = find(ismember(MCNS_synapse_coordinates.pre_root_id, T.root_id(i)));
    out_syn_loc = MCNS_synapse_coordinates(out_syn_idx, :);
    out_syn_loc = out_syn_loc(ismember(out_syn_loc.neuropil, neuropilNames), :);
    out_syn_loc = out_syn_loc{:, 1:3};
    if isempty(out_syn_loc)
        zeroIdx = [zeroIdx i]; %#ok<AGROW>
        continue;
    end
    T.out_syn_loc{i} = out_syn_loc;
end
T(zeroIdx,:) = [];
end

function T = denoise_out_syn(T)
% Pool output synapses per hemisphere, denoise the point cloud (pcdenoise), and
% keep each neuron's synapses that survived. Drops rows with none left.
out_syn_loc_R = [];
out_syn_loc_L = [];
for i = 1:size(T,1)
    if isempty(T.out_syn_loc{i}), continue; end
    if strcmp(T.hemisphere{i},'R')
        out_syn_loc_R = [out_syn_loc_R; T.out_syn_loc{i}]; %#ok<AGROW>
    elseif strcmp(T.hemisphere{i},'L')
        out_syn_loc_L = [out_syn_loc_L; T.out_syn_loc{i}]; %#ok<AGROW>
    end
end

R_denoise = denoise_cloud(out_syn_loc_R);
L_denoise = denoise_cloud(out_syn_loc_L);

zeroIdx = [];
for i = 1:size(T,1)
    if isempty(T.out_syn_loc{i}), continue; end
    arr = T.out_syn_loc{i};
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
