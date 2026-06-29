%% s06_FFP_axo_axonic_compute
% Computes, from scratch, the data behind Figure 2F (axo-axonic input) plus the
% optic-lobe proximity-PROBABILITY reference analysis. The paper optic-lobe panels
% 2G / 2H / S1D are produced by s07_FFP_axo_axonic_RF.m instead, which computes the
% RF distance (the OL synapse-centroid distance, um) binned by W.
%   A) For every FFP cell type, the axo-axonic input from the same type:
%      the fraction of each neuron's total input synapses that arise from
%      central-brain (CB) connections between neurons of the same type,
%      averaged per type (% of total input).                          -> Fig 2F
%   B) High same-type-input outlier types, flagged by the upper Tukey fence
%      (Q3 + 1.5 x IQR) of the per-type values.
%   C) Optic-lobe (OL) proximity probability P (reference), computed ONLY for the outlier
%      types. Each neuron's OL INPUT synapses (the neuron is post) are reduced to
%      a single centroid, and the OL distance between two neurons is the Euclidean
%      distance between their centroids. P is computed per type in a per-neuron (stratified)
%      manner: for each neuron, the OL distances to its connected partners
%      (directional same-type CB synapse count W >= thr) are compared with the
%      distances to its unconnected partners, pooled across all neurons of the
%      type. Two connectivity thresholds are used (W > 0, W >= 5).
% Figures are drawn by Figures/fig_2F_G_FFP_axo_axonic_plot.m, which reads the
% saved data only. Requires the intriangulation helper (Helper_Function/) on path.
clear all; close all; clc;

baseDir = fileparts(fileparts(mfilename('fullpath')));  % portable; no absolute paths

%% ===== Common data =====
% Connection table (no threshold)
opt = detectImportOptions(fullfile(baseDir, 'Codex_Data', 'connections_no_threshold.csv'));
opt = setvartype(opt,'pre_root_id','int64');
opt = setvartype(opt,'post_root_id','int64');
FAFBConnections = readtable(fullfile(baseDir, 'Codex_Data', 'connections_no_threshold.csv'),opt);

% FFP neuron table (RightFFP_NPIs), from fig_1D_E
load(fullfile(baseDir, 'Processed_Data', 'right_neurons_thr0.mat'), 'RightFFP_NPIs')

% Neuropil grouping
FAFBNeuropils = unique(FAFBConnections.neuropil);
FAFBNeuropil_OpticLobeRight = {'LA_R';'AME_R';'ME_R';'LO_R';'LOP_R'};
FAFBNeuropil_OpticLobeLeft  = {'LA_L';'AME_L';'ME_L';'LO_L';'LOP_L'};
FAFBNeuropil_OpticLobe = [FAFBNeuropil_OpticLobeRight; FAFBNeuropil_OpticLobeLeft];
% Central for the proximity step (UNASGD excluded)
FAFBNeuropil_Central = FAFBNeuropils(~ismember(FAFBNeuropils, FAFBNeuropil_OpticLobe));
FAFBNeuropil_Central = FAFBNeuropil_Central(~ismember(FAFBNeuropil_Central,'UNASGD'));

% Per-neuron cell type (use the table column if present, else map via Codex types)
if ismember('type', RightFFP_NPIs.Properties.VariableNames)
    neuronType = string(RightFFP_NPIs.type);
else
    optT = detectImportOptions(fullfile(baseDir, 'Codex_Data', 'consolidated_cell_types.csv'));
    optT = setvartype(optT,'root_id','int64');
    FAFBTypes = readtable(fullfile(baseDir, 'Codex_Data', 'consolidated_cell_types.csv'),optT);
    [tfT, locT]   = ismember(RightFFP_NPIs.root_id, FAFBTypes.root_id);
    neuronType    = strings(height(RightFFP_NPIs),1);
    neuronType(tfT) = string(FAFBTypes.primary_type(locT(tfT)));
end

%% ===== A. Axo-axonic input from the same type (per type, % of total input) =====
% Same-type CB input = synapses where pre and post share the same type and the
% synapse is outside the optic lobe; expressed as a fraction of total input.
root_ids = RightFFP_NPIs.root_id;
N = numel(root_ids);
[typeCats, ~, typeCode] = unique(neuronType);
nTypes = numel(typeCats);

[~, preLoc]  = ismember(FAFBConnections.pre_root_id,  root_ids);
[~, postLoc] = ismember(FAFBConnections.post_root_id, root_ids);
syn = FAFBConnections.syn_count;

isPostTarget = postLoc > 0;
isOptic      = ismember(FAFBConnections.neuropil, FAFBNeuropil_OpticLobe);
isCentralIn  = ~isOptic;                                  % CB input includes UNASGD here
isBoth = (preLoc > 0) & (postLoc > 0);
isSelf = false(height(FAFBConnections),1);
isSelf(isBoth) = typeCode(preLoc(isBoth)) == typeCode(postLoc(isBoth));

acc = @(m) accumarray(postLoc(m), syn(m), [N 1]);
self_optic    = acc(isPostTarget & isOptic      & isSelf);
self_central  = acc(isPostTarget & isCentralIn  & isSelf);   % axo-axonic same-type input
other_optic   = acc(isPostTarget & isOptic      & ~isSelf);
other_central = acc(isPostTarget & isCentralIn  & ~isSelf);
total_in = self_optic + self_central + other_optic + other_central;

valid = total_in > 0;
pct_aa = nan(N,1);
pct_aa(valid) = self_central(valid) ./ total_in(valid) * 100;   % per-neuron %

axoaxonic_pct = nan(nTypes,1);
n_neurons     = zeros(nTypes,1);
for k = 1:nTypes
    rows = (typeCode == k & valid);
    n_neurons(k) = sum(rows);
    if n_neurons(k) > 0
        axoaxonic_pct(k) = mean(pct_aa(rows));
    end
end

okType        = n_neurons > 0 & ~isnan(axoaxonic_pct);
types         = typeCats(okType);
axoaxonic_pct = axoaxonic_pct(okType);
n_neurons     = n_neurons(okType);
nType         = numel(types);

%% ===== B. Tukey outlier types (upper / lower fence) =====
whisker     = 1.5;
q1          = prctile(axoaxonic_pct, 25);
q3          = prctile(axoaxonic_pct, 75);
iqr_        = q3 - q1;
lower_fence = q1 - whisker*iqr_;
upper_fence = q3 + whisker*iqr_;

is_outlier  = (axoaxonic_pct > upper_fence) | (axoaxonic_pct < lower_fence);
direction   = strings(nType,1);
direction(axoaxonic_pct > upper_fence) = "high";
direction(axoaxonic_pct < lower_fence) = "low";

outTypes       = types(is_outlier);
out_axoaxonic  = axoaxonic_pct(is_outlier);
nOut           = numel(outTypes);
fprintf('Outlier fence [%.2f, %.2f]%% -> %d outlier types / %d total\n', ...
    lower_fence, upper_fence, nOut, nType);
fprintf('Outlier types: %s\n', strjoin(cellstr(outTypes), ', '));

%% ===== C. Proximity probability P (outlier types only) =====
% Synapse coordinates (force root_id columns to int64; reconstruct full ids)
optS = detectImportOptions(fullfile(baseDir, 'Codex_Data', 'fafb_v783_princeton_synapse_table.csv'));
optS = setvartype(optS,'pre_root_id_720575940','int64');
optS = setvartype(optS,'post_root_id_720575940','int64');
FAFB_synapse_coordinates = readtable(fullfile(baseDir, 'Codex_Data', 'fafb_v783_princeton_synapse_table.csv'),optS);
FAFB_synapse_coordinates.Properties.VariableNames(11) = "pre_root_id";
FAFB_synapse_coordinates.Properties.VariableNames(12) = "post_root_id";
FAFB_synapse_coordinates.pre_root_id  = FAFB_synapse_coordinates.pre_root_id  + int64(720575940000000000);
FAFB_synapse_coordinates.post_root_id = FAFB_synapse_coordinates.post_root_id + int64(720575940000000000);

% Combined whole right-optic-lobe mesh (Optic_R), from s01
meshDir = fullfile(baseDir, 'Processed_Data', 'optic_lobe_neuropil_mesh');
Optic_R_Faces    = readmatrix(fullfile(meshDir, 'Optic_R_faces.csv')) + 1;   % MATLAB indexing
Optic_R_Vertices = readmatrix(fullfile(meshDir, 'Optic_R_vertices.csv'));

W_thrs = [1 5];
nThr   = numel(W_thrs);
out_P      = nan(nOut, nThr);    % per-neuron pooled (stratified) proximity probability
out_nConn  = nan(nOut, nThr);    % connected-pair count per threshold (directional)
out_nComp  = zeros(nOut, nThr);  % per-neuron comparison counts (pooling weights)
out_nPairs = zeros(nOut, 1);     % valid pair count (centroid distance defined)

for t = 1:nOut
    type_name = char(outTypes(t));
    target_root_ids = root_ids(neuronType == outTypes(t));
    Nn = numel(target_root_ids);
    if Nn < 2, fprintf('[%s] too few neurons, skipped (N=%d)\n', type_name, Nn); continue; end
    fprintf('[%d/%d] proximity %s (N=%d) ...\n', t, nOut, type_name, Nn);

    % 1) Same-type CB connection matrix (directional: row = pre, col = post)
    [~, loc_pre]  = ismember(FAFBConnections.pre_root_id,  target_root_ids);
    [~, loc_post] = ismember(FAFBConnections.post_root_id, target_root_ids);
    is_tt = (loc_pre > 0) & (loc_post > 0);
    tc = FAFBConnections(is_tt, :);
    r_i = loc_pre(is_tt); c_i = loc_post(is_tt); sc = tc.syn_count;
    is_c = ismember(tc.neuropil, FAFBNeuropil_Central);
    rc = [r_i(is_c), c_i(is_c)];
    if isempty(rc)
        W = zeros(Nn, Nn);
    else
        W = accumarray(rc, sc(is_c), [Nn, Nn]);
    end
    W(logical(eye(Nn))) = 0;        % drop autapse; keep direction (no symmetrize)

    % 2) Optic-lobe (Optic_R) centroid distance matrix
    %    Each neuron's OL input synapses (neuron is post) are averaged to a single
    %    centroid; the OL distance between two neurons is the Euclidean distance between centroids.
    cen = nan(Nn,3);
    for i = 1:Nn
        nid = target_root_ids(i);
        in_idx  = (FAFB_synapse_coordinates.post_root_id == nid);   % input synapses only (neuron is post)
        syn_xyz = [FAFB_synapse_coordinates.post_x(in_idx), FAFB_synapse_coordinates.post_y(in_idx), FAFB_synapse_coordinates.post_z(in_idx)];
        if ~isempty(syn_xyz)
            in_roi = intriangulation(Optic_R_Vertices, Optic_R_Faces, syn_xyz);
            syn_R  = syn_xyz(in_roi, :);
            if ~isempty(syn_R)
                cen(i,:) = mean(syn_R, 1);
            end
        end
    end
    Dmat = zeros(Nn);
    for i = 1:Nn
        for j = i+1:Nn
            if any(isnan(cen(i,:))) || any(isnan(cen(j,:))), continue; end
            d = norm(cen(i,:) - cen(j,:));
            Dmat(i,j) = d; Dmat(j,i) = d;
        end
    end

    % 3) Valid pair count + directional connected-pair counts
    mask = triu(true(Nn),1);
    out_nPairs(t) = sum(Dmat(mask) > 0);          % pairs with a defined OL distance
    if out_nPairs(t) < 3, continue; end
    for s = 1:nThr
        out_nConn(t,s) = nnz(W >= W_thrs(s));     % directional connected pairs
    end

    % 4) Per-neuron pooled (stratified) proximity probability P per threshold.
    %    Anchor each neuron i; among its same-type partners (with a defined OL
    %    distance), compare connected (i -> partner, W >= thr) vs unconnected
    %    partner distances. Pool across neurons, weighting each neuron by its
    %    number of comparisons. Ties contribute 0.5 (handled in prob_closer).
    for s = 1:nThr
        sumNum = 0; sumDen = 0;
        for i = 1:Nn
            nb  = [1:i-1, i+1:Nn];
            dij = Dmat(i, nb); wij = W(i, nb);    % i -> partner (i as presynaptic)
            v   = dij > 0;
            dij = dij(v); wij = wij(v);
            conn = wij >= W_thrs(s);
            cd = dij(conn); ud = dij(~conn);
            if isempty(cd) || isempty(ud), continue; end
            den_i  = numel(cd) * numel(ud);
            sumNum = sumNum + prob_closer(cd, ud) * den_i;
            sumDen = sumDen + den_i;
        end
        if sumDen > 0
            out_P(t,s)     = sumNum / sumDen;
            out_nComp(t,s) = sumDen;
        end
    end
end

%% ===== Save outputs =====
D = struct();
D.types         = types;
D.axoaxonic_pct = axoaxonic_pct;   % per-type axo-axonic input from same type (%)
D.n_neurons     = n_neurons;
D.is_outlier    = is_outlier;
D.direction     = direction;
D.q1 = q1; D.q3 = q3; D.iqr = iqr_;
D.lower_fence = lower_fence; D.upper_fence = upper_fence;
D.W_thrs            = W_thrs;
D.out_types         = outTypes;
D.out_axoaxonic_pct = out_axoaxonic;
D.out_P             = out_P;        % per-neuron pooled P (nOut x [W>0, W>=5])
D.out_nConn         = out_nConn;
D.out_nComp         = out_nComp;
D.out_nPairs        = out_nPairs;
save(fullfile(baseDir, 'Processed_Data', 'FFP_axo_axonic_outlier_data.mat'), 'D');

% Human-readable CSVs
PerTypeTable = table(types, axoaxonic_pct, n_neurons, is_outlier, direction, ...
    'VariableNames', {'type','axoaxonic_pct','n_neurons','is_outlier','direction'});
PerTypeTable = sortrows(PerTypeTable, 'axoaxonic_pct', 'descend');
writetable(PerTypeTable, fullfile(baseDir, 'Processed_Data', 'FFP_axo_axonic_perType.csv'));

OutProxTable = table(outTypes, out_axoaxonic, out_P(:,1), out_P(:,2), ...
    out_nConn(:,1), out_nConn(:,2), out_nPairs, ...
    'VariableNames', {'type','axoaxonic_pct','P_W1','P_W5', ...
                      'n_connected_W1','n_connected_W5','n_pairs'});
OutProxTable = sortrows(OutProxTable, 'axoaxonic_pct', 'descend');
writetable(OutProxTable, fullfile(baseDir, 'Processed_Data', 'FFP_axo_axonic_outlier_proximity.csv'));

fprintf('Saved: Processed_Data/FFP_axo_axonic_outlier_data.mat (+ _perType.csv, _outlier_proximity.csv)\n');
fprintf('-> Draw figures with Figures/fig_2F_G_FFP_axo_axonic_plot.m\n');

%% ===== Helper =====
function p = prob_closer(connDist, unconnDist)
% P(connected pair closer than unconnected pair). Ties contribute 0.5.
    connDist = connDist(:); unconnDist = unconnDist(:);
    n1 = numel(connDist); n2 = numel(unconnDist);
    r  = tiedrank([connDist; unconnDist]);
    R1 = sum(r(1:n1));
    p  = (n1*n2 + n1*(n1+1)/2 - R1) / (n1*n2);
end
