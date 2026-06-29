%% s07_FFP_axo_axonic_RF  (compute + CSV save only)
% Optic-lobe (OL) proximity of FFP types as a function of central-brain (CB)
% same-type connection strength W. Same hypothesis as s06, but the readout here
% is the mean RF distance (um) binned by W, rather than a proximity probability.
% RF distance = the OL synapse-centroid distance: each neuron's OL INPUT synapses
% (the neuron is post) are reduced to a single centroid and the pair distance is
% the Euclidean distance between the two centroids. Computed ONLY for the high same-type-input
% (axo-axonic) outlier types, flagged here with the upper Tukey fence (Q3 + 1.5 x IQR).
%
% Hypothesis: the stronger the CB same-type connection (W), the closer the two
% neurons sit in the optic lobe. Directional neuron pairs are binned by W into
% three bins and the mean RF distance (um) is taken per bin.
%
%   bin1  W == 0      : no connection
%   bin2  0 <  W <  5 : weak connection (1-4 synapses)
%   bin3  5 <= W      : strong connection (5+ synapses)
%
% Definitions (same as s06):
%   - W(i,j) = directional CB same-type synapse count, i = pre (anchor)
%   - distance = Euclidean distance (um) between the two neurons' Optic_R
%     input-synapse centroids (input = neuron is post)
%   - per-neuron bin means first, then the per-type value = mean over neurons
%
% Figures are drawn by Figures/fig_2G_H_FFP_axo_axonic_RF_plot.m, which reads the
% saved CSV only. Requires the intriangulation helper (Helper_Function/) on path.
clear all; close all; clc;

baseDir = fileparts(fileparts(mfilename('fullpath')));  % portable; no absolute paths

%% ===== Settings =====
coord_scale = 1e-3;     % FAFB coords -> um (synapse table and mesh share raw units)

binLabels = {'W0','W1to4','W5plus'};                   % CSV column suffixes
binDesc   = {'no connection (W=0)','0<W<5','5<=W'};
nBins     = numel(binLabels);

%% ===== Common data =====
opt = detectImportOptions(fullfile(baseDir, 'Codex_Data', 'connections_no_threshold.csv'));
opt = setvartype(opt,'pre_root_id','int64');
opt = setvartype(opt,'post_root_id','int64');
FAFBConnections = readtable(fullfile(baseDir, 'Codex_Data', 'connections_no_threshold.csv'),opt);

FAFBNeuropils = unique(FAFBConnections.neuropil);
FAFBNeuropil_OpticLobeRight = {'LA_R';'AME_R';'ME_R';'LO_R';'LOP_R'};
FAFBNeuropil_OpticLobeLeft  = {'LA_L';'AME_L';'ME_L';'LO_L';'LOP_L'};
FAFBNeuropil_Central = FAFBNeuropils(~ismember(FAFBNeuropils,FAFBNeuropil_OpticLobeRight));
FAFBNeuropil_Central = FAFBNeuropil_Central(~ismember(FAFBNeuropil_Central,FAFBNeuropil_OpticLobeLeft));
FAFBNeuropil_Central = FAFBNeuropil_Central(~ismember(FAFBNeuropil_Central,'UNASGD'));

% Synapse coordinates (force root_id columns to int64; reconstruct full ids)
optS = detectImportOptions(fullfile(baseDir, 'Codex_Data', 'fafb_v783_princeton_synapse_table.csv'));
optS = setvartype(optS,'pre_root_id_720575940','int64');
optS = setvartype(optS,'post_root_id_720575940','int64');
FAFB_synapse_coordinates = readtable(fullfile(baseDir, 'Codex_Data', 'fafb_v783_princeton_synapse_table.csv'),optS);
FAFB_synapse_coordinates.Properties.VariableNames(11) = "pre_root_id";
FAFB_synapse_coordinates.Properties.VariableNames(12) = "post_root_id";
FAFB_synapse_coordinates.pre_root_id  = FAFB_synapse_coordinates.pre_root_id  + int64(720575940000000000);
FAFB_synapse_coordinates.post_root_id = FAFB_synapse_coordinates.post_root_id + int64(720575940000000000);

% ----- coords to um (distance unit = um) -----
FAFB_synapse_coordinates.pre_x  = FAFB_synapse_coordinates.pre_x  * coord_scale;
FAFB_synapse_coordinates.pre_y  = FAFB_synapse_coordinates.pre_y  * coord_scale;
FAFB_synapse_coordinates.pre_z  = FAFB_synapse_coordinates.pre_z  * coord_scale;
FAFB_synapse_coordinates.post_x = FAFB_synapse_coordinates.post_x * coord_scale;
FAFB_synapse_coordinates.post_y = FAFB_synapse_coordinates.post_y * coord_scale;
FAFB_synapse_coordinates.post_z = FAFB_synapse_coordinates.post_z * coord_scale;

% Combined whole right-optic-lobe mesh (Optic_R), from s01
meshDir = fullfile(baseDir, 'Processed_Data', 'optic_lobe_neuropil_mesh');
Optic_R_Faces    = readmatrix(fullfile(meshDir, 'Optic_R_faces.csv')) + 1;   % MATLAB indexing
Optic_R_Vertices = readmatrix(fullfile(meshDir, 'Optic_R_vertices.csv'));
Optic_R_Vertices = Optic_R_Vertices * coord_scale;   % match the synapse-coord unit (um)

%% ===== FFP neuron table =====
load(fullfile(baseDir, 'Processed_Data', 'right_neurons_thr0.mat'), 'RightFFP_NPIs')
RightFF_NPIs = RightFFP_NPIs;    % alias (this repo names the table RightFFP_NPIs)

% Types only used for verbose console output (independent of compute scope)
want_types = {'MeTu2a'};

%% ===== High same-type-input (axo-axonic) outlier types =====
% central = neuropil that is neither optic lobe (L/R) nor UNASGD (same definition
% as FAFBNeuropil_Central used for the proximity step).
FAFBNeuropil_OpticLobe_All = [FAFBNeuropil_OpticLobeRight; FAFBNeuropil_OpticLobeLeft];

oc_root_ids = RightFF_NPIs.root_id;                        % N x 1
oc_N        = numel(oc_root_ids);
[oc_typeCats, ~, oc_typeCode] = unique(RightFF_NPIs.type); % per-neuron type code
oc_nTypes   = numel(oc_typeCats);

[~, oc_preLoc]  = ismember(FAFBConnections.pre_root_id,  oc_root_ids);
[~, oc_postLoc] = ismember(FAFBConnections.post_root_id, oc_root_ids);
oc_syn = FAFBConnections.syn_count;

oc_isPostTarget = oc_postLoc > 0;                          % post is a target neuron
oc_isOptic      = ismember(FAFBConnections.neuropil, FAFBNeuropil_OpticLobe_All);
oc_isUnasgd     = strcmp(FAFBConnections.neuropil, 'UNASGD');
oc_isCentral    = ~oc_isOptic & ~oc_isUnasgd;             % neither optic nor UNASGD -> central

oc_isBoth = (oc_preLoc > 0) & (oc_postLoc > 0);          % self connection = pre/post both targets and same type
oc_isSelf = false(height(FAFBConnections),1);
oc_isSelf(oc_isBoth) = oc_typeCode(oc_preLoc(oc_isBoth)) == oc_typeCode(oc_postLoc(oc_isBoth));

oc_acc = @(m) accumarray(oc_postLoc(m), oc_syn(m), [oc_N 1]);   % synapse sum per post neuron
oc_self_optic    = oc_acc(oc_isPostTarget & oc_isOptic   & oc_isSelf);
oc_self_central  = oc_acc(oc_isPostTarget & oc_isCentral & oc_isSelf);
oc_other_optic   = oc_acc(oc_isPostTarget & oc_isOptic   & ~oc_isSelf);
oc_other_central = oc_acc(oc_isPostTarget & oc_isCentral & ~oc_isSelf);
oc_total_in = oc_self_optic + oc_self_central + oc_other_optic + oc_other_central;

% Per-neuron self central input fraction (%) -- exclude neurons with no input
oc_valid = oc_total_in > 0;
oc_selfCentralPct = nan(oc_N,1);
oc_selfCentralPct(oc_valid) = oc_self_central(oc_valid) ./ oc_total_in(oc_valid) * 100;

% Per-type mean self central input (%)
oc_mean_pct  = nan(oc_nTypes,1);
oc_n_neurons = zeros(oc_nTypes,1);
for k = 1:oc_nTypes
    rows = (oc_typeCode == k & oc_valid);
    oc_n_neurons(k) = sum(rows);
    if oc_n_neurons(k) > 0
        oc_mean_pct(k) = mean(oc_selfCentralPct(rows));
    end
end

% Outlier types : Tukey rule Q3+1.5*IQR (upper) / Q1-1.5*IQR (lower)
oc_typesWithData = oc_typeCats(oc_n_neurons > 0);
oc_nWithData     = oc_n_neurons(oc_n_neurons > 0);
oc_allVals       = oc_mean_pct(oc_n_neurons > 0);
oc_q1   = prctile(oc_allVals, 25);
oc_q3   = prctile(oc_allVals, 75);
oc_iqr  = oc_q3 - oc_q1;
oc_whisker    = 1.5;
oc_upperFence = oc_q3 + oc_whisker*oc_iqr;
oc_lowerFence = oc_q1 - oc_whisker*oc_iqr;

oc_isOutlier = (oc_allVals > oc_upperFence) | (oc_allVals < oc_lowerFence);
oc_outDir    = repmat({'high'}, numel(oc_allVals), 1);    % above / below
oc_outDir(oc_allVals < oc_lowerFence) = {'low'};

OutlierTable = table(oc_typesWithData(oc_isOutlier), oc_allVals(oc_isOutlier), oc_outDir(oc_isOutlier), oc_nWithData(oc_isOutlier), ...
    'VariableNames', {'type','self_central_pct','direction','n_neurons'});
OutlierTable = sortrows(OutlierTable, 'self_central_pct', 'descend');
writetable(OutlierTable, fullfile(baseDir, 'Processed_Data', 'FFP_axo_axonic_RF_outlier_types.csv'));
fprintf('Outlier fence: Q1=%.2f, Q3=%.2f, IQR=%.2f -> fence [%.2f, %.2f]%%\n', oc_q1, oc_q3, oc_iqr, oc_lowerFence, oc_upperFence);
fprintf('%d outlier types -> FFP_axo_axonic_RF_outlier_types.csv saved.\n', height(OutlierTable));

%% ===== Analysis scope : high (upper-fence) outlier types only =====
analysis_types = OutlierTable.type(strcmp(OutlierTable.direction, 'high'));
if isempty(analysis_types)
    error('No high (upper-fence) outlier type found. Check the fence criterion.');
end
fprintf('OL-distance scope (high outliers) %d types: %s\n\n', ...
    numel(analysis_types), strjoin(cellstr(analysis_types), ', '));

nT = numel(analysis_types);
% Per-bin (nT x nBins) metrics
Dbin_type   = nan(nT,nBins);    % per-type bin mean distance (um)
SEMbin_type = nan(nT,nBins);    % within-type SEM over neurons
nNeuron_bin = zeros(nT,nBins);  % anchor neurons with >=1 partner in the bin
nPair_bin   = zeros(nT,nBins);  % directional neuron pairs in the bin

% connected(W>0) vs unconnected(W==0) split (for the bar graph)
Dconn_type     = nan(nT,1);   SEMconn_type   = nan(nT,1);   nConn_neuron   = zeros(nT,1);   nConn_pair   = zeros(nT,1);
Dunconn_type   = nan(nT,1);   SEMunconn_type = nan(nT,1);   nUnconn_neuron = zeros(nT,1);   nUnconn_pair = zeros(nT,1);

for t = 1:nT
    type_name = char(analysis_types(t));
    target_root_ids = RightFF_NPIs.root_id(strcmp(RightFF_NPIs.type, type_name), :);
    N = numel(target_root_ids);
    if N < 2, fprintf('[%s] too few neurons, skipped (N=%d)\n', type_name, N); continue; end
    fprintf('[%d/%d] %s (N=%d) ...\n', t, nT, type_name, N);

    %% 1) CB same-type connection matrix (directional: row = pre, i = anchor)
    [~, loc_pre]  = ismember(FAFBConnections.pre_root_id,  target_root_ids);
    [~, loc_post] = ismember(FAFBConnections.post_root_id, target_root_ids);
    is_tt = (loc_pre > 0) & (loc_post > 0);
    tc = FAFBConnections(is_tt, :);
    r_i = loc_pre(is_tt); c_i = loc_post(is_tt); sc = tc.syn_count;
    is_c = ismember(tc.neuropil, FAFBNeuropil_Central);
    rc = [r_i(is_c), c_i(is_c)];
    if isempty(rc)
        W = zeros(N, N);
    else
        W = accumarray(rc, sc(is_c), [N, N]);
    end

    %% 2) Optic-lobe (Optic_R) synapse-centroid distance matrix (um)
    cen = nan(N,3);
    for i = 1:N
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
    Dmat = zeros(N);
    for i = 1:N
        for j = i+1:N
            if any(isnan(cen(i,:))) || any(isnan(cen(j,:))), continue; end
            d = norm(cen(i,:) - cen(j,:));   % um
            Dmat(i,j) = d; Dmat(j,i) = d;
        end
    end

    %% 3) Per-neuron bin means -> per-type value (bins + connected/unconnected split)
    dbin_i = nan(N, nBins);     % anchor i bin mean distance
    npair_acc = zeros(1,nBins); % type-wide per-bin pair count
    dconn_i   = nan(N,1);       % anchor i connected(W>0) mean distance
    dunconn_i = nan(N,1);       % anchor i unconnected(W==0) mean distance
    nconn_acc = 0; nunconn_acc = 0;
    for i = 1:N
        nb  = [1:i-1, i+1:N];
        dij = Dmat(i, nb); wij = W(i, nb);    % directional: i sending
        v   = dij > 0;                        % pairs with a defined centroid distance
        dij = dij(v); wij = wij(v);
        if isempty(dij), continue; end

        bmask = binAssign(wij, nBins);        % nBins x numel(wij) logical (each row = one bin)
        for b = 1:nBins
            sel = bmask(b,:);
            if any(sel)
                dbin_i(i,b)  = mean(dij(sel));
                npair_acc(b) = npair_acc(b) + sum(sel);
            end
        end

        % connected(W>0) / unconnected(W==0) split
        cmask = wij > 0;  umask = wij == 0;
        if any(cmask), dconn_i(i)   = mean(dij(cmask));  nconn_acc   = nconn_acc   + sum(cmask); end
        if any(umask), dunconn_i(i) = mean(dij(umask));  nunconn_acc = nunconn_acc + sum(umask); end
    end

    for b = 1:nBins
        ok = ~isnan(dbin_i(:,b));
        nNeuron_bin(t,b) = sum(ok);
        nPair_bin(t,b)   = npair_acc(b);
        if sum(ok) == 0, continue; end
        d = dbin_i(ok,b);
        Dbin_type(t,b)   = mean(d);
        SEMbin_type(t,b) = std(d) / sqrt(numel(d));
    end

    % connected/unconnected per-type value
    okc = ~isnan(dconn_i);   nConn_neuron(t)   = sum(okc); nConn_pair(t)   = nconn_acc;
    if any(okc), Dconn_type(t)   = mean(dconn_i(okc));   SEMconn_type(t)   = std(dconn_i(okc))   / sqrt(sum(okc)); end
    oku = ~isnan(dunconn_i); nUnconn_neuron(t) = sum(oku); nUnconn_pair(t) = nunconn_acc;
    if any(oku), Dunconn_type(t) = mean(dunconn_i(oku)); SEMunconn_type(t) = std(dunconn_i(oku)) / sqrt(sum(oku)); end
end

%% ===== Overall summary (bin trend) =====
fprintf('\n========== Quantification (high outlier types, per-bin mean distance, um) ==========\n');
for b = 1:nBins
    valid = ~isnan(Dbin_type(:,b));
    fprintf('  [%-18s] valid types %3d | mean distance (type mean)=%.3f um | total pairs=%d\n', ...
        binDesc{b}, sum(valid), mean(Dbin_type(valid,b)), sum(nPair_bin(:,b)));
end

% Trend test over types defined in all bins (Friedman) + W0 pairwise
allbin = all(~isnan(Dbin_type), 2);
nComplete = sum(allbin);
fprintf('  Types defined in all %d bins: %d\n', nBins, nComplete);
if nComplete >= 3
    M = Dbin_type(allbin, :);
    try
        p_fried = friedman(M, 1, 'off');
        fprintf('  Friedman p(bin difference) = %.3g\n', p_fried);
    catch ME
        fprintf('  Friedman test failed: %s\n', ME.message);
    end
    for b = 2:nBins
        % closer than W0? (distance decrease = supports hypothesis) : tail=left
        p = signrank(M(:,b), M(:,1), 'tail','left');
        fprintf('  signrank %s < W0 (closer) p=%.3g | median diff = %.3f um\n', ...
            binLabels{b}, p, median(M(:,b)-M(:,1)));
    end
end
fprintf('================================================================\n');

%% ===== CSV save =====
vars = {}; names = {};
vars{end+1} = analysis_types(:);       names{end+1} = 'type';
for b = 1:nBins
    vars{end+1} = Dbin_type(:,b);      names{end+1} = ['Distance_um_'  binLabels{b}];
    vars{end+1} = SEMbin_type(:,b);    names{end+1} = ['SEM_um_'       binLabels{b}];
    vars{end+1} = nNeuron_bin(:,b);    names{end+1} = ['n_neurons_'    binLabels{b}];
    vars{end+1} = nPair_bin(:,b);      names{end+1} = ['n_pairs_'      binLabels{b}];
end
% connected(W>0) / unconnected(W==0) split columns (for the bar graph)
vars{end+1} = Dconn_type;     names{end+1} = 'Distance_um_connected';
vars{end+1} = SEMconn_type;   names{end+1} = 'SEM_um_connected';
vars{end+1} = nConn_neuron;   names{end+1} = 'n_neurons_connected';
vars{end+1} = nConn_pair;     names{end+1} = 'n_pairs_connected';
vars{end+1} = Dunconn_type;   names{end+1} = 'Distance_um_unconnected';
vars{end+1} = SEMunconn_type; names{end+1} = 'SEM_um_unconnected';
vars{end+1} = nUnconn_neuron; names{end+1} = 'n_neurons_unconnected';
vars{end+1} = nUnconn_pair;   names{end+1} = 'n_pairs_unconnected';
ResTable = table(vars{:}, 'VariableNames', names);
writetable(ResTable, fullfile(baseDir, 'Processed_Data', 'FFP_axo_axonic_RF_proximity.csv'));
fprintf('Saved: Processed_Data/FFP_axo_axonic_RF_proximity.csv (high outlier types %d, %d bins: %s)\n', ...
    nT, nBins, strjoin(binLabels, ', '));
fprintf('-> Draw figures with Figures/fig_2G_H_FFP_axo_axonic_RF_plot.m\n');

%% ===== Verbose console output for want_types =====
detail_idx = find(ismember(analysis_types, want_types));
for k = 1:numel(detail_idx)
    t = detail_idx(k);
    fprintf('[%s] ', char(analysis_types(t)));
    for b = 1:nBins
        fprintf('%s=%.3f um(n=%d) ', binLabels{b}, Dbin_type(t,b), nNeuron_bin(t,b));
    end
    fprintf('\n');
end

%% ===== Helper : connection-strength bin assignment =====
% wij (1xM) -> nBins x M logical matrix (each row = one bin)
function bmask = binAssign(wij, nBins)
    wij = wij(:)';
    bmask = false(nBins, numel(wij));
    bmask(1,:) = (wij == 0);                 % W0  : no connection
    bmask(2,:) = (wij >= 1) & (wij < 5);     % 0<W<5
    bmask(3,:) = (wij >= 5);                 % 5<=W
end
