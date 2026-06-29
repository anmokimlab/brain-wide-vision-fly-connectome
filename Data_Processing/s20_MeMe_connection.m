%% s20_MeMe_connection
% Build the neuron-level connectivity used by the MeMe network simulation
% (Figure S5D-H). This is the data-only step: it produces the connection .mat
% that Figures\fig_S5D_E_F_G_H_MeMe_simulation.m loads. No figure is drawn here
% (the connection-graph figure is the L/R-split one in fig_S5A_B_*).
%
% For the simulation cell types (MeMe_e01/e02, Sm07, Dm2, MeTu1, MeTu3c) it
% reads connections_no_threshold.csv and consolidated_cell_types.csv from
% Codex_Data\, removes connectome-error minority-NT cells and sparse type-pair
% edges, and writes to Processed_Data\MeMe_AllSimNeurons_Graph.mat:
%   pair_table   - per (pre,post) neuron pair (signed_syn = sign(NT) * syn_count)
%   neuron_table - root_id / primary_type for every simulation neuron
%   edge_table   - type-pair edge table (G.Edges; mean syn per connected pair, signed)
%   sim_types, center_types, inhib_nts, weight_mode
%
% Only the 'absolute' edge-weight definition is built (the recv-mean variant of
% the original analysis script is not needed by the simulation).
clear all; close all; clc

baseDir = fileparts(fileparts(mfilename('fullpath')));  % resolved from this script's location (portable; no absolute paths)
codexDir = fullfile(baseDir, 'Codex_Data');
outDir   = fullfile(baseDir, 'Processed_Data');
if ~exist(outDir, 'dir'), mkdir(outDir); end

%% --- Settings ---
sim_types = {'MeMe_e01','MeMe_e02','Sm07','Dm2','MeTu1','MeTu3c'};
center_types = {'MeMe_e01','MeMe_e02'};
exclude_neuropils  = {'UNASGD'};
inhib_nts    = {'GABA','GLUT'};
edge_threshold = 0;
drop_isolated  = false;
weight_mode = 'absolute';        % mean over connected (pre,post) pairs of syn(a->b)
nt_minority_warn_frac = 0.10;
edge_min_frac_of_tgt_input = 0.01;

%% --- Load data ---
opt = detectImportOptions(fullfile(codexDir, 'connections_no_threshold.csv'));
opt = setvartype(opt,'pre_root_id','int64');
opt = setvartype(opt,'post_root_id','int64');
FAFBConnections = readtable(fullfile(codexDir, 'connections_no_threshold.csv'), opt);
FAFBConnections(FAFBConnections.syn_count < 5, :) = [];

opt = detectImportOptions(fullfile(codexDir, 'consolidated_cell_types.csv'));
opt = setvartype(opt,'root_id','int64');
FAFBConsolidated_type = readtable(fullfile(codexDir, 'consolidated_cell_types.csv'), opt);

% Remove UNASGD (unassigned neuropil)
FAFBConnections = FAFBConnections(~ismember(FAFBConnections.neuropil, exclude_neuropils), :);

%% --- Node definition (use sim_types as-is) ---
sim_types = sim_types(:);
miss = sim_types(~ismember(sim_types, FAFBConsolidated_type.primary_type));
if ~isempty(miss)
    warning('Types in sim_types not present in consolidated_cell_types: %s', strjoin(miss, ', '));
end
all_nodes = unique(sim_types, 'stable');
fprintf('Total nodes (sim types) : %d\n', numel(all_nodes));

%% --- type -> root ids ---
type_to_root_ids = containers.Map();
for i = 1:numel(all_nodes)
    t = all_nodes{i};
    type_to_root_ids(t) = unique( ...
        FAFBConsolidated_type.root_id(strcmp(FAFBConsolidated_type.primary_type, t)));
    fprintf('   %-10s : %d cells\n', t, numel(type_to_root_ids(t)));
end

%% --- pre_root_id -> nt_type mapping ---
[u_pre_all, ia_pre] = unique(FAFBConnections.pre_root_id);
nt_per_pre = FAFBConnections.nt_type(ia_pre);
pre2nt = containers.Map(num2cell(u_pre_all), nt_per_pre);

%% --- Within-type NT consistency check + remove minority cells ---
fprintf('\n--- NT consistency check (minority-NT cells -> removed) ---\n');
n_removed_total_nt = 0;
for i = 1:numel(all_nodes)
    tp = all_nodes{i};
    type_root_ids = type_to_root_ids(tp);
    if isempty(type_root_ids), continue; end
    nt_per_cell = repmat({''}, numel(type_root_ids), 1);
    for k = 1:numel(type_root_ids)
        if isKey(pre2nt, type_root_ids(k)), nt_per_cell{k} = pre2nt(type_root_ids(k)); end
    end
    has_nt = ~cellfun(@isempty, nt_per_cell);
    if ~any(has_nt), continue; end
    nt_known = nt_per_cell(has_nt);
    uNT_t    = unique(nt_known);
    counts   = zeros(numel(uNT_t),1);
    for u = 1:numel(uNT_t), counts(u) = sum(strcmp(nt_known, uNT_t{u})); end
    [~, max_idx] = max(counts);
    dom_NT_t = uNT_t{max_idx};
    n_tot    = sum(has_nt);
    minority_mask = has_nt & ~strcmp(nt_per_cell, dom_NT_t);
    n_min = sum(minority_mask);
    if n_min == 0, continue; end
    if n_min / n_tot >= nt_minority_warn_frac
        fprintf('   [WARN] %s : NT heterogeneous (minority=%.1f%%) -> KEEP ALL\n', tp, n_min/n_tot*100);
        continue;
    end
    type_to_root_ids(tp) = type_root_ids(~minority_mask);
    n_removed_total_nt = n_removed_total_nt + n_min;
end
fprintf('Total minority-NT cells removed across all types: %d\n\n', n_removed_total_nt);

%% --- Build edges ---
src_list = {}; tgt_list = {};
abs_w = []; sgn = []; dom_nt = {}; n_pair = [];
mean_n_src_per_tgt = []; mean_n_tgt_per_src = [];
n_src_cells_all = []; n_tgt_cells_all = [];
n_src_connected = []; n_tgt_connected = [];
total_syn_all   = []; mean_syn_per_pair = [];

pair_pre = zeros(0,1,'int64');  pair_post = zeros(0,1,'int64');
pair_syn = zeros(0,1);  pair_src_t = {};  pair_tgt_t = {};
pair_nt  = {};  pair_sgn = zeros(0,1);

for i = 1:numel(all_nodes)
    for j = 1:numel(all_nodes)
        if i == j, continue; end
        src = all_nodes{i};  tgt = all_nodes{j};
        source_root_ids = type_to_root_ids(src);  target_root_ids = type_to_root_ids(tgt);
        if isempty(source_root_ids) || isempty(target_root_ids), continue; end
        rows = FAFBConnections( ...
            ismember(FAFBConnections.pre_root_id,  source_root_ids) & ...
            ismember(FAFBConnections.post_root_id, target_root_ids), :);
        if isempty(rows), continue; end

        [uPair, ~, ic] = unique([rows.pre_root_id rows.post_root_id], 'rows');
        pairSyn = accumarray(ic, rows.syn_count, [], @sum);
        n_existing_pairs = numel(pairSyn);
        total_syn    = sum(rows.syn_count);
        n_post_cells = numel(target_root_ids);
        n_pre_cells  = numel(source_root_ids);

        w = mean(pairSyn);                      % absolute weight
        mean_src_per_tgt = n_existing_pairs / n_post_cells;
        mean_tgt_per_src = n_existing_pairs / n_pre_cells;
        n_src_conn = numel(unique(rows.pre_root_id));
        n_tgt_conn = numel(unique(rows.post_root_id));

        nts = rows.nt_type;
        [uNT, ~, icnt] = unique(nts);
        ntTotals = accumarray(icnt, rows.syn_count, [], @sum);
        [~, idxNT] = max(ntTotals);
        domNT = uNT{idxNT};
        s = 1;  if ismember(domNT, inhib_nts), s = -1; end

        src_list{end+1,1} = src;  tgt_list{end+1,1} = tgt; %#ok<SAGROW>
        abs_w(end+1,1) = w;  sgn(end+1,1) = s; %#ok<SAGROW>
        dom_nt{end+1,1} = domNT;  n_pair(end+1,1) = n_existing_pairs; %#ok<SAGROW>
        mean_n_src_per_tgt(end+1,1) = mean_src_per_tgt; %#ok<SAGROW>
        mean_n_tgt_per_src(end+1,1) = mean_tgt_per_src; %#ok<SAGROW>
        n_src_cells_all(end+1,1) = n_pre_cells;  n_tgt_cells_all(end+1,1) = n_post_cells; %#ok<SAGROW>
        n_src_connected(end+1,1) = n_src_conn;   n_tgt_connected(end+1,1) = n_tgt_conn; %#ok<SAGROW>
        total_syn_all(end+1,1) = total_syn;  mean_syn_per_pair(end+1,1) = mean(pairSyn); %#ok<SAGROW>

        pre_root_ids_pair  = uPair(:,1);  post_root_ids_pair = uPair(:,2);
        nt_p = cell(n_existing_pairs, 1);
        for kk = 1:n_existing_pairs
            if isKey(pre2nt, pre_root_ids_pair(kk)), nt_p{kk} = pre2nt(pre_root_ids_pair(kk));
            else, nt_p{kk} = ''; end
        end
        sgn_p = ones(n_existing_pairs, 1);
        sgn_p(ismember(nt_p, inhib_nts)) = -1;

        pair_pre   = [pair_pre;   int64(pre_root_ids_pair)];
        pair_post  = [pair_post;  int64(post_root_ids_pair)];
        pair_syn   = [pair_syn;   pairSyn];
        pair_src_t = [pair_src_t; repmat({src}, n_existing_pairs, 1)];
        pair_tgt_t = [pair_tgt_t; repmat({tgt}, n_existing_pairs, 1)];
        pair_nt    = [pair_nt;    nt_p];
        pair_sgn   = [pair_sgn;   sgn_p];
    end
end

%% --- Remove sparse type-pair edges ---
total_input_per_type  = containers.Map();
total_output_per_type = containers.Map();
for i = 1:numel(all_nodes)
    tp = all_nodes{i};
    type_root_ids = type_to_root_ids(tp);
    if isempty(type_root_ids)
        total_input_per_type(tp)  = 0;  total_output_per_type(tp) = 0;
    else
        total_input_per_type(tp)  = sum(FAFBConnections.syn_count(ismember(FAFBConnections.post_root_id, type_root_ids)));
        total_output_per_type(tp) = sum(FAFBConnections.syn_count(ismember(FAFBConnections.pre_root_id,  type_root_ids)));
    end
end

keep_edge_by_frac  = true(numel(src_list), 1);
removed_type_pairs = {};
for k = 1:numel(src_list)
    total_out_src = total_output_per_type(src_list{k});
    total_in_tgt  = total_input_per_type(tgt_list{k});
    if total_out_src > 0, frac_out = total_syn_all(k) / total_out_src; else, frac_out = 0; end
    if total_in_tgt  > 0, frac_in  = total_syn_all(k) / total_in_tgt;  else, frac_in  = 0; end
    if max(frac_in, frac_out) < edge_min_frac_of_tgt_input
        keep_edge_by_frac(k) = false;
        removed_type_pairs{end+1,1} = [src_list{k} '__' tgt_list{k}]; %#ok<SAGROW>
    end
end
fprintf('Edges removed by fraction filter: %d / %d\n', sum(~keep_edge_by_frac), numel(src_list));

src_list = src_list(keep_edge_by_frac);  tgt_list = tgt_list(keep_edge_by_frac);
abs_w = abs_w(keep_edge_by_frac);  sgn = sgn(keep_edge_by_frac);
dom_nt = dom_nt(keep_edge_by_frac);  n_pair = n_pair(keep_edge_by_frac);
mean_n_src_per_tgt = mean_n_src_per_tgt(keep_edge_by_frac);
mean_n_tgt_per_src = mean_n_tgt_per_src(keep_edge_by_frac);
n_src_cells_all = n_src_cells_all(keep_edge_by_frac);
n_tgt_cells_all = n_tgt_cells_all(keep_edge_by_frac);
n_src_connected = n_src_connected(keep_edge_by_frac);
n_tgt_connected = n_tgt_connected(keep_edge_by_frac);
total_syn_all = total_syn_all(keep_edge_by_frac);
mean_syn_per_pair = mean_syn_per_pair(keep_edge_by_frac);

if ~isempty(removed_type_pairs)
    keep_pair_row = true(numel(pair_pre), 1);
    for r = 1:numel(removed_type_pairs)
        parts = strsplit(removed_type_pairs{r}, '__');
        keep_pair_row = keep_pair_row & ~(strcmp(pair_src_t, parts{1}) & strcmp(pair_tgt_t, parts{2}));
    end
    pair_pre = pair_pre(keep_pair_row);  pair_post = pair_post(keep_pair_row);
    pair_syn = pair_syn(keep_pair_row);  pair_src_t = pair_src_t(keep_pair_row);
    pair_tgt_t = pair_tgt_t(keep_pair_row);  pair_nt = pair_nt(keep_pair_row);
    pair_sgn = pair_sgn(keep_pair_row);
end

signed_w = abs_w .* sgn;

%% --- threshold ---
keep_mask = abs_w >= edge_threshold;
src_list = src_list(keep_mask);  tgt_list = tgt_list(keep_mask);
abs_w = abs_w(keep_mask);  sgn = sgn(keep_mask);
dom_nt = dom_nt(keep_mask);  n_pair = n_pair(keep_mask);
mean_n_src_per_tgt = mean_n_src_per_tgt(keep_mask);
mean_n_tgt_per_src = mean_n_tgt_per_src(keep_mask);
n_src_cells_all = n_src_cells_all(keep_mask);
n_tgt_cells_all = n_tgt_cells_all(keep_mask);
n_src_connected = n_src_connected(keep_mask);
n_tgt_connected = n_tgt_connected(keep_mask);
total_syn_all = total_syn_all(keep_mask);
mean_syn_per_pair = mean_syn_per_pair(keep_mask);
signed_w = signed_w(keep_mask);
fprintf('Edges kept : %d\n', numel(signed_w));

%% --- digraph (absolute weight) -> edge_table ---
G = digraph(src_list, tgt_list, signed_w, all_nodes);
edge_idx = findedge(G, src_list, tgt_list);
nE = height(G.Edges);
AbsW_vec = nan(nE,1);  DomNT_vec = repmat({''}, nE,1);  NPair_vec = nan(nE,1);
MeanSrcPerTgt = nan(nE,1);  MeanTgtPerSrc = nan(nE,1);
NSrcCellsAll = nan(nE,1);  NTgtCellsAll = nan(nE,1);
NSrcConnected = nan(nE,1); NTgtConnected = nan(nE,1);
TotalSyn_vec = nan(nE,1);  MeanSynPerPair = nan(nE,1);
ok = edge_idx > 0;
AbsW_vec(edge_idx(ok)) = abs_w(ok);  DomNT_vec(edge_idx(ok)) = dom_nt(ok);
NPair_vec(edge_idx(ok)) = n_pair(ok);
MeanSrcPerTgt(edge_idx(ok)) = mean_n_src_per_tgt(ok);
MeanTgtPerSrc(edge_idx(ok)) = mean_n_tgt_per_src(ok);
NSrcCellsAll(edge_idx(ok)) = n_src_cells_all(ok);
NTgtCellsAll(edge_idx(ok)) = n_tgt_cells_all(ok);
NSrcConnected(edge_idx(ok)) = n_src_connected(ok);
NTgtConnected(edge_idx(ok)) = n_tgt_connected(ok);
TotalSyn_vec(edge_idx(ok)) = total_syn_all(ok);
MeanSynPerPair(edge_idx(ok)) = mean_syn_per_pair(ok);
G.Edges.AbsWeight = AbsW_vec;  G.Edges.DomNT = DomNT_vec;  G.Edges.NPair = NPair_vec;
G.Edges.MeanSrcPerTgt = MeanSrcPerTgt;  G.Edges.MeanTgtPerSrc = MeanTgtPerSrc;
G.Edges.NSrcCellsAll = NSrcCellsAll;  G.Edges.NTgtCellsAll = NTgtCellsAll;
G.Edges.NSrcConnected = NSrcConnected;  G.Edges.NTgtConnected = NTgtConnected;
G.Edges.TotalSyn = TotalSyn_vec;  G.Edges.MeanSynPerPair = MeanSynPerPair;

if drop_isolated
    deg_total = indegree(G) + outdegree(G);
    iso_mask  = (deg_total == 0) & ~ismember(G.Nodes.Name, center_types);
    if any(iso_mask), G = rmnode(G, G.Nodes.Name(iso_mask)); end
end
edge_table = G.Edges;

%% --- neuron-level pair table ---
pair_table = table( ...
    pair_pre, pair_post, pair_src_t, pair_tgt_t, ...
    pair_syn, pair_nt, pair_sgn, pair_sgn .* pair_syn, ...
    'VariableNames', {'pre_root_id','post_root_id','src_type','tgt_type', ...
                      'syn_count','nt_type','sign','signed_syn'});
fprintf('Pair table   : %d connected (pre,post) pairs\n', height(pair_table));

%% --- neuron-level neuron table ---
neuron_rid  = zeros(0,1,'int64');  neuron_type = {};
for ii = 1:numel(all_nodes)
    type_root_ids = type_to_root_ids(all_nodes{ii});
    neuron_rid  = [neuron_rid;  int64(type_root_ids)];
    neuron_type = [neuron_type; repmat(all_nodes(ii), numel(type_root_ids), 1)];
end
neuron_table = table(neuron_rid, neuron_type, 'VariableNames', {'root_id','primary_type'});
fprintf('Neuron table : %d neurons across %d sim types\n', height(neuron_table), numel(all_nodes));

%% --- save ---
save_path = fullfile(outDir, 'MeMe_AllSimNeurons_Graph.mat');
save(save_path, 'G', 'edge_table', 'all_nodes', 'sim_types', 'center_types', ...
    'weight_mode', 'inhib_nts', 'pair_table', 'neuron_table');
fprintf('Saved : %s\n', save_path);
