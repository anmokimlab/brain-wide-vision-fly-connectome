%% fig_S5A_B_MeMe_connection_graph
%  MeMe / Sm / Dm2 / MeTu simulation cell types, left/right (side) split into
%  separate nodes (e.g. MeMe_e01 -> MeMe_e01_L / MeMe_e01_R). Side comes from
%  the 'side' column of classification.csv.
%
%  It builds the type-pair connectivity (edge weight = mean synapse count over
%  connected (pre,post) neuron pairs; sign from the dominant neurotransmitter,
%  GABA/GLUT = inhibitory = negative) and produces:
%
%   - Figure 1 (paper panel S5A) - the L/R-split connection digraph
%     (orange = excitatory, blue = inhibitory; edge width prop. to |weight|).
%   - Figure 2 (paper panel S5B) - the same edge weights as a signed
%     source(pre) x target(post) connection matrix.
%
%  Reads the Codex tables (connections_no_threshold.csv, consolidated_cell_types.csv,
%  classification.csv) from Codex_Data\. Only the 'absolute' weight definition is
%  drawn (the recv-mean variant of the original analysis script is not a paper panel).
clear all; close all; clc

baseDir = fileparts(fileparts(mfilename('fullpath')));  % resolved from this script's location (portable; no absolute paths)
codexDir = fullfile(baseDir, 'Codex_Data');

%% --- Settings ---
% All neuron types used in the simulation (base type; each split into _L/_R)
sim_types = {'MeMe_e01','MeMe_e02','Sm07','Dm2','MeTu1','MeTu3c'};

% Center types to highlight in the visualization (by base type; both _L/_R highlighted)
center_types = {'MeMe_e01','MeMe_e02'};

exclude_neuropils  = {'UNASGD'};
inhib_nts    = {'GABA','GLUT'};

edge_threshold = 0;          % remove edges with |weight| below this value
drop_isolated  = false;      % keep false to show the used types as-is

% --- Side-split settings ---
%   Warn if the fraction of neurons without a side (center/empty) is at least this value.
side_missing_warn_frac = 0.05;   % 5%

% --- Edge-weight computation mode ---
%   'absolute' : mean_{(a,b) connected} [ syn(a->b) ]
%   'fraction' : mean_{(a,b) connected} [ syn(a->b) / total_input_syn(b) ]
weight_mode = 'absolute';

% --- Within-type minority-NT handling (remove connectome errors) ---
nt_minority_warn_frac = 0.10;

% --- Remove sparse type-pair edges ---
edge_min_frac_of_tgt_input = 0.01;

%% =========================================================
% Layout settings
%   'custom_lr' : L on the left / R on the right, same base type at the same height
%   'force' | 'layered' | 'circle' ... : built-in MATLAB layout
%% =========================================================
layout_mode = 'force';

%% =========================================================
% Drag-snap grid settings
%   drag_snap_grid : if true, nodes snap to grid points when dragged
%   drag_grid_step : grid spacing. If [], determined automatically from the axis range.
%% =========================================================
drag_snap_grid = true;
drag_grid_step = [];     % e.g. set to 0.5 to fix it. [] = automatic.

%% --- Load data ---
opt = detectImportOptions(fullfile(codexDir, 'connections_no_threshold.csv'));
opt = setvartype(opt,'pre_root_id','int64');
opt = setvartype(opt,'post_root_id','int64');
FAFBConnections = readtable(fullfile(codexDir, 'connections_no_threshold.csv'), opt);
FAFBConnections(FAFBConnections.syn_count < 5, :) = [];

opt = detectImportOptions(fullfile(codexDir, 'consolidated_cell_types.csv'));
opt = setvartype(opt,'root_id','int64');
FAFBConsolidated_type = readtable(fullfile(codexDir, 'consolidated_cell_types.csv'), opt);

% --- classification.csv (side information) ---
opt = detectImportOptions(fullfile(codexDir, 'classification.csv'));
opt = setvartype(opt,'root_id','int64');
FAFBClassification = readtable(fullfile(codexDir, 'classification.csv'), opt);

% root_id -> side ('left'/'right'/'center'/'') map
side_raw = lower(strtrim(string(FAFBClassification.side)));
rid2side = containers.Map(num2cell(int64(FAFBClassification.root_id)), ...
                          cellstr(side_raw));

% Remove UNASGD (unassigned neuropil)
FAFBConnections = FAFBConnections( ...
    ~ismember(FAFBConnections.neuropil, exclude_neuropils), :);

%% --- Check base types exist ---
sim_types = sim_types(:);
miss = sim_types(~ismember(sim_types, FAFBConsolidated_type.primary_type));
if ~isempty(miss)
    warning('Types in sim_types not present in consolidated_cell_types: %s', ...
        strjoin(miss, ', '));
end

%% --- Per-side node definition (Type_L / Type_R) ---
sides_use   = {'L','R'};
side_label  = struct('L','left','R','right');

all_nodes  = {};        % split node names (e.g. 'MeMe_e01_L')
node_base  = {};        % base type of each node
node_side  = {};        % side of each node ('L'/'R')
type_to_root_ids  = containers.Map();

n_total_cells   = 0;
n_missing_side  = 0;
missing_by_type = containers.Map();

fprintf('--- Side split (per base type) ---\n');
for i = 1:numel(sim_types)
    t = sim_types{i};
    type_root_ids = unique( ...
        FAFBConsolidated_type.root_id( ...
        strcmp(FAFBConsolidated_type.primary_type, t)));

    side_of_cell = repmat({''}, numel(type_root_ids), 1);
    for k = 1:numel(type_root_ids)
        if isKey(rid2side, type_root_ids(k))
            side_of_cell{k} = rid2side(type_root_ids(k));
        end
    end

    is_left  = strcmp(side_of_cell, 'left');
    is_right = strcmp(side_of_cell, 'right');
    is_lr    = is_left | is_right;

    n_total_cells  = n_total_cells  + numel(type_root_ids);
    n_miss_t       = sum(~is_lr);
    n_missing_side = n_missing_side + n_miss_t;
    missing_by_type(t) = n_miss_t;

    for sidx = 1:numel(sides_use)
        sd  = sides_use{sidx};
        msk = strcmp(side_of_cell, side_label.(sd));
        if ~any(msk), continue; end
        node_name = [t '_' sd];
        all_nodes{end+1,1} = node_name; %#ok<SAGROW>
        node_base{end+1,1} = t;         %#ok<SAGROW>
        node_side{end+1,1} = sd;        %#ok<SAGROW>
        type_to_root_ids(node_name) = type_root_ids(msk);
    end

    fprintf('   %-10s : total=%-4d  L=%-4d  R=%-4d  no-side=%d\n', ...
        t, numel(type_root_ids), sum(is_left), sum(is_right), n_miss_t);
end
fprintf('Total split nodes : %d\n', numel(all_nodes));

%% --- Warn about missing side ---
if n_total_cells > 0
    miss_frac = n_missing_side / n_total_cells;
else
    miss_frac = 0;
end
fprintf('Neurons without left/right side : %d / %d (%.1f%%)\n', ...
    n_missing_side, n_total_cells, miss_frac*100);
if miss_frac >= side_missing_warn_frac
    warning(['Many neurons have no left/right side: %d/%d (%.1f%%) >= %.0f%%.\n' ...
             '       These neurons are excluded from the graph.'], ...
        n_missing_side, n_total_cells, miss_frac*100, side_missing_warn_frac*100);
end
fprintf('\n');

%% --- pre_root_id -> nt_type mapping ---
[u_pre_all, ia_pre] = unique(FAFBConnections.pre_root_id);
nt_per_pre = FAFBConnections.nt_type(ia_pre);
pre2nt = containers.Map(num2cell(u_pre_all), nt_per_pre);

%% --- Within-type (=split node) NT consistency check + remove minority cells ---
fprintf('--- NT consistency check (minority-NT cells -> removed) ---\n');
n_removed_total_nt = 0;
for i = 1:numel(all_nodes)
    tp = all_nodes{i};
    type_root_ids = type_to_root_ids(tp);
    if isempty(type_root_ids), continue; end

    nt_per_cell = repmat({''}, numel(type_root_ids), 1);
    for k = 1:numel(type_root_ids)
        if isKey(pre2nt, type_root_ids(k))
            nt_per_cell{k} = pre2nt(type_root_ids(k));
        end
    end
    has_nt = ~cellfun(@isempty, nt_per_cell);
    if ~any(has_nt), continue; end

    nt_known = nt_per_cell(has_nt);
    uNT_t    = unique(nt_known);
    counts   = zeros(numel(uNT_t),1);
    for u = 1:numel(uNT_t)
        counts(u) = sum(strcmp(nt_known, uNT_t{u}));
    end
    [~, max_idx] = max(counts);
    dom_NT_t = uNT_t{max_idx};
    n_tot    = sum(has_nt);

    minority_mask = has_nt & ~strcmp(nt_per_cell, dom_NT_t);
    n_min = sum(minority_mask);
    if n_min == 0, continue; end

    min_frac = n_min / n_tot;
    if min_frac >= nt_minority_warn_frac
        fprintf('   [WARN] %s : NT heterogeneous (minority=%.1f%%) -> KEEP ALL\n', ...
            tp, min_frac*100);
        continue;
    end

    keep_mask = ~minority_mask;
    type_to_root_ids(tp) = type_root_ids(keep_mask);
    n_removed_total_nt = n_removed_total_nt + n_min;
end
fprintf('Total minority-NT cells removed across all nodes: %d\n\n', n_removed_total_nt);

%% --- fraction mode: precompute total input syn per post_root_id ---
if strcmpi(weight_mode, 'fraction')
    [uPostAll, ~, ic_all] = unique(FAFBConnections.post_root_id);
    totalIn_per_post      = accumarray(ic_all, FAFBConnections.syn_count, [], @sum);
    totalInMap = containers.Map(num2cell(uPostAll), num2cell(totalIn_per_post));
else
    totalInMap = containers.Map('KeyType','int64','ValueType','double');
end

%% --- Build edges ---
src_list = {}; tgt_list = {};
abs_w = []; sgn = []; dom_nt = {}; n_pair = [];
mean_n_src_per_tgt = []; mean_n_tgt_per_src = [];
n_src_cells_all = []; n_tgt_cells_all = [];
n_src_connected = []; n_tgt_connected = [];
total_syn_all   = []; mean_syn_per_pair = [];

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

        % ===== Edge weight =====
        switch lower(weight_mode)
            case 'absolute'
                w = mean(pairSyn);
            case 'fraction'
                post_root_ids_of_pair = uPair(:,2);
                fracs = nan(n_existing_pairs, 1);
                for kk = 1:n_existing_pairs
                    rid_b = post_root_ids_of_pair(kk);
                    if isKey(totalInMap, rid_b)
                        T = totalInMap(rid_b);
                        if T > 0, fracs(kk) = pairSyn(kk) / T; end
                    end
                end
                w = mean(fracs, 'omitnan');
                if isnan(w), w = 0; end
            otherwise
                error('Unknown weight_mode: %s', weight_mode);
        end

        mean_src_per_tgt = n_existing_pairs / n_post_cells;
        mean_tgt_per_src = n_existing_pairs / n_pre_cells;
        n_src_conn = numel(unique(rows.pre_root_id));
        n_tgt_conn = numel(unique(rows.post_root_id));

        % ===== dominant neurotransmitter (-> sign) =====
        nts = rows.nt_type;
        [uNT, ~, icnt] = unique(nts);
        ntTotals = accumarray(icnt, rows.syn_count, [], @sum);
        [~, idxNT] = max(ntTotals);
        domNT = uNT{idxNT};
        s = 1;
        if ismember(domNT, inhib_nts), s = -1; end

        src_list{end+1,1} = src;  tgt_list{end+1,1} = tgt; %#ok<SAGROW>
        abs_w(end+1,1) = w;  sgn(end+1,1) = s; %#ok<SAGROW>
        dom_nt{end+1,1} = domNT;  n_pair(end+1,1) = n_existing_pairs; %#ok<SAGROW>
        mean_n_src_per_tgt(end+1,1) = mean_src_per_tgt; %#ok<SAGROW>
        mean_n_tgt_per_src(end+1,1) = mean_tgt_per_src; %#ok<SAGROW>
        n_src_cells_all(end+1,1) = n_pre_cells;  n_tgt_cells_all(end+1,1) = n_post_cells; %#ok<SAGROW>
        n_src_connected(end+1,1) = n_src_conn;   n_tgt_connected(end+1,1) = n_tgt_conn; %#ok<SAGROW>
        total_syn_all(end+1,1) = total_syn;  mean_syn_per_pair(end+1,1) = mean(pairSyn); %#ok<SAGROW>
    end
end

%% --- Remove sparse type-pair edges ---
total_input_per_type  = containers.Map();
total_output_per_type = containers.Map();
for i = 1:numel(all_nodes)
    tp = all_nodes{i};
    type_root_ids = type_to_root_ids(tp);
    if isempty(type_root_ids)
        total_input_per_type(tp)  = 0;
        total_output_per_type(tp) = 0;
    else
        msk_in  = ismember(FAFBConnections.post_root_id, type_root_ids);
        msk_out = ismember(FAFBConnections.pre_root_id,  type_root_ids);
        total_input_per_type(tp)  = sum(FAFBConnections.syn_count(msk_in));
        total_output_per_type(tp) = sum(FAFBConnections.syn_count(msk_out));
    end
end

keep_edge_by_frac = true(numel(src_list), 1);
for k = 1:numel(src_list)
    total_out_src = total_output_per_type(src_list{k});
    total_in_tgt  = total_input_per_type(tgt_list{k});
    if total_out_src > 0, frac_out = total_syn_all(k) / total_out_src; else, frac_out = 0; end
    if total_in_tgt  > 0, frac_in  = total_syn_all(k) / total_in_tgt;  else, frac_in  = 0; end
    if max(frac_in, frac_out) < edge_min_frac_of_tgt_input
        keep_edge_by_frac(k) = false;
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

%% =========================================================
% Figure 1 (panel S5A) : L/R-split connection digraph  (absolute weight)
%% =========================================================
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

% --- Remove isolates (optional) ---
if drop_isolated
    deg_total = indegree(G) + outdegree(G);
    iso_mask  = (deg_total == 0) & ~ismember(G.Nodes.Name, ...
        strcat(repelem(center_types(:),numel(sides_use)), '_', repmat(sides_use(:),numel(center_types),1)));
    if any(iso_mask)
        G = rmnode(G, G.Nodes.Name(iso_mask));
    end
end

% --- node category (base type / side) ---
node_names = G.Nodes.Name;
nN = numel(node_names);
[~, loc] = ismember(node_names, all_nodes);
nbase = node_base(loc);  nside = node_side(loc);
is_center = ismember(nbase, center_types);
is_left   = strcmp(nside, 'L');
is_right  = strcmp(nside, 'R');

% --- Positions (custom_lr) ---
switch lower(layout_mode)
    case 'custom_lr'
        x_left = -2.5; x_right = 2.5;
        nBase  = numel(sim_types);
        y_span = max(1.2, 0.8 * nBase);
        if nBase == 1, y_levels = 0; else, y_levels = linspace(1, -1, nBase) * y_span; end
        X = zeros(nN,1); Y = zeros(nN,1);
        for i = 1:nN
            bi = find(strcmp(sim_types, nbase{i}), 1);
            if isempty(bi), bi = 1; end
            Y(i) = y_levels(bi);
            if is_left(i), X(i) = x_left; else, X(i) = x_right; end
        end
    otherwise
        X = []; Y = [];
end
use_custom_xy = ~isempty(X);

% --- visualization ---
fig = figure(1);
set(fig, 'Color','w', 'Position',[100 100 1200 900], ...
    'Name','AllSimNeurons connection graph (L/R split)');

if height(G.Edges) > 0
    max_abs = max(abs(G.Edges.Weight));  if max_abs == 0, max_abs = 1; end
    lw = 0.5 + 5.5 * abs(G.Edges.Weight) / max_abs;
    ec = repmat([0.85 0.325 0.098], height(G.Edges), 1);
    ec(G.Edges.Weight < 0, :) = repmat([0.00 0.447 0.741], sum(G.Edges.Weight < 0), 1);
else
    lw = 1;  ec = [0 0 0];
end

nc = repmat([0.85 0.85 0.85], numnodes(G), 1);
nc(is_center, :) = repmat([0.10 0.10 0.10], sum(is_center), 1);
ns = 7 * ones(numnodes(G), 1);
ns(is_center) = 12;

if use_custom_xy
    p = plot(G, 'XData', X, 'YData', Y, ...
        'LineWidth', lw, 'EdgeColor', ec, 'EdgeAlpha', 0.85, ...
        'NodeColor', nc, 'MarkerSize', ns, ...
        'NodeFontSize', 11, 'NodeFontWeight', 'bold', ...
        'ArrowSize', 10, 'ArrowPosition', 0.85, 'Interpreter', 'none');
else
    p = plot(G, 'Layout', layout_mode, ...
        'LineWidth', lw, 'EdgeColor', ec, 'EdgeAlpha', 0.85, ...
        'NodeColor', nc, 'MarkerSize', ns, ...
        'NodeFontSize', 11, 'NodeFontWeight', 'bold', ...
        'ArrowSize', 10, 'ArrowPosition', 0.85, 'Interpreter', 'none');
end

title({ ...
    sprintf('All-sim-neurons connection graph (L/R split)  |  layout: %s  |  weight\\_mode: %s', ...
        layout_mode, weight_mode), ...
    'edge width \propto mean syn over connected (pre,post) pairs   |   orange = excitatory,   blue = inhibitory (GABA / GLUT)'}, ...
    'Interpreter', 'tex');
axis off

if use_custom_xy
    yl = ylim;
    text(-2.5, yl(2)+0.05*range(yl), 'LEFT',  'HorizontalAlignment','center', ...
        'FontWeight','bold','FontSize',13,'Color',[0.2 0.2 0.2]);
    text( 2.5, yl(2)+0.05*range(yl), 'RIGHT', 'HorizontalAlignment','center', ...
        'FontWeight','bold','FontSize',13,'Color',[0.2 0.2 0.2]);
end

hold on
hp = plot(NaN, NaN, '-', 'Color', [0.85 0.325 0.098], 'LineWidth', 3);
hn = plot(NaN, NaN, '-', 'Color', [0.00 0.447 0.741], 'LineWidth', 3);

% --- Edge data-tip (click an edge midpoint to show info) ---
xd = p.XData;  yd = p.YData;
EN_g = G.Edges.EndNodes;
nE_g = height(G.Edges);
mid_x = nan(nE_g, 1);  mid_y = nan(nE_g, 1);
src_label = cell(nE_g, 1);  tgt_label = cell(nE_g, 1);
srcNodeIdx = ones(nE_g, 1);  tgtNodeIdx = ones(nE_g, 1);
for ee = 1:nE_g
    if iscell(EN_g), s_nm = EN_g{ee,1}; t_nm = EN_g{ee,2};
    else,            s_nm = char(EN_g(ee,1)); t_nm = char(EN_g(ee,2)); end
    s_idx = find(strcmp(G.Nodes.Name, s_nm), 1);
    t_idx = find(strcmp(G.Nodes.Name, t_nm), 1);
    if isempty(s_idx) || isempty(t_idx), continue; end
    srcNodeIdx(ee) = s_idx;  tgtNodeIdx(ee) = t_idx;
    mid_x(ee) = (xd(s_idx) + xd(t_idx)) / 2;
    mid_y(ee) = (yd(s_idx) + yd(t_idx)) / 2;
    src_label{ee} = s_nm;  tgt_label{ee} = t_nm;
end

h_mid = scatter(mid_x, mid_y, 24, [0.3 0.3 0.3], 'filled', ...
    'MarkerFaceAlpha', 0.35, 'MarkerEdgeColor', 'none');
h_mid.DataTipTemplate.DataTipRows = [ ...
    dataTipTextRow('src',                  src_label); ...
    dataTipTextRow('tgt',                  tgt_label); ...
    dataTipTextRow('weight (signed)',      G.Edges.Weight); ...
    dataTipTextRow('|w|',                  G.Edges.AbsWeight); ...
    dataTipTextRow('dom NT',               G.Edges.DomNT); ...
    dataTipTextRow('N pair (connected)',   G.Edges.NPair); ...
    dataTipTextRow('total syn',            G.Edges.TotalSyn); ...
    dataTipTextRow('mean syn / pair',      G.Edges.MeanSynPerPair); ...
    dataTipTextRow('fanIn  (src/postCell)',G.Edges.MeanSrcPerTgt); ...
    dataTipTextRow('fanOut (tgt/preCell)', G.Edges.MeanTgtPerSrc); ...
    dataTipTextRow('N src cells (all)',    G.Edges.NSrcCellsAll); ...
    dataTipTextRow('N tgt cells (all)',    G.Edges.NTgtCellsAll); ...
    dataTipTextRow('N src connected',      G.Edges.NSrcConnected); ...
    dataTipTextRow('N tgt connected',      G.Edges.NTgtConnected)];

hm = plot(NaN, NaN, '.', 'Color', [0.3 0.3 0.3], 'MarkerSize', 14);
legend([hp hn hm], ...
    {'Excitatory', 'Inhibitory (GABA/GLUT)', 'Click an edge midpoint for info'}, ...
    'Location', 'southoutside', 'Orientation', 'horizontal', 'Box', 'off');

dcm = datacursormode(fig);
dcm.SnapToDataVertex = 'on';
dcm.Enable = 'off';

% --- Snap grid + node dragging ---
ax_g = ancestor(p, 'axes');
drawnow;
if drag_snap_grid
    gstep = drag_grid_step;
    if isempty(gstep)
        span  = max(range(xlim(ax_g)), range(ylim(ax_g)));
        gstep = niceStep(span / 20);
    end
    drawSnapGrid(ax_g, gstep);
else
    gstep = [];
end
setappdata(fig, 'dcm',      dcm);
setappdata(fig, 'dragArgs', {p, h_mid, srcNodeIdx, tgtNodeIdx, gstep});
enableNodeDragging(fig, p, h_mid, srcNodeIdx, tgtNodeIdx, gstep);
uicontrol(fig, 'Style', 'togglebutton', 'String', 'Mode: Drag nodes', ...
    'Units', 'normalized', 'Position', [0.01 0.95 0.20 0.04], ...
    'FontSize', 9, 'Value', 0, 'Callback', @onToggleMode);

%% =========================================================
% Figure 2 (panel S5B) : connection MATRIX (= figure 1 absolute edge weights)
%   row = source (pre), column = target (post), cell = signed mean syn / pair.
%% =========================================================
N = numel(all_nodes);
M = nan(N, N);  Npair = nan(N, N);  TotSyn = nan(N, N);
NsrcCon = nan(N, N);  NtgtCon = nan(N, N);
for k = 1:numel(src_list)
    i = find(strcmp(all_nodes, src_list{k}), 1);
    j = find(strcmp(all_nodes, tgt_list{k}), 1);
    if ~isempty(i) && ~isempty(j)
        M(i, j) = signed_w(k);  Npair(i, j) = n_pair(k);
        TotSyn(i, j) = total_syn_all(k);
        NsrcCon(i, j) = n_src_connected(k);  NtgtCon(i, j) = n_tgt_connected(k);
    end
end

fig2 = figure(2);
set(fig2, 'Color','w', 'Position',[120 120 950 820], ...
    'Name','Connection matrix (absolute weight, signed)');
ax3 = axes(fig2);
him = imagesc(ax3, M);
set(him, 'AlphaData', ~isnan(M));
set(ax3, 'Color', 'w');
colormap(ax3, diverging_bwo(256));
maxabs = max(abs(M(:)), [], 'omitnan');
if isempty(maxabs) || maxabs == 0 || isnan(maxabs), maxabs = 1; end
caxis(ax3, [-maxabs maxabs]);
cb = colorbar(ax3);
cb.Label.String = 'signed mean syn / connected (pre,post) pair';

ax3.XTick = 1:N;  ax3.YTick = 1:N;
ax3.XTickLabel = all_nodes;  ax3.YTickLabel = all_nodes;
ax3.XTickLabelRotation = 45;
ax3.TickLabelInterpreter = 'none';
ax3.TickDir = 'out';
xlabel(ax3, 'target (post)');
ylabel(ax3, 'source (pre)');
axis(ax3, 'image');
title(ax3, { ...
    'Connection matrix  (= figure 1 edge weights, absolute mode)', ...
    'value = signed mean syn over connected (pre,post) pairs   |   + excit / - inhib (GABA/GLUT)'}, ...
    'Interpreter','none');

for i = 1:N
    for j = 1:N
        if isnan(M(i,j)), continue; end
        v = M(i,j);
        if abs(v) > 0.6*maxabs, txtcol = 'w'; else, txtcol = 'k'; end
        text(ax3, j, i, sprintf('%.1f', v), ...
            'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
            'FontSize', 8, 'Color', txtcol);
    end
end

hold(ax3, 'on');
for g = 0.5:1:(N+0.5)
    line(ax3, [0.5 N+0.5], [g g], 'Color', [0.85 0.85 0.85], 'LineWidth', 0.5);
    line(ax3, [g g], [0.5 N+0.5], 'Color', [0.85 0.85 0.85], 'LineWidth', 0.5);
end
hold(ax3, 'off');

setappdata(fig2, 'mat_M', M);  setappdata(fig2, 'mat_Npair', Npair);
setappdata(fig2, 'mat_TotSyn', TotSyn);
setappdata(fig2, 'mat_NsrcCon', NsrcCon);  setappdata(fig2, 'mat_NtgtCon', NtgtCon);
setappdata(fig2, 'mat_nodes', all_nodes);
dcm3 = datacursormode(fig2);
dcm3.SnapToDataVertex = 'off';
dcm3.UpdateFcn = @matrixDataTip;
dcm3.Enable = 'on';
text(ax3, 0.5, N+1.1, 'Click a cell to show connection count (N pairs) and synapse info', ...
    'HorizontalAlignment','left', 'FontSize', 9, 'Color', [0.35 0.35 0.35], ...
    'Interpreter','none');

%% =========================================================
% Local functions
%% =========================================================
function enableNodeDragging(fig, p, h_mid, srcNodeIdx, tgtNodeIdx, gstep)
    if nargin < 6, gstep = []; end
    ax      = ancestor(p, 'axes');
    dragIdx = [];
    set(fig, 'WindowButtonDownFcn',   @onDown);
    set(fig, 'WindowButtonMotionFcn', @onMove);
    set(fig, 'WindowButtonUpFcn',     @onUp);
    set(fig, 'Pointer', 'arrow');
    function onDown(~,~)
        cp = get(ax, 'CurrentPoint');
        x  = cp(1,1);  y = cp(1,2);
        xd = p.XData;  yd = p.YData;
        xl = xlim(ax); yl = ylim(ax);
        rx = max(range(xl), eps);  ry = max(range(yl), eps);
        dd = sqrt(((xd - x)/rx).^2 + ((yd - y)/ry).^2);
        [dmin, idx] = min(dd);
        if dmin < 0.04
            dragIdx = idx;
            set(fig, 'Pointer', 'fleur');
        else
            dragIdx = [];
        end
    end
    function onMove(~,~)
        if isempty(dragIdx), return; end
        cp = get(ax, 'CurrentPoint');
        nx = cp(1,1);  ny = cp(1,2);
        if ~isempty(gstep) && gstep > 0
            nx = round(nx / gstep) * gstep;
            ny = round(ny / gstep) * gstep;
        end
        xd = p.XData;  yd = p.YData;
        xd(dragIdx) = nx;  yd(dragIdx) = ny;
        p.XData = xd;  p.YData = yd;
        if ~isempty(h_mid) && isvalid(h_mid)
            h_mid.XData = (xd(srcNodeIdx) + xd(tgtNodeIdx)) / 2;
            h_mid.YData = (yd(srcNodeIdx) + yd(tgtNodeIdx)) / 2;
        end
    end
    function onUp(~,~)
        dragIdx = [];
        set(fig, 'Pointer', 'arrow');
    end
end

function onToggleMode(btn, ~)
    fig  = ancestor(btn, 'figure');
    dcm  = getappdata(fig, 'dcm');
    args = getappdata(fig, 'dragArgs');
    if btn.Value
        set(fig, 'WindowButtonDownFcn', '', 'WindowButtonMotionFcn', '', ...
                 'WindowButtonUpFcn', '', 'Pointer', 'arrow');
        dcm.Enable = 'on';
        btn.String = 'Mode: Data tips';
    else
        dcm.Enable = 'off';
        enableNodeDragging(fig, args{:});
        btn.String = 'Mode: Drag nodes';
    end
end

function drawSnapGrid(ax, gstep)
    xl = xlim(ax);  yl = ylim(ax);
    pad = 0.6;
    xmin = xl(1) - pad*range(xl);  xmax = xl(2) + pad*range(xl);
    ymin = yl(1) - pad*range(yl);  ymax = yl(2) + pad*range(yl);
    gx = gstep * (floor(xmin/gstep) : ceil(xmax/gstep));
    gy = gstep * (floor(ymin/gstep) : ceil(ymax/gstep));
    col_minor = [0.90 0.90 0.90];
    col_axis  = [0.70 0.70 0.78];
    hold(ax, 'on');
    for k = 1:numel(gx)
        c = col_minor;  w = 0.5;
        if abs(gx(k)) < gstep/2, c = col_axis; w = 1.2; end
        line(ax, [gx(k) gx(k)], [ymin ymax], 'Color', c, 'LineWidth', w, 'HandleVisibility', 'off');
    end
    for k = 1:numel(gy)
        c = col_minor;  w = 0.5;
        if abs(gy(k)) < gstep/2, c = col_axis; w = 1.2; end
        line(ax, [xmin xmax], [gy(k) gy(k)], 'Color', c, 'LineWidth', w, 'HandleVisibility', 'off');
    end
    ch = ax.Children;
    isGrid = false(numel(ch),1);
    for k = 1:numel(ch)
        isGrid(k) = isa(ch(k), 'matlab.graphics.primitive.Line') && ...
                    strcmp(get(ch(k), 'HandleVisibility'), 'off');
    end
    ax.Children = [ch(~isGrid); ch(isGrid)];
end

function g = niceStep(raw)
    if raw <= 0 || ~isfinite(raw), g = 1; return; end
    p10  = 10^floor(log10(raw));
    frac = raw / p10;
    if     frac < 2,  m = 1;
    elseif frac < 5,  m = 2;
    else,             m = 5;
    end
    g = m * p10;
end

function cmap = diverging_bwo(n)
    if nargin < 1, n = 256; end
    blue   = [0.00 0.447 0.741];
    white  = [1.00 1.000 1.000];
    orange = [0.85 0.325 0.098];
    h = floor(n/2);
    low  = [linspace(blue(1),  white(1),  h)', linspace(blue(2),  white(2),  h)', linspace(blue(3),  white(3),  h)'];
    high = [linspace(white(1), orange(1), n-h)', linspace(white(2), orange(2), n-h)', linspace(white(3), orange(3), n-h)'];
    cmap = [low; high];
end

function txt = matrixDataTip(~, event)
    fig   = ancestor(event.Target, 'figure');
    M       = getappdata(fig, 'mat_M');
    Npair   = getappdata(fig, 'mat_Npair');
    TotSyn  = getappdata(fig, 'mat_TotSyn');
    NsrcCon = getappdata(fig, 'mat_NsrcCon');
    NtgtCon = getappdata(fig, 'mat_NtgtCon');
    nodes   = getappdata(fig, 'mat_nodes');
    pos = event.Position;
    j = round(pos(1));  i = round(pos(2));
    N = numel(nodes);
    if i < 1 || i > N || j < 1 || j > N
        txt = {'(out of range)'};  return;
    end
    if isnan(M(i,j))
        txt = {sprintf('%s  ->  %s', nodes{i}, nodes{j}), 'no connection'};  return;
    end
    txt = { ...
        sprintf('src (pre) : %s', nodes{i}), ...
        sprintf('tgt (post): %s', nodes{j}), ...
        sprintf('weight (signed) : %.2f', M(i,j)), ...
        sprintf('N connected pairs : %d', Npair(i,j)), ...
        sprintf('total syn : %d', round(TotSyn(i,j))), ...
        sprintf('N src connected : %d', NsrcCon(i,j)), ...
        sprintf('N tgt connected : %d', NtgtCon(i,j))};
end
