%% s04_compute_graph_centrality
% Build the whole-brain connectivity graph from the Codex connection table and
% compute node centrality (PageRank, unweighted/weighted betweenness).
% Precursor data processing for Figure 2C, 2D, 2E.
clear all; close all; clc;

baseDir = fileparts(fileparts(mfilename('fullpath')));  % resolved from this script's location (portable; no absolute paths)

opt = detectImportOptions(fullfile(baseDir, 'Codex_Data', 'connections_no_threshold.csv'));
opt = setvartype(opt,'pre_root_id','int64');
opt = setvartype(opt,'post_root_id','int64');
FAFBConnections = readtable(fullfile(baseDir, 'Codex_Data', 'connections_no_threshold.csv'),opt);

FAFBConnections = sortrows(FAFBConnections, {'pre_root_id','post_root_id','syn_count'}, {'ascend','ascend','descend'});


% 2. Check whether each row has the same (pre, post) pair as the previous row
same_as_previous = [false; ...
    FAFBConnections.pre_root_id(2:end) == FAFBConnections.pre_root_id(1:end-1) & ...
    FAFBConnections.post_root_id(2:end) == FAFBConnections.post_root_id(1:end-1)];

% 3. Assign a unique group ID to each (pre, post) pair
group_id = cumsum(~same_as_previous);

% 4. Sum syn_count within each group
[G, pre_group] = findgroups(group_id);
syn_sum = splitapply(@sum, FAFBConnections.syn_count, G);

% 5. Keep only the first row of each group
first_in_group = [true; diff(group_id) ~= 0];
FAFBConnections = FAFBConnections(first_in_group, :);

% 6. Assign the summed syn_count
FAFBConnections.syn_count = syn_sum;
%% Make graph
rootIds = unique([FAFBConnections.pre_root_id; FAFBConnections.post_root_id]);
[~, preIdx] = ismember(FAFBConnections.pre_root_id, rootIds);
[~, postIdx] = ismember(FAFBConnections.post_root_id, rootIds);
weight = FAFBConnections.syn_count;

G = digraph(preIdx, postIdx, weight);

%% Compute centrality

% (1) PageRank: weights not used
pagerank_vals = centrality(G, 'pagerank');

% (2) Betweenness (unweighted)
betweenness_vals_unweighted = centrality(G, 'betweenness');

% (3) Betweenness (weighted: a larger weight is treated as a shorter distance)
edge_costs = 1 ./ (G.Edges.Weight + eps);  % inverse weight, avoid division by zero
betweenness_vals_weighted = centrality(G, 'betweenness', 'Cost', edge_costs);

%% Save result
save(fullfile(baseDir, 'Processed_Data', 'allgraph_thr0.mat'), "G", "postIdx", "preIdx", "rootIds", "weight", ...
    "pagerank_vals", "betweenness_vals_unweighted", "betweenness_vals_weighted");
