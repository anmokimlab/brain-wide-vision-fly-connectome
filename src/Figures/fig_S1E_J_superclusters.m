%% Supercluster output composition pie charts (panels S1E–J)
clear all; close all; clc

baseDir = fileparts(fileparts(mfilename('fullpath')));  % resolved from this script's location (portable; no absolute paths)

% Per-supercluster post-synaptic targets, from Figures/fig_3C_D_E_F_G_FFP_clustering_plot.m
load(fullfile(baseDir, 'Processed_Data', 'FFP_supercluster_targets.mat'))   % SuperCluster_targets
nSC = numel(SuperCluster_targets);

%% 1. Post-synaptic neuron type distribution (by neuron count)
Thr = 0.05;
for i = 1:nSC
    items = collapse_small(SuperCluster_targets(i).post_type_counts, SuperCluster_targets(i).n_post_neurons * Thr);
    draw_pie(items);
end

%% 2. Post-synaptic neuropil distribution (by synapse count)
Thr = 0.05;
for i = 1:nSC
    items = SuperCluster_targets(i).post_neuropil_syn;
    items = collapse_small(items, sum(cell2mat(items(:,2))) * Thr);
    draw_pie(items);
end

%% 3. Post-synaptic neuron type distribution (by synapse count)
Thr = 0.01;
for i = 1:nSC
    items = SuperCluster_targets(i).post_type_syn;
    items = collapse_small(items, sum(cell2mat(items(:,2))) * Thr);
    draw_pie(items);
end

%% Local functions
function items = collapse_small(items, thr)
% Move {label, value} rows with value < thr into a single 'Others' row.
vals = cell2mat(items(:,2));
keep = vals >= thr;
items = [items(keep,:); {'Others', sum(vals(~keep))}];
end

function draw_pie(items)
labels = strrep(items(:,1), '_', ' ');
figure('Color','w');
p = piechart(cell2mat(items(:,2)), labels, 'FaceAlpha',1, 'EdgeColor',[0.50 0.50 0.50]);
p.LabelStyle = "namedata";
p.FontSize = 12;
colororder("gem12")
end
