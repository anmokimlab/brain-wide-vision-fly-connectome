%% fig_2_FFP_optic_central (MCNS)
% MCNS analogue of the FAFB Figures/fig_2C_D_E_plot_FFP_optic_central.m.
% Compares the Optic / Feedforward (FFP) / Central neuron groups on per-type
% fan-in / fan-out statistics. (No graph-centrality panels: MCNS has no centrality.)
clear all; close all; clc

baseDir = fileparts(fileparts(mfilename('fullpath')));  % resolved from this script's location (portable; no absolute paths)

% Per-type fan-in/out summaries from MCNS/Data_Processing/s06_compute_FFP_optic_central_fan_in_out.m
load(fullfile(baseDir, 'Processed_Data', 'FFP_optic_central_fan_in_out.mat'))   % type_Optic_InOut, type_Central_InOut, type_RightFFP_InOut

%% Figure 2: Optic vs Feedforward vs Central box plots (with rank-sum tests)
% {field, ylabel, ylim, ytick}
metrics = {
    'InNeuronNumber',     'InNeuronNumber',     [0 2100], 0:500:2000;
    'InNeuronTypeNumber', 'InNeuronTypeNumber', [0 500],  0:100:500;
    'InNeuronRatio',      'InNeuronRatio',      [0 16],   0:5:20;
    'OutNeuronNumber',    'OutNeuronNumber',    [0 1500], [];
    'OutNeuronTypeNumber','OutNeuronType',      [0 500],  0:100:500;
    'OutNeuronRatio',     'OutNeuronRatio',     [0 16],   0:5:15;
};

for m = 1:size(metrics,1)
    f = metrics{m,1};
    d_optic   = cell2mat(type_Optic_InOut.(f));
    d_ff      = cell2mat(type_RightFFP_InOut.(f));
    d_central = cell2mat(type_Central_InOut.(f));
    draw_metric_boxplot(d_optic, d_ff, d_central, metrics{m,2}, metrics{m,3}, metrics{m,4});
end

%% Figure S2: per-type Top-10 bar charts (FFP / Optic / Central)
N_bar = 10;
% {field, title}
barMetrics = {
    'InNeuronNumber',      'In Neuron Number';
    'InNeuronTypeNumber',  'In Neuron Type Number';
    'OutNeuronNumber',     'Out Neuron Number';
    'OutNeuronTypeNumber', 'Out Neuron Type Number';
};
% {table, name, color}
barGroups = {
    type_RightFFP_InOut, 'FFP',     '#0072BD';
    type_Optic_InOut,    'Optic',   '#EDB120';
    type_Central_InOut,  'Central', '#7E2F8E';
};

for bm = 1:size(barMetrics,1)
    for g = 1:size(barGroups,1)
        T = barGroups{g,1};
        draw_top_bar(T.type, cell2mat(T.(barMetrics{bm,1})), N_bar, barGroups{g,3}, ...
            sprintf('%s %s', barMetrics{bm,2}, barGroups{g,2}));
    end
end

%% ===================== Local functions =====================
function draw_metric_boxplot(d_optic, d_ff, d_central, ylabelStr, yl, ytick)
% 3-group box plot (Optic / Feedforward / Central) with rank-sum stats and colored boxes.
data = {d_optic, d_ff, d_central};
group = repelem(1:3, cellfun(@numel, data));
combinedData = vertcat(data{:});

figure('Color','w');
boxplot(combinedData, group, 'Labels', {'Optic','Feedforward','Central'}, 'Notch','on','Symbol','');
ylabel(ylabelStr); ylim(yl);
set(gca,'Box','off','TickDir','out');
if ~isempty(ytick), set(gca,'YTick',ytick); end

% Color the boxes/medians (findobj returns them in reverse order: 1=Central, 2=FF, 3=Optic)
boxColors = [0.8500 0.3250 0.0980;   % box 1
             0.0000 0.4470 0.7410;   % box 2
             0.9290 0.6940 0.1250];  % box 3
ax = gca;
boxes   = findobj(ax,'Tag','Box');
medians = findobj(ax,'Tag','Median');
for k = 1:min(3,numel(boxes))
    boxes(k).Color   = boxColors(k,:);
    medians(k).Color = boxColors(k,:);
end

pairwise_ranksum(data, {'Optic','Feedforward','Central'}, ylabelStr);
end

function pairwise_ranksum(data, names, label)
% Wilcoxon rank-sum (Bonferroni-corrected) + probabilistic dominance, printed to console.
pairs = {[1 2],[1 3],[2 3]};
correctedAlpha = 0.05 / numel(pairs);
fprintf('=== %s Pairwise Comparison (Bonferroni corrected alpha = %.4f) ===\n\n', label, correctedAlpha);
for i = 1:numel(pairs)
    a = pairs{i}(1); b = pairs{i}(2);
    d1 = data{a}; d2 = data{b};
    p = ranksum(d1, d2);
    p_dominance  = mean(d1 > d2', 'all');
    p_dominance2 = mean(d1 < d2', 'all');
    fprintf('%s vs %s:\n', names{a}, names{b});
    fprintf('  p-value = %.4f --> %s\n', p, ternary(p < correctedAlpha, 'Significant', 'Not significant'));
    fprintf('  P(%s > %s) = %.3f\n', names{a}, names{b}, p_dominance);
    fprintf('  P(%s < %s) = %.3f\n\n', names{a}, names{b}, p_dominance2);
end
end

function draw_top_bar(types, values, n, color, titleStr)
% Top-n cell types by value (descending), one bar chart.
[sorted, idx] = sort(values, 'descend', 'MissingPlacement', 'last');
n = min(n, numel(sorted));
figure('Color','w','InnerPosition',[551,202,621,291],'OuterPosition',[543,194,637,385]);
bar(types(idx(1:n)), sorted(1:n), 'FaceColor', color, 'EdgeColor','flat');
set(gca,'Box','off','TickDir','out');
title(titleStr);
end

function out = ternary(cond, valTrue, valFalse)
if cond
    out = valTrue;
else
    out = valFalse;
end
end
