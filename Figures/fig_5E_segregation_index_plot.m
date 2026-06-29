%% 1. Load data and initialize
clear all; close all; clc

baseDir = fileparts(fileparts(mfilename('fullpath')));  % resolved from this script's location (portable; no absolute paths)

% NPI tables (type, root_id, Right_PostPI, Right_PrePI) from
% Figures/fig_1D_E_postPI_prePI_right_neurons.m
load(fullfile(baseDir, 'Processed_Data', 'right_neurons_thr0.mat'), ...
    'RightFFP_NPIs', 'RightFBP_NPIs', 'RightBDP_real_NPIs')

% Per-neuron segregation index from Data_Processing/s13_segregation_index.ipynb,
% attached to each NPI table by root_id.
RightFFP_NPIs = attach_segregation_index(RightFFP_NPIs, ...
    fullfile(baseDir, 'Processed_Data', 'right_FFP_segregation_index.csv'));
RightFBP_NPIs = attach_segregation_index(RightFBP_NPIs, ...
    fullfile(baseDir, 'Processed_Data', 'right_FBP_segregation_index.csv'));
RightBDP_NPIs = attach_segregation_index(RightBDP_real_NPIs, ...
    fullfile(baseDir, 'Processed_Data', 'right_BDP_segregation_index.csv'));

%% 2. Per-type mean segregation index
Type_FFP = summarize_by_type(RightFFP_NPIs);
Type_FBP = summarize_by_type(RightFBP_NPIs);
Type_BDP = summarize_by_type(RightBDP_NPIs);

FFP_segIdx = Type_FFP.mean_segregation_index;
FBP_segIdx = Type_FBP.mean_segregation_index;
BDP_segIdx = Type_BDP.mean_segregation_index;

% Example bidirectional neurons (per neuron)
LC9_segIdx  = RightBDP_NPIs.segregation_index(strcmp(RightBDP_NPIs.type, 'LC9'));
LT43_segIdx = RightBDP_NPIs.segregation_index(strcmp(RightBDP_NPIs.type, 'LT43'));
LT52_segIdx = RightBDP_NPIs.segregation_index(strcmp(RightBDP_NPIs.type, 'LT52'));

%% 3. Figure 1 (panel 5E): segregation-index box plot
data = {LC9_segIdx, LT43_segIdx, LT52_segIdx, BDP_segIdx, FFP_segIdx, FBP_segIdx};
group = repelem(1:length(data), cellfun(@length, data));
combinedData = vertcat(data{:});

figure(1); set(gcf,'Color','w')
boxplot(combinedData, group, 'Labels', {'LC9','LT43','LT52','BDP','FFP','FBP'}, 'Notch','off','Symbol','');
title('segIdx')
set(gca,'Box','off','TickDir','out')
ylim([0 1])

% Pairwise comparison (BDP / FFP / FBP)
groupNames = {'BDP', 'FFP', 'FBP'};
dataPairs = {
    BDP_segIdx, FFP_segIdx;
    BDP_segIdx, FBP_segIdx;
    FFP_segIdx, FBP_segIdx;
    };
groupPairs = {[1, 2]; [1, 3]; [2, 3]};

alpha = 0.05;
numComparisons = size(dataPairs, 1);
correctedAlpha = alpha / numComparisons;

fprintf('=== segIdx Pairwise Comparison (Bonferroni corrected alpha = %.4f) ===\n\n', correctedAlpha);
for i = 1:numComparisons
    d1 = dataPairs{i, 1};
    d2 = dataPairs{i, 2};
    g1 = groupNames{groupPairs{i}(1)};
    g2 = groupNames{groupPairs{i}(2)};

    p = ranksum(d1, d2);                 % Wilcoxon rank-sum test
    p_dominance  = mean(d1 > d2', 'all'); % probabilistic dominance
    p_dominance2 = mean(d1 < d2', 'all');

    fprintf('%s vs %s:\n', g1, g2);
    fprintf('  p-value = %.4f --> %s\n', p, ...
        ternary(p < correctedAlpha, 'Significant', 'Not significant'));
    fprintf('  P(%s > %s) = %.3f\n', g1, g2, p_dominance);
    fprintf('  P(%s < %s) = %.3f\n\n', g1, g2, p_dominance2);
end

%% Local functions
function T = attach_segregation_index(T, csvPath)
% Read a (root_id, segregation_index) CSV and attach it to table T by root_id.
opts = detectImportOptions(csvPath);
opts = setvartype(opts, 'root_id', 'int64');
S = readtable(csvPath, opts);

T.segregation_index = nan(height(T), 1);
[tf, loc] = ismember(T.root_id, S.root_id);
T.segregation_index(tf) = S.segregation_index(loc(tf));
end

function Tsum = summarize_by_type(T)
% Mean segregation index per cell type.
[types, ~, ic] = unique(T.type);
Tsum = table(types, 'VariableNames', {'type'});
for i = 1:numel(types)
    Tsum.mean_segregation_index(i) = mean(T.segregation_index(ic == i), 'omitnan');
end
end

function out = ternary(cond, valTrue, valFalse)
if cond
    out = valTrue;
else
    out = valFalse;
end
end
