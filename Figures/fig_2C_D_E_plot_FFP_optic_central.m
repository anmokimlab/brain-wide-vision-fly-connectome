clear all; close all ;clc

baseDir = fileparts(fileparts(mfilename('fullpath')));  % resolved from this script's location (portable; no absolute paths)

% Per-type fan-in / fan-out summaries (from s05_compute_FFP_optic_central_fan_in_out.m)
load(fullfile(baseDir, 'Processed_Data', 'FFP_optic_central_fan_in_out.mat'))

%% 1. In neuron number

data1 = cell2mat(type_Optic_InOut.InNeuronNumber);
data2 = cell2mat(type_RightFFP_InOut.InNeuronNumber);
data3 = cell2mat(type_Central_InOut.InNeuronNumber);

% Define the data and groups (the three groups have different lengths)
data = {data1, data2, data3};
group = repelem(1:length(data), cellfun(@length, data));
combinedData = vertcat(data{:});
figure(1);set(gcf,'Color','w')
% Box plot
boxplot(combinedData, group, 'Labels', {'Optic', 'Feedforward', 'Central'},'Notch','on','Symbol','');

ylabel('InNeuronNumber');
ylim([0, 2000]);
set(gca,'Box','off','TickDir','out')

% Recolor the boxes and medians per group
ax = gca;
boxes = findobj(ax, 'Tag', 'Box');         % box handles
boxes(3).Color = [0.9290 0.6940 0.1250];   % Central
boxes(2).Color = [0 0.4470 0.7410];        % Feedforward
boxes(1).Color = [0.8500 0.3250 0.0980];   % Optic
medians = findobj(ax, 'Tag', 'Median');    % median handles
medians(3).Color = [0.9290 0.6940 0.1250];
medians(2).Color = [0 0.4470 0.7410];
medians(1).Color = [0.8500 0.3250 0.0980];

%%%% Nonparametric test
groupNames = {'Optic', 'Feedforward', 'Central'};

% Comparison pairs
dataPairs = {
    data1, data2;
    data1, data3;
    data2, data3
    };

groupPairs = {
    [1, 2];
    [1, 3];
    [2, 3]
    };

% Bonferroni-corrected significance level
alpha = 0.05;
numComparisons = size(dataPairs, 1);
correctedAlpha = alpha / numComparisons;

% Print results
fprintf('=== Number of In Neuron Pairwise Comparison (Bonferroni corrected alpha = %.4f) ===\n\n', correctedAlpha);

for i = 1:numComparisons
    d1 = dataPairs{i, 1};
    d2 = dataPairs{i, 2};
    g1 = groupNames{groupPairs{i}(1)};
    g2 = groupNames{groupPairs{i}(2)};

    % Wilcoxon rank-sum test (nonparametric, compares the medians of two groups)
    p = ranksum(d1, d2);

    % Probabilistic dominance (vectorized)
    p_dominance = mean(d1 > d2', 'all');
    p_dominance2 = mean(d1 < d2', 'all');

    % Print
    fprintf('%s vs %s:\n', g1, g2);
    fprintf('  p-value = %.4f --> %s\n', p, ...
        ternary(p < correctedAlpha, 'Significant', 'Not significant'));
    fprintf('  P(%s > %s) = %.3f\n', g1, g2, p_dominance);
    fprintf('  P(%s < %s) = %.3f\n\n', g1, g2, p_dominance2);
end

%% 2. In neuron type number

data1 = cell2mat(type_Optic_InOut.InNeuronTypeNumber);
data2 = cell2mat(type_RightFFP_InOut.InNeuronTypeNumber);
data3 = cell2mat(type_Central_InOut.InNeuronTypeNumber);

data = {data1, data2, data3};
group = repelem(1:length(data), cellfun(@length, data));
combinedData = vertcat(data{:});
figure(2);set(gcf,'Color','w')
boxplot(combinedData, group, 'Labels', {'Optic', 'Feedforward', 'Central'},'Notch','on','Symbol','');

ylabel('InNeuronTypeNumber');
ylim([0, 500]);
set(gca,'Box','off','TickDir','out')

ax = gca;
boxes = findobj(ax, 'Tag', 'Box');
boxes(3).Color = [0.9290 0.6940 0.1250];   % Central
boxes(2).Color = [0 0.4470 0.7410];        % Feedforward
boxes(1).Color = [0.8500 0.3250 0.0980];   % Optic
medians = findobj(ax, 'Tag', 'Median');
medians(3).Color = [0.9290 0.6940 0.1250];
medians(2).Color = [0 0.4470 0.7410];
medians(1).Color = [0.8500 0.3250 0.0980];

%%%% Nonparametric test
groupNames = {'Optic', 'Feedforward', 'Central'};

dataPairs = {
    data1, data2;
    data1, data3;
    data2, data3
    };

groupPairs = {
    [1, 2];
    [1, 3];
    [2, 3]
    };

alpha = 0.05;
numComparisons = size(dataPairs, 1);
correctedAlpha = alpha / numComparisons;

fprintf('=== Number of In Neuron Type Pairwise Comparison (Bonferroni corrected alpha = %.4f) ===\n\n', correctedAlpha);

for i = 1:numComparisons
    d1 = dataPairs{i, 1};
    d2 = dataPairs{i, 2};
    g1 = groupNames{groupPairs{i}(1)};
    g2 = groupNames{groupPairs{i}(2)};

    p = ranksum(d1, d2);

    p_dominance = mean(d1 > d2', 'all');
    p_dominance2 = mean(d1 < d2', 'all');

    fprintf('%s vs %s:\n', g1, g2);
    fprintf('  p-value = %.4f --> %s\n', p, ...
        ternary(p < correctedAlpha, 'Significant', 'Not significant'));
    fprintf('  P(%s > %s) = %.3f\n', g1, g2, p_dominance);
    fprintf('  P(%s < %s) = %.3f\n\n', g1, g2, p_dominance2);
end

%% 3. In neuron ratio

data1 = cell2mat(type_Optic_InOut.InNeuronRatio);
data2 = cell2mat(type_RightFFP_InOut.InNeuronRatio);
data3 = cell2mat(type_Central_InOut.InNeuronRatio);

data = {data1, data2, data3};
group = repelem(1:length(data), cellfun(@length, data));
combinedData = vertcat(data{:});
figure(3);set(gcf,'Color','w')
boxplot(combinedData, group, 'Labels', {'Optic', 'Feedforward', 'Central'},'Notch','on','Symbol','');

ylabel('InNeuronRatio');
ylim([0,15]);
set(gca,'Box','off','TickDir','out')

ax = gca;
boxes = findobj(ax, 'Tag', 'Box');
boxes(3).Color = [0.9290 0.6940 0.1250];   % Central
boxes(2).Color = [0 0.4470 0.7410];        % Feedforward
boxes(1).Color = [0.8500 0.3250 0.0980];   % Optic
medians = findobj(ax, 'Tag', 'Median');
medians(3).Color = [0.9290 0.6940 0.1250];
medians(2).Color = [0 0.4470 0.7410];
medians(1).Color = [0.8500 0.3250 0.0980];

%%%% Nonparametric test
groupNames = {'Optic', 'Feedforward', 'Central'};

dataPairs = {
    data1, data2;
    data1, data3;
    data2, data3
    };

groupPairs = {
    [1, 2];
    [1, 3];
    [2, 3]
    };

alpha = 0.05;
numComparisons = size(dataPairs, 1);
correctedAlpha = alpha / numComparisons;

fprintf('=== InNeuron Ratio Pairwise Comparison (Bonferroni corrected alpha = %.4f) ===\n\n', correctedAlpha);

for i = 1:numComparisons
    d1 = dataPairs{i, 1};
    d2 = dataPairs{i, 2};
    g1 = groupNames{groupPairs{i}(1)};
    g2 = groupNames{groupPairs{i}(2)};

    p = ranksum(d1, d2);

    p_dominance = mean(d1 > d2', 'all');
    p_dominance2 = mean(d1 < d2', 'all');

    fprintf('%s vs %s:\n', g1, g2);
    fprintf('  p-value = %.4f --> %s\n', p, ...
        ternary(p < correctedAlpha, 'Significant', 'Not significant'));
    fprintf('  P(%s > %s) = %.3f\n', g1, g2, p_dominance);
    fprintf('  P(%s < %s) = %.3f\n\n', g1, g2, p_dominance2);
end

%% 4. Out neuron number

data1 = cell2mat(type_Optic_InOut.OutNeuronNumber);
data2 = cell2mat(type_RightFFP_InOut.OutNeuronNumber);
data3 = cell2mat(type_Central_InOut.OutNeuronNumber);

data = {data1, data2, data3};
group = repelem(1:length(data), cellfun(@length, data));
combinedData = vertcat(data{:});
figure(4);set(gcf,'Color','w')
boxplot(combinedData, group, 'Labels', {'Optic', 'Feedforward', 'Central'},'Notch','on','Symbol','');

ylabel('OutNeuronNumber');
ylim([0, 1200]);
set(gca,'Box','off','TickDir','out')

ax = gca;
boxes = findobj(ax, 'Tag', 'Box');
boxes(3).Color = [0.9290 0.6940 0.1250];   % Central
boxes(2).Color = [0 0.4470 0.7410];        % Feedforward
boxes(1).Color = [0.8500 0.3250 0.0980];   % Optic
medians = findobj(ax, 'Tag', 'Median');
medians(3).Color = [0.9290 0.6940 0.1250];
medians(2).Color = [0 0.4470 0.7410];
medians(1).Color = [0.8500 0.3250 0.0980];

%%%% Nonparametric test
groupNames = {'Optic', 'Feedforward', 'Central'};

dataPairs = {
    data1, data2;
    data1, data3;
    data2, data3
    };

groupPairs = {
    [1, 2];
    [1, 3];
    [2, 3]
    };

alpha = 0.05;
numComparisons = size(dataPairs, 1);
correctedAlpha = alpha / numComparisons;

fprintf('=== OutNeuronNumber Pairwise Comparison (Bonferroni corrected alpha = %.4f) ===\n\n', correctedAlpha);

for i = 1:numComparisons
    d1 = dataPairs{i, 1};
    d2 = dataPairs{i, 2};
    g1 = groupNames{groupPairs{i}(1)};
    g2 = groupNames{groupPairs{i}(2)};

    p = ranksum(d1, d2);

    p_dominance = mean(d1 > d2', 'all');
    p_dominance2 = mean(d1 < d2', 'all');

    fprintf('%s vs %s:\n', g1, g2);
    fprintf('  p-value = %.4f --> %s\n', p, ...
        ternary(p < correctedAlpha, 'Significant', 'Not significant'));
    fprintf('  P(%s > %s) = %.3f\n', g1, g2, p_dominance);
    fprintf('  P(%s < %s) = %.3f\n\n', g1, g2, p_dominance2);
end

%% 5. Out neuron type number

data1 = cell2mat(type_Optic_InOut.OutNeuronTypeNumber);
data2 = cell2mat(type_RightFFP_InOut.OutNeuronTypeNumber);
data3 = cell2mat(type_Central_InOut.OutNeuronTypeNumber);

data = {data1, data2, data3};
group = repelem(1:length(data), cellfun(@length, data));
combinedData = vertcat(data{:});
figure(5);set(gcf,'Color','w')
boxplot(combinedData, group, 'Labels', {'Optic', 'Feedforward', 'Central'},'Notch','on','Symbol','');

ylabel('OutNeuronType');
ylim([0, 600]);
set(gca,'Box','off','TickDir','out')

ax = gca;
boxes = findobj(ax, 'Tag', 'Box');
boxes(3).Color = [0.9290 0.6940 0.1250];   % Central
boxes(2).Color = [0 0.4470 0.7410];        % Feedforward
boxes(1).Color = [0.8500 0.3250 0.0980];   % Optic
medians = findobj(ax, 'Tag', 'Median');
medians(3).Color = [0.9290 0.6940 0.1250];
medians(2).Color = [0 0.4470 0.7410];
medians(1).Color = [0.8500 0.3250 0.0980];

%%%% Nonparametric test
groupNames = {'Optic', 'Feedforward', 'Central'};

dataPairs = {
    data1, data2;
    data1, data3;
    data2, data3
    };

groupPairs = {
    [1, 2];
    [1, 3];
    [2, 3]
    };

alpha = 0.05;
numComparisons = size(dataPairs, 1);
correctedAlpha = alpha / numComparisons;

fprintf('=== OutNeuronType Pairwise Comparison (Bonferroni corrected alpha = %.4f) ===\n\n', correctedAlpha);

for i = 1:numComparisons
    d1 = dataPairs{i, 1};
    d2 = dataPairs{i, 2};
    g1 = groupNames{groupPairs{i}(1)};
    g2 = groupNames{groupPairs{i}(2)};

    p = ranksum(d1, d2);

    p_dominance = mean(d1 > d2', 'all');
    p_dominance2 = mean(d1 < d2', 'all');

    fprintf('%s vs %s:\n', g1, g2);
    fprintf('  p-value = %.4f --> %s\n', p, ...
        ternary(p < correctedAlpha, 'Significant', 'Not significant'));
    fprintf('  P(%s > %s) = %.3f\n', g1, g2, p_dominance);
    fprintf('  P(%s < %s) = %.3f\n\n', g1, g2, p_dominance2);
end

%% 6. Out neuron ratio

data1 = cell2mat(type_Optic_InOut.OutNeuronRatio);
data2 = cell2mat(type_RightFFP_InOut.OutNeuronRatio);
data3 = cell2mat(type_Central_InOut.OutNeuronRatio);

data = {data1, data2, data3};
group = repelem(1:length(data), cellfun(@length, data));
combinedData = vertcat(data{:});
figure(6);set(gcf,'Color','w')
boxplot(combinedData, group, 'Labels', {'Optic', 'Feedforward', 'Central'},'Notch','on','Symbol','');

ylabel('OutNeuronRatio');
ylim([0, 15]);
set(gca,'Box','off','TickDir','out')

ax = gca;
boxes = findobj(ax, 'Tag', 'Box');
boxes(3).Color = [0.9290 0.6940 0.1250];   % Central
boxes(2).Color = [0 0.4470 0.7410];        % Feedforward
boxes(1).Color = [0.8500 0.3250 0.0980];   % Optic
medians = findobj(ax, 'Tag', 'Median');
medians(3).Color = [0.9290 0.6940 0.1250];
medians(2).Color = [0 0.4470 0.7410];
medians(1).Color = [0.8500 0.3250 0.0980];

%%%% Nonparametric test
groupNames = {'Optic', 'Feedforward', 'Central'};

dataPairs = {
    data1, data2;
    data1, data3;
    data2, data3
    };

groupPairs = {
    [1, 2];
    [1, 3];
    [2, 3]
    };

alpha = 0.05;
numComparisons = size(dataPairs, 1);
correctedAlpha = alpha / numComparisons;

fprintf('=== OutNeuronRatio Pairwise Comparison (Bonferroni corrected alpha = %.4f) ===\n\n', correctedAlpha);

for i = 1:numComparisons
    d1 = dataPairs{i, 1};
    d2 = dataPairs{i, 2};
    g1 = groupNames{groupPairs{i}(1)};
    g2 = groupNames{groupPairs{i}(2)};

    p = ranksum(d1, d2);

    p_dominance = mean(d1 > d2', 'all');
    p_dominance2 = mean(d1 < d2', 'all');

    fprintf('%s vs %s:\n', g1, g2);
    fprintf('  p-value = %.4f --> %s\n', p, ...
        ternary(p < correctedAlpha, 'Significant', 'Not significant'));
    fprintf('  P(%s > %s) = %.3f\n', g1, g2, p_dominance);
    fprintf('  P(%s < %s) = %.3f\n\n', g1, g2, p_dominance2);
end

%% 7. Betweenness (unweighted)

data1 = cell2mat(type_Optic_InOut.betweenness_unweighted);
data2 = cell2mat(type_RightFFP_InOut.betweenness_unweighted);
data3 = cell2mat(type_Central_InOut.betweenness_unweighted);

data = {data1, data2, data3};
group = repelem(1:length(data), cellfun(@length, data));
combinedData = vertcat(data{:});
figure(7);set(gcf,'Color','w')
boxplot(combinedData, group, 'Labels', {'Optic', 'Feedforward', 'Central'},'Notch','on','Symbol','');

ylabel('betweenness\_unweighted');
ylim([0, 1.2e7]);
set(gca,'Box','off','TickDir','out')

ax = gca;
boxes = findobj(ax, 'Tag', 'Box');
boxes(3).Color = [0.9290 0.6940 0.1250];   % Central
boxes(2).Color = [0 0.4470 0.7410];        % Feedforward
boxes(1).Color = [0.8500 0.3250 0.0980];   % Optic
medians = findobj(ax, 'Tag', 'Median');
medians(3).Color = [0.9290 0.6940 0.1250];
medians(2).Color = [0 0.4470 0.7410];
medians(1).Color = [0.8500 0.3250 0.0980];

%%%% Nonparametric test
groupNames = {'Optic', 'Feedforward', 'Central'};

dataPairs = {
    data1, data2;
    data1, data3;
    data2, data3
    };

groupPairs = {
    [1, 2];
    [1, 3];
    [2, 3]
    };

alpha = 0.05;
numComparisons = size(dataPairs, 1);
correctedAlpha = alpha / numComparisons;

fprintf('=== Betweeness unweighted Comparison (Bonferroni corrected alpha = %.4f) ===\n\n', correctedAlpha);

for i = 1:numComparisons
    d1 = dataPairs{i, 1};
    d2 = dataPairs{i, 2};
    g1 = groupNames{groupPairs{i}(1)};
    g2 = groupNames{groupPairs{i}(2)};

    p = ranksum(d1, d2);

    p_dominance = mean(d1 > d2', 'all');
    p_dominance2 = mean(d1 < d2', 'all');

    fprintf('%s vs %s:\n', g1, g2);
    fprintf('  p-value = %.4f --> %s\n', p, ...
        ternary(p < correctedAlpha, 'Significant', 'Not significant'));
    fprintf('  P(%s > %s) = %.3f\n', g1, g2, p_dominance);
    fprintf('  P(%s < %s) = %.3f\n\n', g1, g2, p_dominance2);
end

%% 8. PageRank

data1 = cell2mat(type_Optic_InOut.pagerank);
data2 = cell2mat(type_RightFFP_InOut.pagerank);
data3 = cell2mat(type_Central_InOut.pagerank);

data = {data1, data2, data3};
group = repelem(1:length(data), cellfun(@length, data));
combinedData = vertcat(data{:});
figure(8);set(gcf,'Color','w')
boxplot(combinedData, group, 'Labels', {'Optic', 'Feedforward', 'Central'},'Notch','on','Symbol','');

ylabel('pagerank');
ylim([0, 10e-5]);
set(gca,'Box','off','TickDir','out')

ax = gca;
boxes = findobj(ax, 'Tag', 'Box');
boxes(3).Color = [0.9290 0.6940 0.1250];   % Central
boxes(2).Color = [0 0.4470 0.7410];        % Feedforward
boxes(1).Color = [0.8500 0.3250 0.0980];   % Optic
medians = findobj(ax, 'Tag', 'Median');
medians(3).Color = [0.9290 0.6940 0.1250];
medians(2).Color = [0 0.4470 0.7410];
medians(1).Color = [0.8500 0.3250 0.0980];

%%%% Nonparametric test
groupNames = {'Optic', 'Feedforward', 'Central'};

dataPairs = {
    data1, data2;
    data1, data3;
    data2, data3
    };

groupPairs = {
    [1, 2];
    [1, 3];
    [2, 3]
    };

alpha = 0.05;
numComparisons = size(dataPairs, 1);
correctedAlpha = alpha / numComparisons;

fprintf('=== Pagerank Comparison (Bonferroni corrected alpha = %.4f) ===\n\n', correctedAlpha);

for i = 1:numComparisons
    d1 = dataPairs{i, 1};
    d2 = dataPairs{i, 2};
    g1 = groupNames{groupPairs{i}(1)};
    g2 = groupNames{groupPairs{i}(2)};

    p = ranksum(d1, d2);

    p_dominance = mean(d1 > d2', 'all');
    p_dominance2 = mean(d1 < d2', 'all');

    fprintf('%s vs %s:\n', g1, g2);
    fprintf('  p-value = %.4f --> %s\n', p, ...
        ternary(p < correctedAlpha, 'Significant', 'Not significant'));
    fprintf('  P(%s > %s) = %.3f\n', g1, g2, p_dominance);
    fprintf('  P(%s < %s) = %.3f\n\n', g1, g2, p_dominance2);
end

%% Simple ternary helper
function out = ternary(cond, valTrue, valFalse)
if cond
    out = valTrue;
else
    out = valFalse;
end
end
