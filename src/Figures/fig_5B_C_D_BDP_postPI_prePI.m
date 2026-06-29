%% 1. Load data and initialize
clear all; close all; clc

baseDir = fileparts(fileparts(mfilename('fullpath')));  % resolved from this script's location (portable; no absolute paths)

SeeType = 1;   % 1: aggregate per cell type, 0: per neuron

% right_neurons_thr0.mat is produced by fig_1D_E_postPI_prePI_right_neurons.m and
% contains RightFFP_*, RightFBP_*, and the RightBDP_real_* (bidirectional) tables.
load(fullfile(baseDir, 'Processed_Data', 'right_neurons_thr0.mat'))

%% 2. Select PostPI / PrePI (per type or per neuron)
if SeeType
    FFP_postPI = RightFFP_by_type.Mean_Right_PostPI;
    FBP_postPI = RightFBP_by_type.Mean_Right_PostPI;
    BDP_postPI  = RightBDP_real_type.Mean_Right_PostPI;

    FFP_prePI = RightFFP_by_type.Mean_Right_PrePI;
    FBP_prePI = RightFBP_by_type.Mean_Right_PrePI;
    BDP_prePI  = RightBDP_real_type.Mean_Right_PrePI;
else
    FFP_postPI = RightFFP_NPIs.Right_PostPI;
    FBP_postPI = RightFBP_NPIs.Right_PostPI;
    BDP_postPI  = RightBDP_real_NPIs.Right_PostPI;

    FFP_prePI = RightFFP_NPIs.Right_PrePI;
    FBP_prePI = RightFBP_NPIs.Right_PrePI;
    BDP_prePI  = RightBDP_real_NPIs.Right_PrePI;
end

% Example bidirectional neurons (always per neuron)
LC9_idx  = strcmp(RightBDP_real_NPIs.type, 'LC9');
LT43_idx = strcmp(RightBDP_real_NPIs.type, 'LT43');
LT52_idx = strcmp(RightBDP_real_NPIs.type, 'LT52');

LC9_postPI  = RightBDP_real_NPIs.Right_PostPI(LC9_idx);
LC9_prePI   = RightBDP_real_NPIs.Right_PrePI(LC9_idx);
LT43_postPI = RightBDP_real_NPIs.Right_PostPI(LT43_idx);
LT43_prePI  = RightBDP_real_NPIs.Right_PrePI(LT43_idx);
LT52_postPI = RightBDP_real_NPIs.Right_PostPI(LT52_idx);
LT52_prePI  = RightBDP_real_NPIs.Right_PrePI(LT52_idx);

%% 3. Figure 1 (panel 5B): bias-direction scatter
figure(1); set(gcf,'Color','w'); hold on;

% Diamond guide lines (|Post-Pre| = 2 boundary)
plot(linspace(-2,0,10), linspace(0,2,10),  'Color','#EAEBEB','LineWidth',1);
plot(linspace(-2,0,10), linspace(0,-2,10), 'Color','#EAEBEB','LineWidth',1);
plot(linspace(0,2,10),  linspace(2,0,10),  'Color','#EAEBEB','LineWidth',1);
plot(linspace(0,2,10),  linspace(-2,0,10), 'Color','#EAEBEB','LineWidth',1);

% Group means with SD error bars
groups = {RightFFP_by_type, RightFBP_by_type, RightBDP_real_type};
colors = [0.0000, 0.4470, 0.7410;   % Blue  (FFP)
          0.4660, 0.6740, 0.1880;   % Green (FBP)
          0.6500, 0.6500, 0.6500];  % Gray  (BDP)

for i = 1:length(groups)
    data = groups{i};
    x_raw = data.Mean_Right_PostPI + data.Mean_Right_PrePI;
    y_raw = data.Mean_Right_PostPI - data.Mean_Right_PrePI;

    mX = mean(x_raw);
    mY = mean(y_raw);
    eX = std(x_raw);   % standard deviation (use /sqrt(n) for SEM)
    eY = std(y_raw);

    scatter(mX, mY, 100, 'filled', 'square', 'MarkerFaceColor', colors(i,:));
    errorbar(mX, mY, eY, eY, eX, eX, ...
        'Color', colors(i,:), 'LineStyle', 'none', 'LineWidth', 1.2, 'CapSize', 10);
end

% Highlight the example bidirectional neuron types
x_raw_LC9 = RightBDP_real_type.Mean_Right_PostPI(strcmp(RightBDP_real_type.type,'LC9')) + ...
            RightBDP_real_type.Mean_Right_PrePI(strcmp(RightBDP_real_type.type,'LC9'));
y_raw_LC9 = RightBDP_real_type.Mean_Right_PostPI(strcmp(RightBDP_real_type.type,'LC9')) - ...
            RightBDP_real_type.Mean_Right_PrePI(strcmp(RightBDP_real_type.type,'LC9'));
scatter(x_raw_LC9, y_raw_LC9, 50, 'filled', 'o', 'MarkerFaceColor', '#dc421c');

x_raw_LT43 = RightBDP_real_type.Mean_Right_PostPI(strcmp(RightBDP_real_type.type,'LT43')) + ...
             RightBDP_real_type.Mean_Right_PrePI(strcmp(RightBDP_real_type.type,'LT43'));
y_raw_LT43 = RightBDP_real_type.Mean_Right_PostPI(strcmp(RightBDP_real_type.type,'LT43')) - ...
             RightBDP_real_type.Mean_Right_PrePI(strcmp(RightBDP_real_type.type,'LT43'));
scatter(x_raw_LT43, y_raw_LT43, 50, 'filled', 'o', 'MarkerFaceColor', '#5cb599');

x_raw_LT52 = RightBDP_real_type.Mean_Right_PostPI(strcmp(RightBDP_real_type.type,'LT52')) + ...
             RightBDP_real_type.Mean_Right_PrePI(strcmp(RightBDP_real_type.type,'LT52'));
y_raw_LT52 = RightBDP_real_type.Mean_Right_PostPI(strcmp(RightBDP_real_type.type,'LT52')) - ...
             RightBDP_real_type.Mean_Right_PrePI(strcmp(RightBDP_real_type.type,'LT52'));
scatter(x_raw_LT52, y_raw_LT52, 50, 'filled', 'o', 'MarkerFaceColor', '#b5287a');

grid on;
xlabel('Post + Pre (B)');
ylabel('Post - Pre (D)');
set(gca,'TickDir','out')
axis equal
xlim([-2 2]); ylim([-2 2]);

%% 4. Figure 2 (panel 5C): PostPI box plot
data = {LC9_postPI, LT43_postPI, LT52_postPI, BDP_postPI, FFP_postPI, FBP_postPI};
group = repelem(1:length(data), cellfun(@length, data));
combinedData = vertcat(data{:});

figure(2); set(gcf,'Color','w')
boxplot(combinedData, group, 'Labels', {'LC9','LT43','LT52','BDP','FFP','FBP'}, 'Notch','off','Symbol','');
title('postPI')
set(gca,'Box','off','TickDir','out')
ylim([-1 1])

% Pairwise comparison (BDP / FFP / FBP)
groupNames = {'BDP', 'FFP', 'FBP'};
dataPairs = {
    BDP_postPI, FFP_postPI;
    BDP_postPI, FBP_postPI;
    FFP_postPI, FBP_postPI;
    };
groupPairs = {[1, 2]; [1, 3]; [2, 3]};

alpha = 0.05;
numComparisons = size(dataPairs, 1);
correctedAlpha = alpha / numComparisons;

fprintf('=== PostPI Pairwise Comparison (Bonferroni corrected alpha = %.4f) ===\n\n', correctedAlpha);
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

%% 5. Figure 3 (panel 5D): PrePI box plot
data = {LC9_prePI, LT43_prePI, LT52_prePI, BDP_prePI, FFP_prePI, FBP_prePI};
group = repelem(1:length(data), cellfun(@length, data));
combinedData = vertcat(data{:});

figure(3); set(gcf,'Color','w')
boxplot(combinedData, group, 'Labels', {'LC9','LT43','LT52','BDP','FFP','FBP'}, 'Notch','off','Symbol','');
title('prePI')
set(gca,'Box','off','TickDir','out')
ylim([-1 1])

% Pairwise comparison (BDP / FFP / FBP)
groupNames = {'BDP', 'FFP', 'FBP'};
dataPairs = {
    BDP_prePI, FFP_prePI;
    BDP_prePI, FBP_prePI;
    FFP_prePI, FBP_prePI;
    };
groupPairs = {[1, 2]; [1, 3]; [2, 3]};

alpha = 0.05;
numComparisons = size(dataPairs, 1);
correctedAlpha = alpha / numComparisons;

fprintf('=== PrePI Pairwise Comparison (Bonferroni corrected alpha = %.4f) ===\n\n', correctedAlpha);
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
function out = ternary(cond, valTrue, valFalse)
if cond
    out = valTrue;
else
    out = valFalse;
end
end
