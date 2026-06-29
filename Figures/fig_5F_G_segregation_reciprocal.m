%% 1. Load data and initialize
clear all; close all; clc

% Project root resolved from this script's location (Figures/ -> project root),
% so the analysis runs after download without editing any paths.
baseDir = fileparts(fileparts(mfilename('fullpath')));

% Per-type input/output reciprocity (weighted Jaccard) + example neurons (Want),
% saved by Data_Processing/s14_BDP_reciprocal.m
load(fullfile(baseDir, 'Processed_Data', 'BDP_reciprocity.mat'))

% Attach the per-type mean segregation index from the s13 CSVs
% (Data_Processing/s13_segregation_index.ipynb).
Type_FFP    = attach_mean_segregation(Type_FFP,    fullfile(baseDir, 'Processed_Data', 'right_FFP_segregation_index.csv'));
Type_FBP    = attach_mean_segregation(Type_FBP,    fullfile(baseDir, 'Processed_Data', 'right_FBP_segregation_index.csv'));
Type_BDP    = attach_mean_segregation(Type_BDP,    fullfile(baseDir, 'Processed_Data', 'right_BDP_segregation_index.csv'));
Type_Others = attach_mean_segregation(Type_Others, fullfile(baseDir, 'Processed_Data', 'right_Others_segregation_index.csv'));

% Reciprocity per group
FFP_Reci    = Type_FFP.weighted_Jaccard;
FBP_Reci    = Type_FBP.weighted_Jaccard;
BDP_Reci    = Type_BDP.weighted_Jaccard;
Others_Reci = Type_Others.weighted_Jaccard;

% Segregation index per group
FFP_segIdx    = Type_FFP.mean_segregation_index;
FBP_segIdx    = Type_FBP.mean_segregation_index;
BDP_segIdx    = Type_BDP.mean_segregation_index;
Others_segIdx = Type_Others.mean_segregation_index;

% Example BDP neurons (per neuron reciprocity)
LC9_Reci  = Want(1).weighted_Jaccard;
LT43_Reci = Want(2).weighted_Jaccard;
LT52_Reci = Want(3).weighted_Jaccard;

%% 2. Figure 1 (panel 5F): segregation index vs reciprocity
figure(1); set(gcf,'Color','w'); hold on;

scatter(FFP_segIdx, FFP_Reci, 'filled', ...
    'MarkerFaceAlpha',0.5, 'MarkerFaceColor',[0 0.4470 0.7410], ...
    'DisplayName', sprintf('Feedforward (n=%d)', numel(FFP_Reci)));
scatter(FBP_segIdx, FBP_Reci, 'filled', ...
    'MarkerFaceAlpha',0.5, 'MarkerFaceColor',[0.4660 0.6740 0.1880], ...
    'DisplayName', sprintf('Feedback (n=%d)', numel(FBP_Reci)));
scatter(BDP_segIdx, BDP_Reci, 'filled', ...
    'MarkerFaceAlpha',0.5, 'MarkerFaceColor','#5cb599', ...
    'DisplayName', sprintf('Bidirectional (n=%d)', numel(BDP_Reci)));
scatter(Others_segIdx, Others_Reci, 'filled', ...
    'MarkerFaceAlpha',0.5, 'MarkerFaceColor',[0.8 0.8 0.8], ...
    'DisplayName', sprintf('Others (n=%d)', numel(Others_Reci)));

% Linear fit + prediction interval over all groups
xAll = [FFP_segIdx; FBP_segIdx; BDP_segIdx; Others_segIdx];
yAll = [FFP_Reci;   FBP_Reci;   BDP_Reci;   Others_Reci];
mask = isfinite(xAll) & isfinite(yAll);
xAll = xAll(mask); yAll = yAll(mask);

mdl  = fitlm(xAll, yAll);
xfit = linspace(min(xAll), max(xAll), 200)';
[yfit, ~] = predict(mdl, xfit);
[~, yPI]  = predict(mdl, xfit, 'Prediction', 'observation');

plot(xfit, yfit, 'k-', 'LineWidth', 2, 'DisplayName', 'Linear fit (all)');
patch([xfit; flipud(xfit)], [yPI(:,1); flipud(yPI(:,2))], [0 0 0], ...
    'FaceAlpha', 0.06, 'EdgeColor','none', 'DisplayName', '95% PI (new obs)');

% Fit statistics annotation
b0 = mdl.Coefficients.Estimate(1);
b1 = mdl.Coefficients.Estimate(2);
R2 = mdl.Rsquared.Ordinary;
[rp, p_r] = corr(xAll, yAll, 'Type','Pearson', 'Rows','complete');
p_slope   = mdl.Coefficients.pValue(2);
txt = sprintf('y = %.3f + %.3f x   R^2=%.3f   r=%.3f (p_r=%.3g)   p_{slope}=%.3g', ...
    b0, b1, R2, rp, p_r, p_slope);
text(min(xAll) + 0.02*range(xAll), max(yAll) - 0.05*range(yAll), txt, ...
    'FontSize',10, 'BackgroundColor','w', 'Margin',2);

% Highlight the example BDP types
LC9_idx  = strcmp(Type_BDP.type, 'LC9');
LT43_idx = strcmp(Type_BDP.type, 'LT43');
LT52_idx = strcmp(Type_BDP.type, 'LT52');
scatter(Type_BDP.mean_segregation_index(LC9_idx),  Type_BDP.weighted_Jaccard(LC9_idx),  'filled', 'MarkerFaceColor','#dc421c');
scatter(Type_BDP.mean_segregation_index(LT43_idx), Type_BDP.weighted_Jaccard(LT43_idx), 'filled', 'MarkerFaceColor','#b26720');
scatter(Type_BDP.mean_segregation_index(LT52_idx), Type_BDP.weighted_Jaccard(LT52_idx), 'filled', 'MarkerFaceColor','#772e86');

hold off;
xlabel('Segregation index'); ylabel('Reciprocal %');
set(gca,'TickDir','out'); box off;
legend('Location','bestoutside');
xlim([0 1]); ylim([0 0.4]); axis square

%% 3. Figure 2 (panel 5G): reciprocity box plot
data = {LC9_Reci, LT43_Reci, LT52_Reci, BDP_Reci, FFP_Reci, FBP_Reci};
group = repelem(1:length(data), cellfun(@length, data));
combinedData = vertcat(data{:});

figure(2); set(gcf,'Color','w')
boxplot(combinedData, group, 'Labels', {'LC9','LT43','LT52','BDP','FFP','FBP'}, 'Notch','off','Symbol','');
title('Reci')
set(gca,'Box','off','TickDir','out')
ylim([0 0.3])

% Pairwise comparison (BDP / FFP / FBP)
groupNames = {'BDP', 'FFP', 'FBP'};
dataPairs = {
    BDP_Reci, FFP_Reci;
    BDP_Reci, FBP_Reci;
    FFP_Reci, FBP_Reci;
    };
groupPairs = {[1, 2]; [1, 3]; [2, 3]};

alpha = 0.05;
numComparisons = size(dataPairs, 1);
correctedAlpha = alpha / numComparisons;

fprintf('=== Reciprocity Pairwise Comparison (Bonferroni corrected alpha = %.4f) ===\n\n', correctedAlpha);
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
function T = attach_mean_segregation(T, csvPath)
% Mean segregation index per type, averaged over each type's member root_ids.
opts = detectImportOptions(csvPath);
opts = setvartype(opts, 'root_id', 'int64');
S = readtable(csvPath, opts);

T.mean_segregation_index = nan(height(T), 1);
for i = 1:height(T)
    m = ismember(S.root_id, T.root_id{i});
    T.mean_segregation_index(i) = mean(S.segregation_index(m), 'omitnan');
end
end

function out = ternary(cond, valTrue, valFalse)
if cond
    out = valTrue;
else
    out = valFalse;
end
end
