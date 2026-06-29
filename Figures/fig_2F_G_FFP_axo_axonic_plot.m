%% fig_2F_G_FFP_axo_axonic_plot
% Draws the FFP axo-axonic / optic-lobe proximity panels from the data saved by
% Data_Processing/s06_FFP_axo_axonic_compute.m (no recomputation here).
%
% Figure 1 : per-type axo-axonic input (% of total input) boxplot for ALL FFP
%            types, with the upper Tukey fence; also DEFINES the high
%            same-type-input outlier types.                           -> Fig 2F
% Figures 2-5 : optic-lobe proximity PROBABILITY P readout (per-type P bars;
%            P vs axo-axonic input at W > 0 and at W >= 5; per-type P boxplot vs
%            the 0.5 chance line). Reference analysis only -- NOT paper panels.
%            The paper optic-lobe panels 2G / 2H / S1D are produced instead by
%            Figures/fig_2G_H_FFP_axo_axonic_RF_plot.m, which computes the RF
%            distance (the OL synapse-centroid distance, um) binned by W.
%
% Proximity probability P is the per-neuron pooled (stratified) value computed
% in s06; distance = OL synapse centroid distance; W = directional same-type CB
% synapse count.
clear; close all; clc;

baseDir = fileparts(fileparts(mfilename('fullpath')));
load(fullfile(baseDir, 'Processed_Data', 'FFP_axo_axonic_outlier_data.mat'), 'D');

cGray   = [0.60 0.60 0.60];   % non-outlier types (light gray)
cDark   = [0.25 0.25 0.25];   % outlier types (dark gray)
cW1     = [0.55 0.70 0.90];   % W > 0  (light blue)
cW5     = [0.10 0.30 0.65];   % W >= 5 (dark blue)

%% ===== Figure 1 : axo-axonic input boxplot (all types) + upper fence =====
% Defines the outlier types via the upper Tukey fence (Q3 + 1.5*IQR).
vals   = D.axoaxonic_pct;
isOut  = D.is_outlier;
outIdx = find(isOut);

figure('Color','w','Name','Fig1 (2F) - axo-axonic input boxplot (fence definition)'); hold on;
boxplot(vals, 'Symbol','', 'Colors', cDark, 'Widths',0.5);
rng(0);
jit = (rand(numel(vals),1)-0.5)*0.28;
scatter(1+jit(~isOut), vals(~isOut), 16, cGray, 'filled', 'MarkerFaceAlpha',0.5);
scatter(1+jit( isOut), vals( isOut), 28, cDark, 'filled', 'MarkerFaceAlpha',0.9);
yline(D.upper_fence, '--', sprintf('upper fence = %.2f%%', D.upper_fence), 'Color',[0.4 0.4 0.4]);
for ii = 1:numel(outIdx)
    k = outIdx(ii);
    text(1.12, double(vals(k)), char(D.types(k)), 'FontSize',7, ...
        'Interpreter','none', 'VerticalAlignment','middle');
end
set(gca, 'XTick',1, 'XTickLabel',{'All FFP types'}, 'TickLabelInterpreter','none', ...
    'Box','off', 'TickDir','out', 'FontSize',12);
ylabel('Axo-axonic input from the same type (% of total input)');
title(sprintf('n = %d types, %d above upper fence', numel(vals), sum(isOut)));
hold off;

%% ===== Figure 2 (reference, not a paper panel) : per-type proximity probability bars =====
% Outlier types, ordered by P (W > 0) descending; mean across types as h-lines.
P      = D.out_P;                 % nOut x [W>0, W>=5]
types  = cellstr(D.out_types);
[~, ord] = sort(P(:,1), 'descend', 'MissingPlacement','last');
keepB  = any(~isnan(P(ord,:)), 2);
ord    = ord(keepB);
Pb     = P(ord, :);
tb     = types(ord);
nB     = numel(ord);
xs     = 1:nB;

Pmean1 = mean(P(~isnan(P(:,1)),1));
Pmean5 = mean(P(~isnan(P(:,2)),2));

figure('Color','w','Name','Fig2 (reference) - per-type proximity probability bars', ...
    'Position',[120 120 max(700,nB*40) 480]); hold on;
hb = bar(xs, Pb, 'grouped', 'EdgeColor','none');
hb(1).FaceColor = cW1; hb(2).FaceColor = cW5;
yline(0.5, '--', 'chance (0.5)', 'Color',[0.4 0.4 0.4]);
yline(Pmean1, '-', sprintf('mean W>0 = %.3f', Pmean1),  'Color',cW1, 'LineWidth',1.3);
yline(Pmean5, '-', sprintf('mean W\\geq5 = %.3f', Pmean5), 'Color',cW5, 'LineWidth',1.3);
set(gca, 'XTick',xs, 'XTickLabel',tb, 'XTickLabelRotation',90, 'TickLabelInterpreter','none', ...
    'Box','off', 'TickDir','out', 'FontSize',11);
xlim([0 nB+1]); ylim([0 1]);
ylabel('Proximity probability P');
legend({'W>0','W\geq5'}, 'Location','northeast', 'Box','off');
title('Per-type proximity probability (outlier types)');
hold off;

%% ===== Figure 3 (reference, not a paper panel) : P (W > 0) vs axo-axonic input =====
x3 = D.out_axoaxonic_pct;
y3 = D.out_P(:, 1);                % W > 0
keep3 = ~isnan(x3) & ~isnan(y3);
x3 = x3(keep3); y3 = y3(keep3);
lab3 = cellstr(D.out_types(keep3));

figure('Color','w','Name','Fig3 (reference) - P (W>0) vs axo-axonic input', ...
    'Position',[160 160 620 560]); hold on;
scatter(x3, y3, 50, cW5, 'filled', 'MarkerFaceAlpha',0.75);
text(x3, y3, lab3, 'FontSize',8, 'VerticalAlignment','bottom', ...
     'HorizontalAlignment','left', 'Interpreter','none', 'Color',[0.3 0.3 0.3]);
yline(0.5, '--', 'chance (0.5)', 'Color',[0.4 0.4 0.4]);
if numel(x3) >= 3
    pf = polyfit(x3, y3, 1); xx = linspace(min(x3), max(x3), 100);
    plot(xx, polyval(pf, xx), '-', 'Color',[0.85 0.33 0.10], 'LineWidth',1.3);
    r_s = corr(x3, y3, 'type','Spearman'); r_p = corr(x3, y3, 'type','Pearson');
    title(sprintf('P (W>0) vs axo-axonic input  |  Spearman \\rho=%.2f, Pearson r=%.2f', r_s, r_p));
else
    title('P (W>0) vs axo-axonic input');
end
xlabel('Axo-axonic input from the same type (% of total input)');
ylabel('Proximity probability P (W > 0)');
ylim([0 1]); set(gca, 'Box','off', 'TickDir','out', 'FontSize',12);
hold off;

%% ===== Figure 4 (reference, not a paper panel) : P (W >= 5) vs axo-axonic input =====
x4 = D.out_axoaxonic_pct;
y4 = D.out_P(:, 2);                % W >= 5
keep4 = ~isnan(x4) & ~isnan(y4);
x4 = x4(keep4); y4 = y4(keep4);
lab4 = cellstr(D.out_types(keep4));

figure('Color','w','Name','Fig4 (reference) - P (W>=5) vs axo-axonic input', ...
    'Position',[200 200 620 560]); hold on;
scatter(x4, y4, 50, cW5, 'filled', 'MarkerFaceAlpha',0.75);
text(x4, y4, lab4, 'FontSize',8, 'VerticalAlignment','bottom', ...
     'HorizontalAlignment','left', 'Interpreter','none', 'Color',[0.3 0.3 0.3]);
yline(0.5, '--', 'chance (0.5)', 'Color',[0.4 0.4 0.4]);
if numel(x4) >= 3
    pf = polyfit(x4, y4, 1); xx = linspace(min(x4), max(x4), 100);
    plot(xx, polyval(pf, xx), '-', 'Color',[0.85 0.33 0.10], 'LineWidth',1.3);
    r_s = corr(x4, y4, 'type','Spearman'); r_p = corr(x4, y4, 'type','Pearson');
    title(sprintf('P (W\\geq5) vs axo-axonic input  |  Spearman \\rho=%.2f, Pearson r=%.2f', r_s, r_p));
else
    title('P (W\geq5) vs axo-axonic input');
end
xlabel('Axo-axonic input from the same type (% of total input)');
ylabel('Proximity probability P (W \geq 5)');
ylim([0 1]); set(gca, 'Box','off', 'TickDir','out', 'FontSize',12);
hold off;

%% ===== Figure 5 (reference, not a paper panel) : per-type P boxplot vs 0.5 =====
a1 = D.out_P(~isnan(D.out_P(:,1)), 1);    % W > 0
a5 = D.out_P(~isnan(D.out_P(:,2)), 2);    % W >= 5
data5 = [a1; a5];
grp5  = [ones(numel(a1),1); 2*ones(numel(a5),1)];

p1 = NaN; p5 = NaN;                        % one-sided Wilcoxon signed-rank vs 0.5
if numel(a1) >= 2, p1 = signrank(a1, 0.5, 'tail','right'); end
if numel(a5) >= 2, p5 = signrank(a5, 0.5, 'tail','right'); end

figure('Color','w','Name','Fig5 (reference) - per-type P boxplot', ...
    'Position',[240 240 460 520]); hold on;
boxplot(data5, grp5, 'Labels', {'W>0','W\geq5'}, 'Symbol','', 'Widths',0.5);
rng(0);
xj1 = 1 + (rand(numel(a1),1)-0.5)*0.30;
xj5 = 2 + (rand(numel(a5),1)-0.5)*0.30;
scatter(xj1, a1, 26, cW1, 'filled', 'MarkerFaceAlpha',0.6);
scatter(xj5, a5, 26, cW5, 'filled', 'MarkerFaceAlpha',0.6);
yline(0.5, '--', 'chance (0.5)', 'Color',[0.4 0.4 0.4]);
ylim([0 1]); ylabel('Proximity probability P');
title(sprintf('signrank vs 0.5:  W>0 p=%.2g,  W\\geq5 p=%.2g', p1, p5));
set(gca, 'Box','off', 'TickDir','out', 'FontSize',12);
text(1, 0.04, sprintf('med=%.3f\nn=%d', median(a1), numel(a1)), 'HorizontalAlignment','center', 'FontSize',8);
text(2, 0.04, sprintf('med=%.3f\nn=%d', median(a5), numel(a5)), 'HorizontalAlignment','center', 'FontSize',8);
hold off;

fprintf('Drawn: Fig1 (2F) / Fig2-5 (proximity-probability reference, not paper panels).\n');
