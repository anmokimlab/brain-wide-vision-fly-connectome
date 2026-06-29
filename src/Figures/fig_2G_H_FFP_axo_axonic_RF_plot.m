%% fig_2G_H_FFP_axo_axonic_RF_plot
% Draws the FFP optic-lobe proximity panels from the data saved by
% Data_Processing/s07_FFP_axo_axonic_RF.m (no recomputation here). The CSV holds,
% for each high same-type-input (axo-axonic) outlier type, the mean RF distance
% (the OL synapse-centroid distance, um) binned by central-brain same-type
% connection strength W, plus the connected(W>0) / unconnected(W==0) split.
%
% Figure 1 (Fig 2G)  : per-type 2D scatter, x = connected (W>0) RF distance,
%                      y = unconnected (W==0) RF distance. Points above the y=x
%                      line have connected partners closer than unconnected ones.
% Figure 2 (Fig 2H)  : per-type mean RF distance by W bin (boxplot) with
%                      Kruskal-Wallis omnibus and adjacent-bin ranksum stars.
% Figure 3 (Fig S1D) : connected (W>0) vs unconnected (W=0) grouped bar chart,
%                      sorted by D_unconnected / D_connected; error bar = std.
%
% RF distance = OL synapse-centroid distance (um): each neuron's OL synapses are
% reduced to a single centroid, and the pair distance is the Euclidean distance
% between the two centroids. W = directional same-type CB synapse count (pre -> partner).
clear; close all; clc

baseDir = fileparts(fileparts(mfilename('fullpath')));

binLabels = {'W0','W1to4','W5plus'};
binDisp   = {'W=0','0<W<5','5\leqW'};
nBins     = numel(binLabels);

%% ===== Load data =====
T = readtable(fullfile(baseDir, 'Processed_Data', 'FFP_axo_axonic_RF_proximity.csv'), 'TextType','string');
nT = height(T);
fprintf('Loaded %d high outlier types: %s\n', nT, strjoin(cellstr(T.type), ', '));

D = nan(nT, nBins);   % per-bin mean distance
S = nan(nT, nBins);   % SEM
for b = 1:nBins
    D(:,b) = T.(['Distance_um_' binLabels{b}]);
    S(:,b) = T.(['SEM_um_'      binLabels{b}]);
end

% Per-bin colors (weak -> strong : light -> dark blue)
cols = [0.78 0.84 0.92;
        0.50 0.66 0.86;
        0.22 0.45 0.74;
        0.06 0.24 0.52];
xpos = 1:nBins;

%% ===== Figure 1 (Fig 2G) : connected distance(x) vs unconnected distance(y) =====
% Per type: x = connected (W>0) RF distance, y = unconnected (W=0) RF distance.
% Above y=x => unconnected is farther (= connected is closer).
all_types = cellstr(T.type);
Dc_f1 = T.Distance_um_connected;     % x : connected
Du_f1 = T.Distance_um_unconnected;   % y : unconnected
v_f1  = ~isnan(Dc_f1) & ~isnan(Du_f1);
Dc_f1 = Dc_f1(v_f1); Du_f1 = Du_f1(v_f1);
types_f1 = all_types(v_f1);

% Per-bin type mean +/- SEM (kept for Figure 2 / console summary)
mBin = nan(1,nBins); seBin = nan(1,nBins); nValidBin = zeros(1,nBins);
for b = 1:nBins
    v = ~isnan(D(:,b));
    nValidBin(b) = sum(v);
    mBin(b)  = mean(D(v,b));
    seBin(b) = std(D(v,b)) / sqrt(sum(v));
end

figure('Color','w','Name','Fig1 (2G) - RF distance: connected vs unconnected (2D scatter)','Position',[120 120 560 560]);
hold on
% y=x reference line
axlim = [0, max([Dc_f1; Du_f1])*1.08];
plot(axlim, axlim, '--', 'Color',[0.6 0.6 0.6], 'LineWidth',1);
% per-type points
scatter(Dc_f1, Du_f1, 55, cols(3,:), 'filled', 'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.8);
% type labels
for g = 1:numel(types_f1)
    text(Dc_f1(g)+axlim(2)*0.01, Du_f1(g), types_f1{g}, 'FontSize',7, 'Interpreter','none', 'Color',[0.4 0.4 0.4]);
end
axis equal
xlim(axlim); ylim(axlim);
set(gca, 'TickDir','out');
xlabel('Connected (W>0) RF distance (\mum)');
ylabel('Unconnected (W=0) RF distance (\mum)');
title('RF distance: connected vs unconnected per type  (above y=x : connected closer)');
box off

%% ===== Figure 2 (Fig 2H) : per-bin boxplot + stats =====
% long-form for boxplot
data = []; grp = [];
for b = 1:nBins
    v = ~isnan(D(:,b));
    data = [data; D(v,b)];
    grp  = [grp;  b*ones(sum(v),1)];
end
figure('Color','w','Name','Fig2 (2H) - RF distance by W bin (boxplot)','Position',[150 150 560 560]);
boxplot(data, grp, 'Labels', binDisp, 'Symbol','', 'Widths',0.6); hold on
rng(0);
for b = 1:nBins
    v = ~isnan(D(:,b));
    xj = b + (rand(sum(v),1)-0.5)*0.30;
    scatter(xj, D(v,b), 30, cols(b,:), 'filled', 'MarkerFaceAlpha',0.7, 'MarkerEdgeColor',[0.3 0.3 0.3]);
end
ylabel('Mean RF distance (\mum)');
xlabel('Connection strength W (central brain)');

% Stats (all unpaired): omnibus = Kruskal-Wallis over all types with a value per bin
p_kw = NaN;
if numel(unique(grp)) >= 3
    try, p_kw = kruskalwallis(data, grp, 'off'); catch, end
end
title(sprintf('RF distance by W bin  (Kruskal-Wallis p=%.2g)', p_kw));

% Adjacent-bin significance: ranksum (unpaired, tail=left=closer) + Bonferroni stars
nAdj   = nBins - 1;                 % adjacent comparisons (Bonferroni factor)
ymax   = max(data);
baseY  = ymax*1.05;
stepY  = ymax*0.09;
tick   = ymax*0.02;
for k = 1:nAdj
    a = k; b = k+1;                 % adjacent pair (a = previous, b = next)
    va = ~isnan(D(:,a)); vb = ~isnan(D(:,b));
    if sum(va) < 2 || sum(vb) < 2, continue; end
    p_raw = ranksum(D(vb,b), D(va,a), 'tail','left');
    p_adj = min(p_raw*nAdj, 1);     % Bonferroni
    yk = baseY + (k-1)*stepY;
    plot([a a b b], [yk-tick yk yk yk-tick], 'k-', 'LineWidth',1);
    text((a+b)/2, yk+tick*0.4, pstars(p_adj), 'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', 'FontSize',11, 'FontWeight','bold');
end
set(gca,'TickDir','out'); box off;
ylim([0, baseY + nAdj*stepY]);

%% ===== Figure 3 (Fig S1D) : connected(W>0) vs unconnected(W=0) grouped bar =====
% Per-type RF distance split by connectivity, sorted by D_unconnected/D_connected
% descending. Error bar = std (compute saved SEM = std/sqrt(n) -> std = SEM*sqrt(n)).
Dc_cu  = T.Distance_um_connected;    Du_cu  = T.Distance_um_unconnected;
STDc   = T.SEM_um_connected   .* sqrt(T.n_neurons_connected);
STDu   = T.SEM_um_unconnected .* sqrt(T.n_neurons_unconnected);
types_cu = cellstr(T.type);

valid_cu = ~isnan(Dc_cu) & ~isnan(Du_cu) & Dc_cu > 0;
Dc_cu = Dc_cu(valid_cu); Du_cu = Du_cu(valid_cu);
STDc  = STDc(valid_cu);  STDu = STDu(valid_cu);
types_cu = types_cu(valid_cu);

ratio_cu = Du_cu ./ Dc_cu;                      % unconnected / connected
[ratio_cu, ord_cu] = sort(ratio_cu, 'descend'); % largest first
Dc_cu = Dc_cu(ord_cu); Du_cu = Du_cu(ord_cu);
STDc  = STDc(ord_cu);  STDu = STDu(ord_cu);
types_cu = types_cu(ord_cu);
nT_cu = numel(types_cu);
fprintf('Types with both connected and unconnected %d (ratio descending)\n', nT_cu);

colU = [0.70 0.70 0.70];   % unconnected (gray)
colC = [0.10 0.30 0.65];   % connected   (blue)
y_cu  = [Dc_cu, Du_cu];     % nT_cu x 2  (col1 connected, col2 unconnected)
err_cu = [STDc, STDu];      % error bar = std

figure('Color','w','Name','Fig3 (S1D) - RF distance connected vs unconnected', ...
    'Position',[100 100 max(900, nT_cu*42) 520]);
hb = bar(y_cu, 'grouped', 'EdgeColor','none'); hold on
hb(1).FaceColor = colC;
hb(2).FaceColor = colU;

ngroups = nT_cu; nbars = 2;
groupwidth = min(0.8, nbars/(nbars+1.5));
for k = 1:nbars
    xk = (1:ngroups) - groupwidth/2 + (2*k-1)*groupwidth/(2*nbars);
    errorbar(xk, y_cu(:,k), err_cu(:,k), 'k', 'linestyle','none', 'LineWidth',0.7, 'CapSize',3);
end

% x tick labels = "type (ratio)"  e.g. MeTu2a (2.51)   (ratio = D_unconnected/D_connected)
xtl_cu = arrayfun(@(g) sprintf('%s (%.2f)', types_cu{g}, ratio_cu(g)), 1:nT_cu, 'UniformOutput', false);
set(gca, 'XTick',1:nT_cu, 'XTickLabel',xtl_cu, 'XTickLabelRotation',90, ...
    'TickLabelInterpreter','none', 'TickDir','out');
xlim([0.5 nT_cu+0.5]);
ylabel('Mean RF distance (\mum)');
xlabel('FF type (ratio = D_{unconnected}/D_{connected}, sorted descending)');
title('RF distance: unconnected (W=0) vs connected (W>0)  | error bar = std');
legend({'connected (W>0)','unconnected (W=0)'}, 'Location','northeast', 'Box','off');
box off

%% ===== Console summary : bin trend =====
fprintf('\n========== Summary (high outlier types, per-bin mean RF distance) ==========\n');
for b = 1:nBins
    fprintf('  [%-8s] %3d types | mean distance = %.3f um | total pairs = %d\n', ...
        binLabels{b}, nValidBin(b), mBin(b), sum(T.(['n_pairs_' binLabels{b}])));
end
fprintf('  Kruskal-Wallis p(bin difference, unpaired) = %.3g\n', p_kw);
% pairwise ranksum vs W0 (unpaired, tail=left : smaller distance = closer)
v1 = ~isnan(D(:,1));
for b = 2:nBins
    vb = ~isnan(D(:,b));
    if sum(vb) >= 2 && sum(v1) >= 2
        p = ranksum(D(vb,b), D(v1,1), 'tail','left');
        fprintf('  %s vs W0 : median diff %+.3f um | ranksum(closer) p = %.3g (n=%d vs W0 n=%d)\n', ...
            binLabels{b}, median(D(vb,b))-median(D(v1,1)), p, sum(vb), sum(v1));
    end
end
% adjacent-bin pairwise (vs previous bin): is the stronger bin closer (distance decrease, tail=left)
% all unpaired ranksum + Bonferroni (same basis as the Figure 2 stars)
nAdj = nBins - 1;
fprintf('  --- adjacent-bin pairwise (vs previous, closer = distance decrease, Bonferroni) ---\n');
for b = 2:nBins
    vb = ~isnan(D(:,b)); va = ~isnan(D(:,b-1));
    if sum(vb) >= 2 && sum(va) >= 2
        % (a) unpaired ranksum : all types in each bin
        p_raw = ranksum(D(vb,b), D(va,b-1), 'tail','left');
        p_adj = min(p_raw*nAdj, 1);
        fprintf('  %s vs %s : median diff %+.3f um | ranksum p=%.3g, Bonferroni p=%.3g %s (n=%d vs n=%d)\n', ...
            binLabels{b}, binLabels{b-1}, median(D(vb,b))-median(D(va,b-1)), p_raw, p_adj, pstars(p_adj), sum(vb), sum(va));
        % (b) paired signrank : types defined in both bins (type as its own control -> removes between-type variance)
        vpair = vb & va;
        if sum(vpair) >= 2
            ps_raw = signrank(D(vpair,b), D(vpair,b-1), 'tail','left');
            ps_adj = min(ps_raw*nAdj, 1);
            fprintf('      |- paired signrank : median diff %+.3f um | p=%.3g, Bonferroni p=%.3g %s (n=%d pairs)\n', ...
                median(D(vpair,b)-D(vpair,b-1)), ps_raw, ps_adj, pstars(ps_adj), sum(vpair));
        else
            fprintf('      |- paired signrank : too few common types (n=%d pairs)\n', sum(vpair));
        end
    else
        fprintf('  %s vs %s : too few comparable types (n=%d vs n=%d)\n', ...
            binLabels{b}, binLabels{b-1}, sum(vb), sum(va));
    end
end
fprintf('====================================================================\n');

%% ===== Console summary : connected vs unconnected (Figure 3) =====
fprintf('\n========== connected vs unconnected (ratio descending) ==========\n');
for g = 1:nT_cu
    fprintf('  %-14s | unconn=%.2f um  conn=%.2f um  ratio=%.3f\n', ...
        types_cu{g}, Du_cu(g), Dc_cu(g), ratio_cu(g));
end
fprintf('  overall mean ratio = %.3f (>1 means connected side is closer)\n', mean(ratio_cu));
if nT_cu >= 2
    p_cu = signrank(ratio_cu, 1, 'tail','right');   % H1: ratio > 1
    fprintf('  signrank(ratio > 1) p = %.3g | median ratio = %.3f\n', p_cu, median(ratio_cu));
end
fprintf('==================================================================\n');
fprintf('Drawn: Fig1 (2G) / Fig2 (2H) / Fig3 (S1D).\n');

%% ===== Helper : p-value -> stars =====
function s = pstars(p)
    if p < 0.001
        s = '***';
    elseif p < 0.01
        s = '**';
    elseif p < 0.05
        s = '*';
    else
        s = 'n.s.';
    end
end
