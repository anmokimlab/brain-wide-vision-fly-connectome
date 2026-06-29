%% =================================================================
%  BDP reciprocity within the optic lobe (OL) vs central brain (CB)  (neuron-level)
%
%  For each right-side real-bidirectional (BDP) neuron, computes the weighted
%  Jaccard index (WJI) between its input-partner and output-partner neurons
%  WITHIN the same compartment, then averages the per-neuron values per cell type:
%    WJI_OL_total : OL-input partners vs OL-output partners (per neuron, syn-weighted)
%    WJI_CB_total : CB-input partners vs CB-output partners (per neuron, syn-weighted)
%
%  Per neuron, the WJI is Sum(min)/Sum(max) over the union of partner neuron ids,
%  with weights = synapse counts to/from each partner, restricted to connections in
%  the given compartment (OL or CB; UNASGD is never part of either set). This is the
%  same neuron-level convention as s14_BDP_reciprocal.m. A neuron with no connections
%  in a compartment yields NaN (omitted from the per-type mean).
%
%  Filtering: per-type root_ids resolved from consolidated_cell_types.primary_type and
%  intersected with RightBDP_real_NPIs.root_id (right side only); partners are
%  individual neurons (no cell-type grouping); 1-hop only.
%
%  Produces two figures:
%    Figure 1 (panel S3A) — boxplot of WJI_OL_total vs WJI_CB_total (signed-rank p)
%    Figure 2 (panel S3B) — scatter of WJI_OL_total vs WJI_CB_total across BDP types
% =================================================================

clear all; close all; clc

baseDir = fileparts(fileparts(mfilename('fullpath')));  % resolved from this script's location (portable; no absolute paths)

%% --- Load data ---
% BDP (real-bidirectional) classification, from fig_1D_E_postPI_prePI_right_neurons.m
load(fullfile(baseDir, 'Processed_Data', 'right_neurons_thr0.mat'), 'RightBDP_real_type', 'RightBDP_real_NPIs')

opt = detectImportOptions(fullfile(baseDir, 'Codex_Data', 'connections_no_threshold.csv'));
opt = setvartype(opt,'pre_root_id','int64');
opt = setvartype(opt,'post_root_id','int64');
FAFBConnections = readtable(fullfile(baseDir, 'Codex_Data', 'connections_no_threshold.csv'), opt);

opt = detectImportOptions(fullfile(baseDir, 'Codex_Data', 'consolidated_cell_types.csv'));
opt = setvartype(opt,'root_id','int64');
FAFBConsolidated_type = readtable(fullfile(baseDir, 'Codex_Data', 'consolidated_cell_types.csv'), opt);

%% --- Neuropil sets ---
opticLobeNPL = {'ME_R','AME_R','LO_R','LOP_R','LA_R', ...
                'ME_L','AME_L','LO_L','LOP_L','LA_L'};
excludeNPL = {'UNASGD'};
allNPL = unique(FAFBConnections.neuropil);
centralBrainNPL = setdiff(setdiff(allNPL, opticLobeNPL), excludeNPL);

typeList = RightBDP_real_type.type(:);

%% --- Compute neuron-level reciprocal WJI (OL / CB), averaged per type ---
fprintf('Computing neuron-level reciprocal WJI (OL / CB) per type ...\n');
rightRoots = int64(RightBDP_real_NPIs.root_id);

N = numel(typeList);
WJI_OL_total = nan(N,1);
WJI_CB_total = nan(N,1);

for i = 1:N
    selfType = typeList{i};

    % Resolve this type's right-side member neurons
    root_ids = unique(FAFBConsolidated_type.root_id( ...
        strcmpi(FAFBConsolidated_type.primary_type, selfType)));
    root_ids = intersect(int64(root_ids), rightRoots);

    % Per-neuron WJI within each compartment, then average across neurons
    jOL = reciprocity_within_npl(root_ids, FAFBConnections, opticLobeNPL);
    jCB = reciprocity_within_npl(root_ids, FAFBConnections, centralBrainNPL);
    WJI_OL_total(i) = mean(jOL, 'omitnan');
    WJI_CB_total(i) = mean(jCB, 'omitnan');

    if mod(i,10)==0 || i==N
        fprintf('  %d/%d done\n', i, N);
    end
end

%% --- Attach columns to RightBDP_real_type ---
RightBDP_real_type.WJI_OL_total = WJI_OL_total;
RightBDP_real_type.WJI_CB_total = WJI_CB_total;

fprintf('Done.\n');
disp(RightBDP_real_type(1:min(10,height(RightBDP_real_type)),:));

%% --- Pull vectors & stats ---
ol_t = RightBDP_real_type.WJI_OL_total;
cb_t = RightBDP_real_type.WJI_CB_total;
typenames = RightBDP_real_type.type;
N = numel(ol_t);

% Paired (signed-rank) and unpaired (rank-sum)
valid = ~isnan(ol_t) & ~isnan(cb_t);
[p_sr, h_sr, stats_sr] = signrank(ol_t(valid), cb_t(valid));
diff_vec = ol_t(valid) - cb_t(valid);

ol_valid = ol_t(~isnan(ol_t));
cb_valid = cb_t(~isnan(cb_t));
[p_rs, h_rs, stats_rs] = ranksum(ol_valid, cb_valid);

med_ol = median(ol_t,'omitnan');  med_cb = median(cb_t,'omitnan');
mu_ol  = mean(ol_t,'omitnan');    mu_cb  = mean(cb_t,'omitnan');

fprintf('\nPopulation summary (n = %d types, NaN omitted):\n', N);
fprintf('  OL total  median = %.3f, mean = %.3f\n', med_ol, mu_ol);
fprintf('  CB total  median = %.3f, mean = %.3f\n', med_cb, mu_cb);

fprintf('\nWilcoxon signed-rank (OL_total vs CB_total, paired)\n');
fprintf('  n pairs = %d\n', sum(valid));
fprintf('  median (OL - CB) = %.3f\n', median(diff_vec));
fprintf('  p = %.4g, h = %d', p_sr, h_sr);
if isfield(stats_sr,'zval')
    fprintf(', z = %.3f, signedrank = %d\n', stats_sr.zval, stats_sr.signedrank);
else
    fprintf(', signedrank = %d\n', stats_sr.signedrank);
end

fprintf('\nWilcoxon ranksum (OL_total vs CB_total, unpaired)\n');
fprintf('  n_OL = %d, n_CB = %d\n', numel(ol_valid), numel(cb_valid));
fprintf('  p = %.4g, h = %d', p_rs, h_rs);
if isfield(stats_rs,'zval')
    fprintf(', z = %.3f, ranksum = %d\n', stats_rs.zval, stats_rs.ranksum);
else
    fprintf(', ranksum = %d\n', stats_rs.ranksum);
end

colScatter = [0.30 0.55 0.75];
nLabel     = 5;   % number of top types to label

%% --- Figure 1 (panel S3A): boxplot of OL vs CB reciprocal WJI ---
figure(1); set(gcf,'Color','w'); hold on; box on; grid on;
boxplot([ol_t, cb_t], 'Labels',{'OL_{total}','CB_{total}'}, ...
        'Symbol','o', 'Whisker',1.5, 'Widths',0.55);
set(gca,'TickLabelInterpreter','tex');
plot([1 2],[mu_ol mu_cb], 'd','MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',9);  % mean
ylim([0 0.5]); ylabel('WJI');

% signed-rank significance bracket
yTop = 0.45;
plot([1 1 2 2],[yTop-0.02 yTop yTop yTop-0.02],'k-','LineWidth',1);
if p_sr >= 0.05, sigStr = 'n.s.'; else, sigStr = '*'; end
text(1.5, yTop+0.02, sprintf('%s (p = %.3f)', sigStr, p_sr), 'HorizontalAlignment','center','FontSize',10);

legend({'mean'},'Location','southwest','AutoUpdate','off');
title(sprintf('S3A  |  Wilcoxon signed-rank, n = %d pairs', sum(valid)));
set(gca,'TickDir','out','Box','off')

%% --- Figure 2 (panel S3B): scatter of OL vs CB reciprocal WJI ---
figure(2); set(gcf,'Color','w'); hold on; box on; grid on;
scatter(ol_t, cb_t, 40, colScatter, 'filled', 'MarkerFaceAlpha',0.6, 'MarkerEdgeColor','none');
plot([0 1],[0 1],':','Color',[0.5 0.5 0.5]);              % y = x
xline(med_ol,'--','Color',[0.4 0.4 0.4]);
yline(med_cb,'--','Color',[0.4 0.4 0.4]);
plot(med_ol, med_cb,'kx','MarkerSize',14,'LineWidth',2.5);    % median marker
% Label the most reciprocal types (by OL+CB sum-score)
labelTopN_byScore(ol_t, cb_t, typenames, ol_t + cb_t, nLabel);
xlim([0 0.5]); ylim([0 0.5]); axis square;
xlabel('WJI_{OL,total}','Interpreter','tex');
ylabel('WJI_{CB,total}','Interpreter','tex');
title(sprintf('S3B  |  median = (%.2f, %.2f),  n = %d', med_ol, med_cb, sum(valid)));
set(gca,'TickDir','out','Box','off')

%% =================================================================
%  Local functions
% =================================================================

function labelTopN_byScore(xv, yv, names, score, n)
% Label top-n points by an externally-supplied score.
    valid = ~isnan(xv) & ~isnan(yv) & ~isnan(score);
    idx_valid = find(valid);
    [~,order] = sort(score(idx_valid),'descend');
    take = idx_valid(order(1:min(n,numel(order))));
    for k = 1:numel(take)
        i = take(k);
        text(xv(i)+0.012, yv(i), strrep(names{i},'_','\_'), ...
             'FontSize',8,'Color',[0.2 0.2 0.2]);
    end
end

function Jw = reciprocity_within_npl(root_ids, FAFBConnections, nplSet)
% Per-neuron weighted Jaccard between a neuron's input-partner and output-partner
% (individual neuron) synapse-weight vectors, restricted to connections in nplSet.
% Returns NaN for neurons with no connections in the compartment.
    root_ids = root_ids(:);
    Jw = nan(numel(root_ids), 1);
    for j = 1:numel(root_ids)
        rid = root_ids(j);
        In  = FAFBConnections(FAFBConnections.post_root_id == rid & ...
                              ismember(FAFBConnections.neuropil, nplSet), :);
        Out = FAFBConnections(FAFBConnections.pre_root_id  == rid & ...
                              ismember(FAFBConnections.neuropil, nplSet), :);
        if isempty(In) && isempty(Out)
            Jw(j) = NaN;   % neuron has no connections in this compartment
            continue;
        end
        [inIDs,  ~, ic_in]  = unique(In.pre_root_id);
        [outIDs, ~, ic_out] = unique(Out.post_root_id);
        inW  = accumarray(ic_in,  In.syn_count,  [numel(inIDs)  1]);
        outW = accumarray(ic_out, Out.syn_count, [numel(outIDs) 1]);
        Jw(j) = weighted_jaccard(inIDs, inW, outIDs, outW);
    end
end

function Jw = weighted_jaccard(idA, wA, idB, wB)
% Sum(min) / Sum(max) over the union of partner ids (wA, wB >= 0).
    allIDs = union(idA, idB);
    wAfull = zeros(size(allIDs));
    wBfull = zeros(size(allIDs));
    [liaA, locA] = ismember(allIDs, idA);
    [liaB, locB] = ismember(allIDs, idB);
    wAfull(liaA) = wA(locA(liaA));
    wBfull(liaB) = wB(locB(liaB));
    den = sum(max(wAfull, wBfull));
    if den == 0
        Jw = 0;
    else
        Jw = sum(min(wAfull, wBfull)) / den;
    end
end
