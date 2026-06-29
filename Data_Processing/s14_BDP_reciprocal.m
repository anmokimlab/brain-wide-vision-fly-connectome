%% s14 — Input/output partner reciprocity (weighted Jaccard) for FFP / FBP / BDP
% For every right-hemisphere FFP, FBP, real-bidirectional (BDP), and other
% bidirectional (Others = optic + central BD) neuron, this script measures how much
% its input-partner set overlaps with its output-partner set using the
% synapse-weighted Jaccard index (Sum min / Sum max over the union of partner neuron
% ids). The per-neuron values are averaged per cell type, and the three example BDP
% types (LC9 / LT43 / LT52) are kept per neuron in `Want`.
%
% It saves Type_FFP, Type_FBP, Type_BDP, Type_Others, and Want to
% Processed_Data/BDP_reciprocity.mat, read by
% Figures/fig_5F_G_segregation_reciprocal.m.

clear all; clc; close all

baseDir = fileparts(fileparts(mfilename('fullpath')));  % resolved from this script's location (portable; no absolute paths)

% Codex connectivity (force root_id columns to int64)
opt = detectImportOptions(fullfile(baseDir, 'Codex_Data', 'connections_no_threshold.csv'));
opt = setvartype(opt, 'pre_root_id', 'int64');
opt = setvartype(opt, 'post_root_id', 'int64');
FAFBConnections = readtable(fullfile(baseDir, 'Codex_Data', 'connections_no_threshold.csv'), opt);

% FFP / FBP / BDP neuron classification, saved by
% Figures/fig_1D_E_postPI_prePI_right_neurons.m
load(fullfile(baseDir, 'Processed_Data', 'right_neurons_thr0.mat'), ...
    'RightFFP_NPIs', 'RightFBP_NPIs', 'RightBDP_real_NPIs', ...
    'RightBDP_optic_NPIs', 'RightBDP_central_NPIs')

% Other bidirectional neurons = optic + central BD (everything but real BD)
RightOthers_NPIs = [RightBDP_optic_NPIs; RightBDP_central_NPIs];

%% Per-neuron reciprocity (weighted Jaccard), aggregated per cell type
RightFFP_NPIs.weighted_Jaccard = reciprocity_for_neurons(RightFFP_NPIs.root_id, FAFBConnections);
RightFBP_NPIs.weighted_Jaccard = reciprocity_for_neurons(RightFBP_NPIs.root_id, FAFBConnections);
RightBDP_real_NPIs.weighted_Jaccard = reciprocity_for_neurons(RightBDP_real_NPIs.root_id, FAFBConnections);
RightOthers_NPIs.weighted_Jaccard = reciprocity_for_neurons(RightOthers_NPIs.root_id, FAFBConnections);

Type_FFP = summarize_reciprocity(RightFFP_NPIs);
Type_FBP = summarize_reciprocity(RightFBP_NPIs);
Type_BDP = summarize_reciprocity(RightBDP_real_NPIs);
Type_Others = summarize_reciprocity(RightOthers_NPIs);

%% Example BDP neuron types (kept per neuron)
WantToSee = {'LC9', 'LT43', 'LT52'};
Want = struct('type', {}, 'root_ids', {}, 'weighted_Jaccard', {});
for k = 1:numel(WantToSee)
    rid = RightBDP_real_NPIs.root_id(strcmp(RightBDP_real_NPIs.type, WantToSee{k}));
    Want(k).type = WantToSee{k};
    Want(k).root_ids = rid(:);
    Want(k).weighted_Jaccard = reciprocity_for_neurons(rid(:), FAFBConnections);
end

%% Save
save(fullfile(baseDir, 'Processed_Data', 'BDP_reciprocity.mat'), ...
    'Type_FFP', 'Type_FBP', 'Type_BDP', 'Type_Others', 'Want')

%% Local functions
function Jw = reciprocity_for_neurons(root_ids, FAFBConnections)
% Weighted Jaccard between each neuron's input-partner and output-partner
% synapse-weight vectors.
root_ids = root_ids(:);
Jw = zeros(numel(root_ids), 1);
for j = 1:numel(root_ids)
    rid = root_ids(j);

    In  = FAFBConnections(FAFBConnections.post_root_id == rid, :);   % upstream partners
    Out = FAFBConnections(FAFBConnections.pre_root_id  == rid, :);   % downstream partners

    [inIDs,  ~, ic_in]  = unique(In.pre_root_id);
    [outIDs, ~, ic_out] = unique(Out.post_root_id);
    inW  = accumarray(ic_in,  In.syn_count,  [numel(inIDs)  1]);
    outW = accumarray(ic_out, Out.syn_count, [numel(outIDs) 1]);

    Jw(j) = weighted_jaccard(inIDs, inW, outIDs, outW);
end
end

function Jw = weighted_jaccard(idA, wA, idB, wB)
% Sum(min) / Sum(max) over the union of ids (wA, wB >= 0).
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

function Tsum = summarize_reciprocity(T)
% Per cell type: member root_ids and the mean weighted Jaccard.
[types, ~, ic] = unique(T.type);
Tsum = table(types, 'VariableNames', {'type'});
for i = 1:numel(types)
    idx = (ic == i);
    Tsum.root_id{i} = T.root_id(idx);
    Tsum.weighted_Jaccard(i) = mean(T.weighted_Jaccard(idx), 'omitnan');
end
end
