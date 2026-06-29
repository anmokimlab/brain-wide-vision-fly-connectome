%% s14_BDP_reciprocal (MCNS)
% MCNS analogue of the FAFB Data_Processing/s14_BDP_reciprocal.m.
%
% For every right-hemisphere FFP, FBP, real-bidirectional (BDP), and other
% bidirectional (Others = optic + central BD) neuron, this script measures how much
% its input-partner set overlaps with its output-partner set using the
% synapse-weighted Jaccard index (Sum min / Sum max over the union of partner neuron
% ids). The per-neuron values are averaged per cell type.
%
% Unlike the FAFB version (which scans the connection table per neuron), the MCNS
% all-connections table is large, so each group's input/output partner weights are
% pre-aggregated once with groupsummary into lookup Maps before the per-neuron
% Jaccard. The three example BDP types (LC9 / LT43 / LT52) are kept per neuron in
% `Want`. It saves Type_FFP, Type_FBP, Type_BDP, Type_Others, and Want to
% MCNS/Processed_Data/BDP_reciprocity.mat.

clear all; clc; close all

baseDir = fileparts(fileparts(mfilename('fullpath')));  % resolved from this script's location (portable; no absolute paths)

% MCNS connectivity (force root_id columns to int64)
opt = detectImportOptions(fullfile(baseDir, 'MCNS_Data', 'male-cns-v0.9-all-connections.csv'));
opt = setvartype(opt, 'pre_root_id', 'int64');
opt = setvartype(opt, 'post_root_id', 'int64');
MCNSConnections = readtable(fullfile(baseDir, 'MCNS_Data', 'male-cns-v0.9-all-connections.csv'), opt);

% FFP / FBP / BDP neuron classification, saved by
% MCNS/Figures/fig_1_postPI_prePI_right_neurons.m
load(fullfile(baseDir, 'Processed_Data', 'right_neurons_thr0.mat'), ...
    'RightFFP_NPIs', 'RightFBP_NPIs', 'RightBDP_real_NPIs', ...
    'RightBDP_optic_NPIs', 'RightBDP_central_NPIs')

% Other bidirectional neurons = optic + central BD (everything but real BD)
RightOthers_NPIs = [RightBDP_optic_NPIs; RightBDP_central_NPIs];

%% Per-neuron reciprocity (weighted Jaccard), aggregated per cell type
RightFFP_NPIs.weighted_Jaccard      = reciprocity_for_neurons(RightFFP_NPIs.root_id, MCNSConnections);
RightFBP_NPIs.weighted_Jaccard      = reciprocity_for_neurons(RightFBP_NPIs.root_id, MCNSConnections);
RightBDP_real_NPIs.weighted_Jaccard = reciprocity_for_neurons(RightBDP_real_NPIs.root_id, MCNSConnections);
RightOthers_NPIs.weighted_Jaccard   = reciprocity_for_neurons(RightOthers_NPIs.root_id, MCNSConnections);

Type_FFP    = summarize_reciprocity(RightFFP_NPIs);
Type_FBP    = summarize_reciprocity(RightFBP_NPIs);
Type_BDP    = summarize_reciprocity(RightBDP_real_NPIs);
Type_Others = summarize_reciprocity(RightOthers_NPIs);

%% Example BDP neuron types (kept per neuron)
WantToSee = {'LC9', 'LT43', 'LT52'};
Want = struct('type', {}, 'root_ids', {}, 'weighted_Jaccard', {});
for k = 1:numel(WantToSee)
    rid = RightBDP_real_NPIs.root_id(strcmp(RightBDP_real_NPIs.type, WantToSee{k}));
    Want(k).type = WantToSee{k};
    Want(k).root_ids = rid(:);
    Want(k).weighted_Jaccard = reciprocity_for_neurons(rid(:), MCNSConnections);
end

%% Save
save(fullfile(baseDir, 'Processed_Data', 'BDP_reciprocity.mat'), ...
    'Type_FFP', 'Type_FBP', 'Type_BDP', 'Type_Others', 'Want')

%% Local functions
function Jw = reciprocity_for_neurons(root_ids, Connections)
% Weighted Jaccard between each neuron's input-partner and output-partner
% synapse-weight vectors. The input/output partner weights for the whole group are
% pre-aggregated once (per target, per partner) into lookup Maps for speed.
root_ids = root_ids(:);

Conn_In  = Connections(ismember(Connections.post_root_id, root_ids), :);   % upstream partners
Conn_Out = Connections(ismember(Connections.pre_root_id,  root_ids), :);   % downstream partners

Input_Map  = aggregate_id_map(Conn_In,  'post_root_id', 'pre_root_id');    % target -> [partner_id, weight]
Output_Map = aggregate_id_map(Conn_Out, 'pre_root_id',  'post_root_id');

Jw = zeros(numel(root_ids), 1);
for j = 1:numel(root_ids)
    rid = root_ids(j);
    if isKey(Input_Map, rid) && isKey(Output_Map, rid)
        A = Input_Map(rid);    % [partner_id, weight]
        B = Output_Map(rid);
        Jw(j) = weighted_jaccard(A(:,1), A(:,2), B(:,1), B(:,2));
    else
        Jw(j) = 0;   % a neuron missing inputs or outputs has zero overlap
    end
end
end

function M = aggregate_id_map(ConnTable, GroupKey, PartnerKey)
% target -> [partner_id, summed synapse weight] Map, summing syn_count over neuropils.
if isempty(ConnTable)
    M = containers.Map('KeyType','int64','ValueType','any');
    return
end
G = groupsummary(ConnTable, {GroupKey, PartnerKey}, 'sum', 'syn_count');
[g_idx, target_list] = findgroups(G.(GroupKey));
pack_func  = @(p, s) {[p, s]};
ValuesList = splitapply(pack_func, G.(PartnerKey), G.sum_syn_count, g_idx);
M = containers.Map(target_list, ValuesList);
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
