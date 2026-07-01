%% s06_compute_FFP_optic_central_fan_in_out (MCNS)
% MCNS analogue of the FAFB Data_Processing/s05_compute_FFP_optic_central_fan_in_out.m
% (without the graph-centrality columns, which have no MCNS equivalent here).
%
% For the right FFP projection neurons and all right-hemisphere optic / central
% neurons, summarize their input/output partners grouped by partner type
% (partner-neuron count, partner-type count, synapse count), then aggregate per
% cell type. Precursor data processing for MCNS Figure 2.
clear all; close all; clc;

baseDir = fileparts(fileparts(mfilename('fullpath')));  % resolved from this script's location (portable; no absolute paths)

% RightFFP_NPIs from MCNS/Figures/fig_1_postPI_prePI_right_neurons.m
load(fullfile(baseDir, 'Processed_Data', 'right_neurons_thr0.mat'), 'RightFFP_NPIs')

opt = detectImportOptions(fullfile(baseDir, 'MCNS_Data', 'male-cns-v0.9-classification.csv'));
opt = setvartype(opt,'root_id','int64');
MCNSClassification = readtable(fullfile(baseDir, 'MCNS_Data', 'male-cns-v0.9-classification.csv'),opt);

opt = detectImportOptions(fullfile(baseDir, 'MCNS_Data', 'male-cns-v0.9-primary-types.csv'));
opt = setvartype(opt,'root_id','int64');
MCNSConsolidatedTypes = readtable(fullfile(baseDir, 'MCNS_Data', 'male-cns-v0.9-primary-types.csv'),opt);

opt = detectImportOptions(fullfile(baseDir, 'MCNS_Data', 'male-cns-v0.9-all-connections.csv'));
opt = setvartype(opt,'pre_root_id','int64');
opt = setvartype(opt,'post_root_id','int64');
MCNSConnections = readtable(fullfile(baseDir, 'MCNS_Data', 'male-cns-v0.9-all-connections.csv'),opt);

%% Per-neuron input/output partner-type tables (right FFP / optic / central)
RightFFP_InOut = compute_partner_in_out(RightFFP_NPIs.root_id, MCNSConnections, MCNSConsolidatedTypes);

Optic_root_id = MCNSClassification.root_id(strcmp(MCNSClassification.super_class,'optic') & ~strcmp(MCNSClassification.side,'left'));
Optic_InOut = compute_partner_in_out(Optic_root_id, MCNSConnections, MCNSConsolidatedTypes);

Central_root_id = MCNSClassification.root_id(strcmp(MCNSClassification.super_class,'central') & ~strcmp(MCNSClassification.side,'left'));
Central_InOut = compute_partner_in_out(Central_root_id, MCNSConnections, MCNSConsolidatedTypes);

%% Per-neuron fan-in / fan-out summary, then aggregate per cell type
RightFFP_InOut = add_neuron_stats(RightFFP_InOut);
Optic_InOut    = add_neuron_stats(Optic_InOut);
Central_InOut  = add_neuron_stats(Central_InOut);

type_RightFFP_InOut = aggregate_by_type(RightFFP_InOut);
type_Optic_InOut    = aggregate_by_type(Optic_InOut);
type_Central_InOut  = aggregate_by_type(Central_InOut);

%% Save the per-type fan-in / fan-out summaries
save(fullfile(baseDir, 'Processed_Data', 'FFP_optic_central_fan_in_out.mat'), ...
    'type_Optic_InOut', 'type_Central_InOut', 'type_RightFFP_InOut');

%% ===================== Local functions =====================
function T = compute_partner_in_out(rootIds, Connections, MCNSConsolidatedTypes)
% For each neuron in rootIds, tabulate input and output partners grouped by partner
% cell type. Returns a table with columns root_id, type, In, Out, where In/Out are
% n x 4 cells: {partner type, partner root_ids, partner-neuron count, synapse count}.
rootIds = rootIds(:);
% Keep only connections involving this group (speeds up the per-neuron filtering)
Connections(~(ismember(Connections.pre_root_id, rootIds) | ismember(Connections.post_root_id, rootIds)), :) = [];

n = numel(rootIds);
tempType = cell(n,1); tempIn = cell(n,1); tempOut = cell(n,1);
parfor i = 1:n
    rid = rootIds(i);
    idx_consol = MCNSConsolidatedTypes.root_id == rid;
    if any(idx_consol), tempType{i} = MCNSConsolidatedTypes.primary_type{idx_consol}; else, tempType{i} = ''; end

    InConnections  = Connections(Connections.post_root_id == rid, :);
    OutConnections = Connections(Connections.pre_root_id  == rid, :);
    tempIn{i}  = partner_stats(InConnections,  'pre_root_id',  MCNSConsolidatedTypes);
    tempOut{i} = partner_stats(OutConnections, 'post_root_id', MCNSConsolidatedTypes);
end

T = table(rootIds, tempType, tempIn, tempOut, 'VariableNames', {'root_id','type','In','Out'});
end

function S = partner_stats(Connections, idCol, MCNSConsolidatedTypes)
% Group a neuron's partners (column idCol of Connections) by cell type:
% {type, partner root_ids, partner-neuron count, synapse count}, sorted by synapse count.
partners = unique(Connections.(idCol));
if isempty(partners), S = {}; return; end

types = cell(numel(partners),1);
for j = 1:numel(partners)
    ic = MCNSConsolidatedTypes.root_id == partners(j);
    if any(ic), types{j} = MCNSConsolidatedTypes.primary_type{ic}; else, types{j} = num2str(partners(j)); end
end

[uTypes,~,ic] = unique(types);
S = cell(numel(uTypes),4);
for j = 1:numel(uTypes)
    pid = partners(ic == j);
    S{j,1} = uTypes{j};
    S{j,2} = pid;
    S{j,3} = numel(pid);
    S{j,4} = sum(Connections.syn_count(ismember(Connections.(idCol), pid)));
end
S = sortrows(S, 4, 'descend');
end

function T = add_neuron_stats(T)
% Per neuron: total partner count, partner-type count, and their ratio (in and out).
for i = 1:size(T,1)
    inStats = T.In{i}; outStats = T.Out{i};
    if ~isempty(inStats)
        T.InNeuronNumber{i}     = sum(cell2mat(inStats(:,3)));
        T.InNeuronTypeNumber{i} = size(inStats,1);
        T.InNeuronRatio{i}      = T.InNeuronNumber{i} / T.InNeuronTypeNumber{i};
    end
    if ~isempty(outStats)
        T.OutNeuronNumber{i}     = sum(cell2mat(outStats(:,3)));
        T.OutNeuronTypeNumber{i} = size(outStats,1);
        T.OutNeuronRatio{i}      = T.OutNeuronNumber{i} / T.OutNeuronTypeNumber{i};
    end
end
end

function Tt = aggregate_by_type(T)
% Mean per-neuron fan-in / fan-out per cell type.
[types,~,ic] = unique(T.type);
Tt = table(types, 'VariableNames', {'type'});
for i = 1:numel(types)
    idx = (ic == i);
    Tt.InNeuronNumber{i}      = mean(cell2mat(T.InNeuronNumber(idx)),     'omitmissing');
    Tt.InNeuronTypeNumber{i}  = mean(cell2mat(T.InNeuronTypeNumber(idx)), 'omitmissing');
    Tt.InNeuronRatio{i}       = mean(cell2mat(T.InNeuronRatio(idx)),      'omitmissing');
    Tt.OutNeuronNumber{i}     = mean(cell2mat(T.OutNeuronNumber(idx)),    'omitmissing');
    Tt.OutNeuronTypeNumber{i} = mean(cell2mat(T.OutNeuronTypeNumber(idx)),'omitmissing');
    Tt.OutNeuronRatio{i}      = mean(cell2mat(T.OutNeuronRatio(idx)),     'omitmissing');
end
end
