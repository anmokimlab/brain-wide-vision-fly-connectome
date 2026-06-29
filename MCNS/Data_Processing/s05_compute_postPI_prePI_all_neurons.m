%% s05_compute_postPI_prePI_all_neurons (MCNS)
% NPI means Neuron Polarity Index including postPI and prePI.
% Compute per-neuron synapse distributions across neuropil groups
% (right optic lobe, left optic lobe, central brain) for the male-CNS dataset.
% These In/Out synapse counts are the basis for the postPI / prePI indices
% computed later for every neuron.
clear all; close all; clc

baseDir = fileparts(fileparts(mfilename('fullpath')));  % resolved from this script's location (portable; no absolute paths)

% Make sure the output folder exists
outDir = fullfile(baseDir, 'Processed_Data');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

opt = detectImportOptions(fullfile(baseDir, 'MCNS_Data', 'male-cns-v0.9-primary-types.csv'));
opt = setvartype(opt,'root_id','int64');
MCNSConsolidatedTypes = readtable(fullfile(baseDir, 'MCNS_Data', 'male-cns-v0.9-primary-types.csv'),opt);

opt = detectImportOptions(fullfile(baseDir, 'MCNS_Data', 'male-cns-v0.9-all-connections.csv'));
opt = setvartype(opt,'pre_root_id','int64');
opt = setvartype(opt,'post_root_id','int64');
MCNSConnections = readtable(fullfile(baseDir, 'MCNS_Data', 'male-cns-v0.9-all-connections.csv'),opt);

opt = detectImportOptions(fullfile(baseDir, 'MCNS_Data', 'male-cns-v0.9-classification.csv'));
opt = setvartype(opt,'root_id','int64');
MCNSClassification = readtable(fullfile(baseDir, 'MCNS_Data', 'male-cns-v0.9-classification.csv'),opt);

%% Define neuropil groups
MCNSNeuropils=unique(MCNSConnections.neuropil);
MCNSNeuropil_OpticLobeRight={'AME(R)','ME(R)','LO(R)','LOP(R)','LA(R)','Optic-unspecified(R)'};
MCNSNeuropil_OpticLobeLeft={'AME(L)','ME(L)','LO(L)','LOP(L)','LA(L)','Optic-unspecified(L)'};
MCNSNeuropil_Central=MCNSNeuropils(~ismember(MCNSNeuropils,MCNSNeuropil_OpticLobeRight));
MCNSNeuropil_Central=MCNSNeuropil_Central(~ismember(MCNSNeuropil_Central,MCNSNeuropil_OpticLobeLeft));
MCNSNeuropil_Central=MCNSNeuropil_Central(~ismember(MCNSNeuropil_Central,'NotPrimary'));

%% Build the table
MCNSNPIs=table(MCNSConsolidatedTypes.root_id,MCNSConsolidatedTypes.primary_type,MCNSConsolidatedTypes.flywireType, ...
    'VariableNames',{'root_id','type','flywireType'});
%%
for i=1:1:size(MCNSNPIs,1)
    root_id=MCNSNPIs.root_id(i);
    %%%
    idx_classification=MCNSClassification.root_id==root_id;
    MCNSNPIs.superclass{i}=MCNSClassification.super_class{idx_classification};
    MCNSNPIs.side{i}=MCNSClassification.side{idx_classification};

    %%% outgoing (pre) connections
    pre_idx=MCNSConnections.pre_root_id==root_id;
    pre_connections=MCNSConnections(pre_idx,:);

    Out_Synapse_Optic_R=sum(pre_connections.syn_count(ismember(pre_connections.neuropil,MCNSNeuropil_OpticLobeRight)));
    Out_Synapse_Optic_L=sum(pre_connections.syn_count(ismember(pre_connections.neuropil,MCNSNeuropil_OpticLobeLeft)));
    Out_Synapse_Central=sum(pre_connections.syn_count(ismember(pre_connections.neuropil,MCNSNeuropil_Central)));

    %%% incoming (post) connections
    post_idx=MCNSConnections.post_root_id==root_id;
    post_connections=MCNSConnections(post_idx,:);

    In_Synapse_Optic_R=sum(post_connections.syn_count(ismember(post_connections.neuropil,MCNSNeuropil_OpticLobeRight)));
    In_Synapse_Optic_L=sum(post_connections.syn_count(ismember(post_connections.neuropil,MCNSNeuropil_OpticLobeLeft)));
    In_Synapse_Central=sum(post_connections.syn_count(ismember(post_connections.neuropil,MCNSNeuropil_Central)));


    MCNSNPIs.In_Synapse_Optic_R(i)=In_Synapse_Optic_R;
    MCNSNPIs.Out_Synapse_Optic_R(i)=Out_Synapse_Optic_R;
    MCNSNPIs.In_Synapse_Central(i)=In_Synapse_Central;
    MCNSNPIs.Out_Synapse_Central(i)=Out_Synapse_Central;
    MCNSNPIs.In_Synapse_Optic_L(i)=In_Synapse_Optic_L;
    MCNSNPIs.Out_Synapse_Optic_L(i)=Out_Synapse_Optic_L;

    i/size(MCNSNPIs,1)*100
end

%% Save result
save(fullfile(outDir, 'MCNS_NPI_thr0.mat'),'MCNSNPIs');
