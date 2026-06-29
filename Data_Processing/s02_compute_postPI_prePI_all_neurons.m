%% s02_compute_postPI_prePI_all_neurons
% NPI means Neuron Polarity Index including postPI and prePI
% Compute per-neuron synapse distributions across neuropil groups
% (right optic lobe, left optic lobe, central brain). These In/Out synapse
% counts are the basis for the postPI / prePI indices computed later for
% every neuron.
clear all; close all; clc

baseDir = fileparts(fileparts(mfilename('fullpath')));  % resolved from this script's location (portable; no absolute paths)

% Make sure the output folder exists
outDir = fullfile(baseDir, 'Processed_Data');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

opt = detectImportOptions(fullfile(baseDir, 'Codex_Data', 'connections_no_threshold.csv'));
opt = setvartype(opt,'pre_root_id','int64');
opt = setvartype(opt,'post_root_id','int64');
FAFBConnections = readtable(fullfile(baseDir, 'Codex_Data', 'connections_no_threshold.csv'),opt);

opt = detectImportOptions(fullfile(baseDir, 'Codex_Data', 'classification.csv'));
opt = setvartype(opt,'root_id','int64');
FAFBClassification = readtable(fullfile(baseDir, 'Codex_Data', 'classification.csv'),opt);

opt = detectImportOptions(fullfile(baseDir, 'Codex_Data', 'consolidated_cell_types.csv'));
opt = setvartype(opt,'root_id','int64');
FAFBConsolidatedTypes = readtable(fullfile(baseDir, 'Codex_Data', 'consolidated_cell_types.csv'),opt);


%% Define neuropil groups
FAFBNeuropils=unique(FAFBConnections.neuropil);
FAFBNeuropil_OpticLobeRight={'AME_R','ME_R','LO_R','LOP_R','LA_R'};
FAFBNeuropil_OpticLobeLeft={'AME_L','ME_L','LO_L','LOP_L','LA_L'};
FAFBNeuropil_Central=FAFBNeuropils(~ismember(FAFBNeuropils,FAFBNeuropil_OpticLobeRight));
FAFBNeuropil_Central=FAFBNeuropil_Central(~ismember(FAFBNeuropil_Central,FAFBNeuropil_OpticLobeLeft));
FAFBNeuropil_Central=FAFBNeuropil_Central(~ismember(FAFBNeuropil_Central,'UNASGD'));

%% Build the table
FAFBNPIs=table(FAFBConsolidatedTypes.root_id,FAFBConsolidatedTypes.primary_type,'VariableNames',{'root_id','type'});
%%
for i=1:1:size(FAFBNPIs,1)
    root_id=FAFBNPIs.root_id(i);
    %%%
    idx_classification=FAFBClassification.root_id==root_id;
    FAFBNPIs.superclass{i}=FAFBClassification.super_class{idx_classification};
    FAFBNPIs.side{i}=FAFBClassification.side{idx_classification};

    %%% outgoing (pre) connections
    pre_idx=FAFBConnections.pre_root_id==root_id;
    pre_connections=FAFBConnections(pre_idx,:);

    Out_Synapse_Optic_R=sum(pre_connections.syn_count(ismember(pre_connections.neuropil,FAFBNeuropil_OpticLobeRight)));
    Out_Synapse_Optic_L=sum(pre_connections.syn_count(ismember(pre_connections.neuropil,FAFBNeuropil_OpticLobeLeft)));
    Out_Synapse_Central=sum(pre_connections.syn_count(ismember(pre_connections.neuropil,FAFBNeuropil_Central)));

    %%% incoming (post) connections
    post_idx=FAFBConnections.post_root_id==root_id;
    post_connections=FAFBConnections(post_idx,:);

    In_Synapse_Optic_R=sum(post_connections.syn_count(ismember(post_connections.neuropil,FAFBNeuropil_OpticLobeRight)));
    In_Synapse_Optic_L=sum(post_connections.syn_count(ismember(post_connections.neuropil,FAFBNeuropil_OpticLobeLeft)));
    In_Synapse_Central=sum(post_connections.syn_count(ismember(post_connections.neuropil,FAFBNeuropil_Central)));


    FAFBNPIs.In_Synapse_Optic_R(i)=In_Synapse_Optic_R;
    FAFBNPIs.Out_Synapse_Optic_R(i)=Out_Synapse_Optic_R;
    FAFBNPIs.In_Synapse_Central(i)=In_Synapse_Central;
    FAFBNPIs.Out_Synapse_Central(i)=Out_Synapse_Central;
    FAFBNPIs.In_Synapse_Optic_L(i)=In_Synapse_Optic_L;
    FAFBNPIs.Out_Synapse_Optic_L(i)=Out_Synapse_Optic_L;

    i/size(FAFBNPIs,1)*100
end

%% Save result
save(fullfile(outDir, 'FAFB_NPI_thr0.mat'),'FAFBNPIs');
