%% FF/FB/BI 별 히스토그램도 여기서 그림
clear all; close all; clc

opt = detectImportOptions('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\connections_no_threshold.csv');
opt = setvartype(opt,'pre_root_id','int64');
opt = setvartype(opt,'post_root_id','int64');
FAFBConnections = readtable('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\connections_no_threshold.csv',opt);

opt = detectImportOptions('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\classification.csv');
opt = setvartype(opt,'root_id','int64');
FAFBClassification = readtable('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\classification.csv',opt);

opt = detectImportOptions('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\consolidated_cell_types.csv');
opt = setvartype(opt,'root_id','int64');
FAFBConsolidatedTypes = readtable('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\consolidated_cell_types.csv',opt);

%% Neuropil 데이터 만들기

FAFBNeuropils=unique(FAFBConnections.neuropil);
FAFBNeuropil_OpticLobeRight={'AME_R','ME_R','LO_R','LOP_R','LA_R'};
FAFBNeuropil_OpticLobeLeft={'AME_L','ME_L','LO_L','LOP_L','LA_L'};
FAFBNeuropil_Central=FAFBNeuropils(~ismember(FAFBNeuropils,FAFBNeuropil_OpticLobeRight));
FAFBNeuropil_Central=FAFBNeuropil_Central(~ismember(FAFBNeuropil_Central,FAFBNeuropil_OpticLobeLeft));
FAFBNeuropil_Central=FAFBNeuropil_Central(~ismember(FAFBNeuropil_Central,'UNASGD'));

%% table 만들기
FAFBNPIs=table(FAFBConsolidatedTypes.root_id,FAFBConsolidatedTypes.primary_type,'VariableNames',{'root_id','type'});
%%
for i=1:1:size(FAFBNPIs,1)
    root_id=FAFBNPIs.root_id(i);
    %%%
    idx_classification=FAFBClassification.root_id==root_id;
    FAFBNPIs.superclass{i}=FAFBClassification.super_class{idx_classification};
    FAFBNPIs.side{i}=FAFBClassification.side{idx_classification};
 
    %%% find is pre?
    pre_idx=FAFBConnections.pre_root_id==root_id;
    pre_connections=FAFBConnections(pre_idx,:);

    Out_Synapse_Optic_R=sum(pre_connections.syn_count(ismember(pre_connections.neuropil,FAFBNeuropil_OpticLobeRight)));
    Out_Synapse_Optic_L=sum(pre_connections.syn_count(ismember(pre_connections.neuropil,FAFBNeuropil_OpticLobeLeft)));
    Out_Synapse_Central=sum(pre_connections.syn_count(ismember(pre_connections.neuropil,FAFBNeuropil_Central)));

    % Out_Partner_Optic_R=sum(ismember(pre_connections.neuropil,FAFBNeuropil_OpticLobeRight));
    % Out_Partner_Optic_L=sum(ismember(pre_connections.neuropil,FAFBNeuropil_OpticLobeLeft));
    % Out_Partner_Central=sum(ismember(pre_connections.neuropil,FAFBNeuropil_Central));

    post_idx=FAFBConnections.post_root_id==root_id;
    post_connections=FAFBConnections(post_idx,:);

    In_Synapse_Optic_R=sum(post_connections.syn_count(ismember(post_connections.neuropil,FAFBNeuropil_OpticLobeRight)));
    In_Synapse_Optic_L=sum(post_connections.syn_count(ismember(post_connections.neuropil,FAFBNeuropil_OpticLobeLeft)));
    In_Synapse_Central=sum(post_connections.syn_count(ismember(post_connections.neuropil,FAFBNeuropil_Central)));

    % In_Partner_Optic_R=sum(ismember(post_connections.neuropil,FAFBNeuropil_OpticLobeRight));
    % In_Partner_Optic_L=sum(ismember(post_connections.neuropil,FAFBNeuropil_OpticLobeLeft));
    % In_Partner_Central=sum(ismember(post_connections.neuropil,FAFBNeuropil_Central));


    FAFBNPIs.In_Synapse_Optic_R(i)=In_Synapse_Optic_R;
    FAFBNPIs.Out_Synapse_Optic_R(i)=Out_Synapse_Optic_R;
    FAFBNPIs.In_Synapse_Central(i)=In_Synapse_Central;
    FAFBNPIs.Out_Synapse_Central(i)=Out_Synapse_Central;
    FAFBNPIs.In_Synapse_Optic_L(i)=In_Synapse_Optic_L;
    FAFBNPIs.Out_Synapse_Optic_L(i)=Out_Synapse_Optic_L;


    % FAFBNPIs.In_Partner_Optic_R(i)=In_Partner_Optic_R;
    % FAFBNPIs.Out_Partner_Optic_R(i)=Out_Partner_Optic_R;
    % FAFBNPIs.In_Partner_Central(i)=In_Partner_Central;
    % FAFBNPIs.Out_Partner_Central(i)=Out_Partner_Central;
    % FAFBNPIs.In_Partner_Optic_L(i)=In_Partner_Optic_L;
    % FAFBNPIs.Out_Partner_Optic_L(i)=Out_Partner_Optic_L;



    i/size(FAFBNPIs,1)*100
end
% save("FAFB_NPI_Thr0.mat","FAFBNPIs");