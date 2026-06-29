%% s09_FBP_out_opticlobes (MCNS)
% MCNS analogue of the FAFB Data_Processing/s10_FBP_out_opticlobes.m.
%
% For every FBP cell type, compute the fraction (%) of its output synapses
% landing in each optic-lobe region, collapsing left/right into AME / ME / LO /
% LOP (stored in the variable RightFBP_OutNeuropils). These values were written
% by hand to MCNS/Processed_Data/FBP_output_neuropils.xlsx, which was then used
% to classify each FBP neuron type by the optic-lobe region it targets; the
% resulting target grouping (Me / Lo / Lop / Multi) is reused by the downstream
% MCNS Figure 4 code.
clear all; close all; clc

baseDir = fileparts(fileparts(mfilename('fullpath')));  % resolved from this script's location (portable; no absolute paths)

% FFP / FBP / BDP neuron classification (provides RightFBP_by_type), saved by
% MCNS/Figures/fig_1_postPI_prePI_right_neurons.m
load(fullfile(baseDir, 'Processed_Data', 'right_neurons_thr0.mat'))

opt = detectImportOptions(fullfile(baseDir, 'MCNS_Data', 'male-cns-v0.9-primary-types.csv'));
opt = setvartype(opt,'root_id','int64');
MCNSConsolidatedTypes = readtable(fullfile(baseDir, 'MCNS_Data', 'male-cns-v0.9-primary-types.csv'),opt);

opt = detectImportOptions(fullfile(baseDir, 'MCNS_Data', 'male-cns-v0.9-all-connections.csv'));
opt = setvartype(opt,'pre_root_id','int64');
opt = setvartype(opt,'post_root_id','int64');
MCNSConnections = readtable(fullfile(baseDir, 'MCNS_Data', 'male-cns-v0.9-all-connections.csv'),opt);

%% List the unique FBP cell types and summarize their output neuropils
RightFBP_type={};
for i=1:1:size(RightFBP_by_type,1)
    RightFBP_type{i}=RightFBP_by_type.type{i};

end
RightFBP_type=unique(RightFBP_type)';

% seeConnection_by_region (local function, see bottom of file): for each cell type
% it returns the in/out partner types and neuropils with their total synapse counts.
RightFBP_Neuropils=seeConnection_by_region(RightFBP_type,MCNSConnections,MCNSConsolidatedTypes);

%% Convert per-neuropil output synapses into percentages
for i=1:1:size(RightFBP_Neuropils,1)
    temp=RightFBP_Neuropils.OutNeuropils{i,1};
    totalSynapse=sum(cell2mat(temp(:,2)));
    for j=1:1:size(temp,1)
        temp{j,3}=temp{j,2}/totalSynapse*100;
    end
    RightFBP_Neuropils.OutNeuropils{i,1}=temp;
end

%% Collapse left/right optic-lobe regions into AME / ME / LO / LOP per type
RightFBP_OutNeuropils = {};
RightFBP_OutNeuropils(:,1) = RightFBP_Neuropils.type;

for i = 1:size(RightFBP_OutNeuropils,1)
    temp = RightFBP_Neuropils.OutNeuropils{i,1};
    idx_aMeL = strcmp(temp(:,1),'AME(L)');
    idx_aMeR = strcmp(temp(:,1),'AME(R)');
    idx_MeL = strcmp(temp(:,1),'ME(L)');
    idx_MeR = strcmp(temp(:,1),'ME(R)');
    idx_LoL = strcmp(temp(:,1),'LO(L)');
    idx_LoR = strcmp(temp(:,1),'LO(R)');
    idx_LoPL = strcmp(temp(:,1),'LOP(L)');
    idx_LoPR = strcmp(temp(:,1),'LOP(R)');

    % handle the case where a region is absent for this type
    if any(idx_aMeL)
        val_aMeL = temp{idx_aMeL, 3};
    else
        val_aMeL = 0;
    end
    if any(idx_aMeR)
        val_aMeR = temp{idx_aMeR, 3};
    else
        val_aMeR = 0;
    end
    if any(idx_MeL)
        val_MeL = temp{idx_MeL, 3};
    else
        val_MeL = 0;
    end
    if any(idx_MeR)
        val_MeR = temp{idx_MeR, 3};
    else
        val_MeR = 0;
    end
    if any(idx_LoL)
        val_LoL = temp{idx_LoL, 3};
    else
        val_LoL = 0;
    end
    if any(idx_LoR)
        val_LoR = temp{idx_LoR, 3};
    else
        val_LoR = 0;
    end
    if any(idx_LoPL)
        val_LoPL = temp{idx_LoPL, 3};
    else
        val_LoPL = 0;
    end
    if any(idx_LoPR)
        val_LoPR = temp{idx_LoPR, 3};
    else
        val_LoPR = 0;
    end

    RightFBP_OutNeuropils{i,2} = val_aMeL + val_aMeR;
    RightFBP_OutNeuropils{i,3} = val_MeL + val_MeR;
    RightFBP_OutNeuropils{i,4} = val_LoL + val_LoR;
    RightFBP_OutNeuropils{i,5} = val_LoPL + val_LoPR;
end

% The per-type percentages in RightFBP_OutNeuropils (columns: type, AME, ME, LO,
% LOP) were copied by hand into MCNS/Processed_Data/FBP_output_neuropils.xlsx,
% which was then used to classify each FBP type's target optic lobe (Me / Lo /
% Lop / Multi), reused by the downstream MCNS Figure 4 code.


%% Local function
function [WantSee] = seeConnection_by_region(WantSeetypes,FAFBConnections,FAFB_consolidated_cell_types)
% For each requested cell type, collect its input/output connections and
% summarize them by partner cell type and by neuropil (total synapse counts).

WantSee=table(WantSeetypes,'VariableNames',{'type'});

for i=1:1:size(WantSee,1)
    WantrootidsConsoltype=FAFB_consolidated_cell_types.root_id(strcmpi(FAFB_consolidated_cell_types.primary_type,WantSee.type{i}));
    Wantrootids=unique(WantrootidsConsoltype);
    WantSee.root_ids{i}=Wantrootids;
    if isempty(WantrootidsConsoltype)
        WantSee.root_ids{i}=sscanf(WantSee.type{i},'%ld');
    end

end

for i=1:1:size(WantSee,1)
    Want_root_ids=WantSee.root_ids{i};
    In_Connections=FAFBConnections(ismember(FAFBConnections.post_root_id,Want_root_ids),:);
    Out_Connections=FAFBConnections(ismember(FAFBConnections.pre_root_id,Want_root_ids),:);
    for j=1:1:size(In_Connections,1)
        In_root_ids=In_Connections.pre_root_id(j);
        idx_Consol=find(FAFB_consolidated_cell_types.root_id==In_root_ids);
        if ~isempty(idx_Consol)
            consolType=FAFB_consolidated_cell_types.primary_type{idx_Consol};

        else
            consolType='';

        end
        In_Connections.consoltype{j}=consolType;
        if ~isempty(consolType)
            In_Connections.type{j}=consolType;
        else
            In_Connections.type{j}=num2str(In_root_ids);
        end

    end

    for j=1:1:size(Out_Connections,1)
        Out_root_ids=Out_Connections.post_root_id(j);

        idx_Consol=find(FAFB_consolidated_cell_types.root_id==Out_root_ids);
        if ~isempty(idx_Consol)
            consolType=FAFB_consolidated_cell_types.primary_type{idx_Consol};

        else
            consolType='';

        end

        if ~isempty(consolType)
            Out_Connections.type{j}=consolType;
        else
            Out_Connections.type{j}=num2str(Out_root_ids);
        end
    end

    if isempty(In_Connections)
        WantSee.InNeuronTypes{i}={};
    else
        [uniqueInTypes,~,ic_in]=unique(In_Connections.type);

        for j=1:1:size(uniqueInTypes,1)
            uniqueInTypes{j,2}=sum(In_Connections.syn_count(ic_in==j));
            uniqueInTypes{j,3}=unique(In_Connections.pre_root_id(ic_in==j));
            temp=In_Connections(ic_in==j,:);
            temp=sortrows(temp,"syn_count","descend");
            temp=sortrows(temp,"post_root_id","ascend");
            temp=sortrows(temp,"pre_root_id","ascend");
            uniqueInTypes{j,4}=temp;
        end
        for j=1:1:size(uniqueInTypes,1)
            uniqueInTypes{j,5}=uniqueInTypes{j,2}/sum(cell2mat(uniqueInTypes(:,2)))*100;
        end
        WantSee.InNeuronTypes{i}=sortrows(uniqueInTypes,2,'descend');
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%% OUT
    if isempty(Out_Connections)
        WantSee.OutNeuronTypes{i}={};
    else
        [uniqueOutTypes,~,ic_out]=unique(Out_Connections.type);

        for j=1:1:size(uniqueOutTypes,1)
            uniqueOutTypes{j,2}=sum(Out_Connections.syn_count(ic_out==j));
            uniqueOutTypes{j,3}=unique(Out_Connections.post_root_id(ic_out==j));
            temp=Out_Connections(ic_out==j,:);
            temp=sortrows(temp,"syn_count","descend");
            temp=sortrows(temp,"post_root_id","ascend");
            temp=sortrows(temp,"pre_root_id","ascend");
            uniqueOutTypes{j,4}=temp;
        end
        for j=1:1:size(uniqueOutTypes,1)
            uniqueOutTypes{j,5}=uniqueOutTypes{j,2}/sum(cell2mat(uniqueOutTypes(:,2)))*100;
        end
        WantSee.OutNeuronTypes{i}=sortrows(uniqueOutTypes,2,'descend');
    end

    [unique_in_neuropil,~,ic_in]=unique(In_Connections.neuropil);

    for j=1:1:size(unique_in_neuropil,1)
        idx=ic_in==j;
        unique_in_neuropil{j,2}=sum(In_Connections.syn_count(idx));

    end
    WantSee.InNeuropils{i}=unique_in_neuropil;

    [unique_out_neuropil,~,ic_out]=unique(Out_Connections.neuropil);

    for j=1:1:size(unique_out_neuropil,1)
        idx=ic_out==j;
        unique_out_neuropil{j,2}=sum(Out_Connections.syn_count(idx));

    end
    WantSee.OutNeuropils{i}=unique_out_neuropil;

    if ~isempty(In_Connections)
        In_Connections = sortrows(In_Connections,"syn_count","descend");
        In_Connections = sortrows(In_Connections,"post_root_id","ascend");
        In_Connections = sortrows(In_Connections,"pre_root_id","ascend");
        In_Connections = sortrows(In_Connections,"type","ascend");

        WantSee.InConnnections{i}=In_Connections;

    end

    if ~isempty(Out_Connections)
        Out_Connections = sortrows(Out_Connections,"syn_count","descend");
        Out_Connections = sortrows(Out_Connections,"post_root_id","ascend");
        Out_Connections = sortrows(Out_Connections,"pre_root_id","ascend");
        Out_Connections = sortrows(Out_Connections,"type","ascend");

        WantSee.OutConnections{i}=Out_Connections;

    end

end
end
