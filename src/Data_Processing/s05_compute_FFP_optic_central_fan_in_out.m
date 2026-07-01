%% s05_compute_FFP_optic_central_fan_in_out
% For the right FFP projection neurons and all right-hemisphere optic / central
% neurons, summarize their input/output partners grouped by partner type
% (partner-neuron count, partner-type count, synapse count) and attach the graph
% centrality values from s04. The results are then aggregated per cell type.
% Precursor data processing for Figure 2C, 2D, 2E.
clear all; close all; clc;

baseDir = fileparts(fileparts(mfilename('fullpath')));  % resolved from this script's location (portable; no absolute paths)

% RightFFP_NPIs (from fig_1D_E) and the graph centrality values (from s04)
load(fullfile(baseDir, 'Processed_Data', 'right_neurons_thr0.mat'))
load(fullfile(baseDir, 'Processed_Data', 'allgraph_thr0.mat'))

opt = detectImportOptions(fullfile(baseDir, 'Codex_Data', 'classification.csv'));
opt = setvartype(opt,'root_id','int64');
FAFBClassification = readtable(fullfile(baseDir, 'Codex_Data', 'classification.csv'),opt);

opt = detectImportOptions(fullfile(baseDir, 'Codex_Data', 'consolidated_cell_types.csv'));
opt = setvartype(opt,'root_id','int64');
FAFBTypes = readtable(fullfile(baseDir, 'Codex_Data', 'consolidated_cell_types.csv'),opt);

WantToSee_root_id=RightFFP_NPIs.root_id;
RightFFP_InOut=table(WantToSee_root_id,'VariableNames',{'root_id'});

opt = detectImportOptions(fullfile(baseDir, 'Codex_Data', 'connections_no_threshold.csv'));
opt = setvartype(opt,'pre_root_id','int64');
opt = setvartype(opt,'post_root_id','int64');
FAFBConnections = readtable(fullfile(baseDir, 'Codex_Data', 'connections_no_threshold.csv'),opt);

FAFBConnections(~(ismember(FAFBConnections.pre_root_id,WantToSee_root_id)|ismember(FAFBConnections.post_root_id,WantToSee_root_id)),:)=[];
%%
nRows = size(RightFFP_InOut, 1);
tempType = cell(nRows, 1);
tempIn = cell(nRows, 1);
tempOut = cell(nRows, 1);

parfor i = 1:nRows
    current_root_id = RightFFP_InOut.root_id(i);
    idx_consol = FAFBTypes.root_id == current_root_id;
    % process the connection info
    InConnections = FAFBConnections(FAFBConnections.post_root_id == current_root_id, :);
    % InConnections(~ismember(InConnections.neuropil,{'LA_R', 'ME_R','AME_R', 'LO_R','LOP_R'}),:)=[];
    OutConnections = FAFBConnections(FAFBConnections.pre_root_id == current_root_id, :);

    % update primary_type
    if any(idx_consol)
        tempType{i} = FAFBTypes.primary_type{idx_consol};
    else
        tempType{i} = '';
    end
    % Inputs
    Inputs=unique(InConnections.pre_root_id);
    InputsTypes=cell(size(Inputs));
    for j=1:1:size(Inputs,1)
        idx_consol_Input=FAFBTypes.root_id == Inputs(j,1);

        if any(idx_consol_Input)
            InputsTypes{j,1} = FAFBTypes.primary_type{idx_consol_Input};
        else
            InputsTypes{j,1} = num2str(Inputs(j,1));
        end
    end
    [uniqueInputs,~,ic]= unique(InputsTypes);
    for j=1:1:size(uniqueInputs,1)
        idx=j==ic;
        uniqueInputs{j,2}=Inputs(idx);
        uniqueInputs{j,3}=sum(idx);
        % --- aggregate by partner type ---
        % (1) Find the partner-neuron IDs of this type.
        partner_ids_of_this_type = Inputs(idx);
        
        % (2) Find the connection rows in InConnections matching those partner IDs.
        idx_synapses = ismember(InConnections.pre_root_id, partner_ids_of_this_type);
        
        % (3) Sum only the syn_count of those connections.
        uniqueInputs{j,4} = sum(InConnections.syn_count(idx_synapses));
    end
    if ~isempty(uniqueInputs)
        tempIn{i}=sortrows(uniqueInputs,4,'descend');
    end
    % Outputs
    Outputs=unique(OutConnections.post_root_id);
    OutputsTypes=cell(size(Outputs));
    for j=1:1:size(Outputs,1)
        idx_consol_Output=FAFBTypes.root_id == Outputs(j,1);

        if any(idx_consol_Output)
            OutputsTypes{j,1} = FAFBTypes.primary_type{idx_consol_Output};
        else
            OutputsTypes{j,1} = num2str(Outputs(j,1));
        end
    end
    [uniqueOutputs,~,ic]= unique(OutputsTypes);
    for j=1:1:size(uniqueOutputs,1)
        idx=j==ic;
        uniqueOutputs{j,2}=Outputs(idx);
        uniqueOutputs{j,3}=sum(idx);
        % --- aggregate by partner type ---
        partner_ids_of_this_type = Outputs(idx);
        idx_synapses = ismember(OutConnections.post_root_id, partner_ids_of_this_type);
        uniqueOutputs{j,4} = sum(OutConnections.syn_count(idx_synapses));
    end
    if ~isempty(uniqueOutputs)
        tempOut{i}=sortrows(uniqueOutputs,4,'descend');
    end
end

% merge parfor results
RightFFP_InOut.type = tempType;
RightFFP_InOut.In = tempIn;
RightFFP_InOut.Out = tempOut;
%% optic...
Optic_root_id=FAFBClassification.root_id(strcmp(FAFBClassification.super_class,'optic')&~strcmp(FAFBClassification.side,'left'));

Optic_InOut=table(Optic_root_id,'VariableNames',{'root_id'});

opt = detectImportOptions(fullfile(baseDir, 'Codex_Data', 'connections_no_threshold.csv'));
opt = setvartype(opt,'pre_root_id','int64');
opt = setvartype(opt,'post_root_id','int64');
FAFBConnections = readtable(fullfile(baseDir, 'Codex_Data', 'connections_no_threshold.csv'),opt);

FAFBConnections(~(ismember(FAFBConnections.pre_root_id,Optic_root_id)|ismember(FAFBConnections.post_root_id,Optic_root_id)),:)=[];
%%
nRows = size(Optic_InOut, 1);
tempType = cell(nRows, 1);
tempIn = cell(nRows, 1);
tempOut = cell(nRows, 1);

parfor i = 1:nRows
    current_root_id = Optic_InOut.root_id(i);
    idx_consol = FAFBTypes.root_id == current_root_id;
    % process the connection info
    InConnections = FAFBConnections(FAFBConnections.post_root_id == current_root_id, :);
    % InConnections(~ismember(InConnections.neuropil,{'LA_R', 'ME_R','AME_R', 'LO_R','LOP_R'}),:)=[];
    OutConnections = FAFBConnections(FAFBConnections.pre_root_id == current_root_id, :);

    % update primary_type
    if any(idx_consol)
        tempType{i} = FAFBTypes.primary_type{idx_consol};
    else
        tempType{i} = '';
    end
    % Inputs
    Inputs=unique(InConnections.pre_root_id);
    InputsTypes=cell(size(Inputs));
    for j=1:1:size(Inputs,1)
        idx_consol_Input=FAFBTypes.root_id == Inputs(j,1);

        if any(idx_consol_Input)
            InputsTypes{j,1} = FAFBTypes.primary_type{idx_consol_Input};
        else
            InputsTypes{j,1} = num2str(Inputs(j,1));
        end
    end
    [uniqueInputs,~,ic]= unique(InputsTypes);
    for j=1:1:size(uniqueInputs,1)
        idx=j==ic;
        uniqueInputs{j,2}=Inputs(idx);
        uniqueInputs{j,3}=sum(idx);
        % --- aggregate by partner type ---
        % (1) Find the partner-neuron IDs of this type.
        partner_ids_of_this_type = Inputs(idx);
        
        % (2) Find the connection rows in InConnections matching those partner IDs.
        idx_synapses = ismember(InConnections.pre_root_id, partner_ids_of_this_type);
        
        % (3) Sum only the syn_count of those connections.
        uniqueInputs{j,4} = sum(InConnections.syn_count(idx_synapses));
    end
    if ~isempty(uniqueInputs)
        tempIn{i}=sortrows(uniqueInputs,4,'descend');
    end
    % Outputs
    Outputs=unique(OutConnections.post_root_id);
    OutputsTypes=cell(size(Outputs));
    for j=1:1:size(Outputs,1)
        idx_consol_Output=FAFBTypes.root_id == Outputs(j,1);

        if any(idx_consol_Output)
            OutputsTypes{j,1} = FAFBTypes.primary_type{idx_consol_Output};
        else
            OutputsTypes{j,1} = num2str(Outputs(j,1));
        end
    end
    [uniqueOutputs,~,ic]= unique(OutputsTypes);
    for j=1:1:size(uniqueOutputs,1)
        idx=j==ic;
        uniqueOutputs{j,2}=Outputs(idx);
        uniqueOutputs{j,3}=sum(idx);
        % --- aggregate by partner type ---
        partner_ids_of_this_type = Outputs(idx);
        idx_synapses = ismember(OutConnections.post_root_id, partner_ids_of_this_type);
        uniqueOutputs{j,4} = sum(OutConnections.syn_count(idx_synapses));
    end
    if ~isempty(uniqueOutputs)
        tempOut{i}=sortrows(uniqueOutputs,4,'descend');
    end
end

% merge parfor results
Optic_InOut.type = tempType;
Optic_InOut.In = tempIn;
Optic_InOut.Out = tempOut;

%% central

Central_root_id=FAFBClassification.root_id(strcmp(FAFBClassification.super_class,'central')&~strcmp(FAFBClassification.side,'left'));

Central_InOut=table(Central_root_id,'VariableNames',{'root_id'});

opt = detectImportOptions(fullfile(baseDir, 'Codex_Data', 'connections_no_threshold.csv'));
opt = setvartype(opt,'pre_root_id','int64');
opt = setvartype(opt,'post_root_id','int64');
FAFBConnections = readtable(fullfile(baseDir, 'Codex_Data', 'connections_no_threshold.csv'),opt);

FAFBConnections(~(ismember(FAFBConnections.pre_root_id,Central_root_id)|ismember(FAFBConnections.post_root_id,Central_root_id)),:)=[];
%%
nRows = size(Central_InOut, 1);
tempType = cell(nRows, 1);
tempIn = cell(nRows, 1);
tempOut = cell(nRows, 1);

parfor i = 1:nRows
    current_root_id = Central_InOut.root_id(i);
    idx_consol = FAFBTypes.root_id == current_root_id;
    % process the connection info
    InConnections = FAFBConnections(FAFBConnections.post_root_id == current_root_id, :);
    % InConnections(~ismember(InConnections.neuropil,{'LA_R', 'ME_R','AME_R', 'LO_R','LOP_R'}),:)=[];
    OutConnections = FAFBConnections(FAFBConnections.pre_root_id == current_root_id, :);

    % update primary_type
    if any(idx_consol)
        tempType{i} = FAFBTypes.primary_type{idx_consol};
    else
        tempType{i} = '';
    end
    % Inputs
    Inputs=unique(InConnections.pre_root_id);
    InputsTypes=cell(size(Inputs));
    for j=1:1:size(Inputs,1)
        idx_consol_Input=FAFBTypes.root_id == Inputs(j,1);

        if any(idx_consol_Input)
            InputsTypes{j,1} = FAFBTypes.primary_type{idx_consol_Input};
        else
            InputsTypes{j,1} = num2str(Inputs(j,1));
        end
    end
    [uniqueInputs,~,ic]= unique(InputsTypes);
    for j=1:1:size(uniqueInputs,1)
        idx=j==ic;
        uniqueInputs{j,2}=Inputs(idx);
        uniqueInputs{j,3}=sum(idx);
        % --- aggregate by partner type ---
        % (1) Find the partner-neuron IDs of this type.
        partner_ids_of_this_type = Inputs(idx);
        
        % (2) Find the connection rows in InConnections matching those partner IDs.
        idx_synapses = ismember(InConnections.pre_root_id, partner_ids_of_this_type);
        
        % (3) Sum only the syn_count of those connections.
        uniqueInputs{j,4} = sum(InConnections.syn_count(idx_synapses));
    end
    if ~isempty(uniqueInputs)
        tempIn{i}=sortrows(uniqueInputs,4,'descend');
    end
    % Outputs
    Outputs=unique(OutConnections.post_root_id);
    OutputsTypes=cell(size(Outputs));
    for j=1:1:size(Outputs,1)
        idx_consol_Output=FAFBTypes.root_id == Outputs(j,1);

        if any(idx_consol_Output)
            OutputsTypes{j,1} = FAFBTypes.primary_type{idx_consol_Output};
        else
            OutputsTypes{j,1} = num2str(Outputs(j,1));
        end
    end
    [uniqueOutputs,~,ic]= unique(OutputsTypes);
    for j=1:1:size(uniqueOutputs,1)
        idx=j==ic;
        uniqueOutputs{j,2}=Outputs(idx);
        uniqueOutputs{j,3}=sum(idx);
        % --- aggregate by partner type ---
        partner_ids_of_this_type = Outputs(idx);
        idx_synapses = ismember(OutConnections.post_root_id, partner_ids_of_this_type);
        uniqueOutputs{j,4} = sum(OutConnections.syn_count(idx_synapses));
    end
    if ~isempty(uniqueOutputs)
        tempOut{i}=sortrows(uniqueOutputs,4,'descend');
    end
end

% merge parfor results
Central_InOut.type = tempType;
Central_InOut.In = tempIn;
Central_InOut.Out = tempOut;

%%
for i=1:1:size(RightFFP_InOut,1)
        inStats=RightFFP_InOut.In{i};
        outStats=RightFFP_InOut.Out{i};
        if ~isempty(inStats)
            RightFFP_InOut.InNeuronNumber{i}=sum(cell2mat(inStats(:,3)));
            RightFFP_InOut.InNeuronTypeNumber{i}=size(inStats,1);
            RightFFP_InOut.InNeuronRatio{i}=RightFFP_InOut.InNeuronNumber{i}/RightFFP_InOut.InNeuronTypeNumber{i};

        end
        if ~isempty(outStats)
            RightFFP_InOut.OutNeuronNumber{i}=sum(cell2mat(outStats(:,3)));
            RightFFP_InOut.OutNeuronTypeNumber{i}=size(outStats,1);
            RightFFP_InOut.OutNeuronRatio{i}=RightFFP_InOut.OutNeuronNumber{i}/RightFFP_InOut.OutNeuronTypeNumber{i};

        end
    target_root_id=RightFFP_InOut.root_id(i);
    idx = find(rootIds == target_root_id);
    RightFFP_InOut.betweenness_unweighted(i)=betweenness_vals_unweighted(idx);
    RightFFP_InOut.betweenness_weighted(i)=betweenness_vals_weighted(idx);
    RightFFP_InOut.pagerank(i)=pagerank_vals(idx);


end

for i=1:1:size(Central_InOut,1)
        inStats=Central_InOut.In{i};
        outStats=Central_InOut.Out{i};
        if ~isempty(inStats)
            Central_InOut.InNeuronNumber{i}=sum(cell2mat(inStats(:,3)));
            Central_InOut.InNeuronTypeNumber{i}=size(inStats,1);
            Central_InOut.InNeuronRatio{i}=Central_InOut.InNeuronNumber{i}/Central_InOut.InNeuronTypeNumber{i};
        end
        if ~isempty(outStats)
            Central_InOut.OutNeuronNumber{i}=sum(cell2mat(outStats(:,3)));
            Central_InOut.OutNeuronTypeNumber{i}=size(outStats,1);
            Central_InOut.OutNeuronRatio{i}=Central_InOut.OutNeuronNumber{i}/Central_InOut.OutNeuronTypeNumber{i};
        end

    target_root_id=Central_InOut.root_id(i);
    idx = find(rootIds == target_root_id);
    if ~isempty(idx)
        Central_InOut.betweenness_unweighted(i)=betweenness_vals_unweighted(idx);
        Central_InOut.betweenness_weighted(i)=betweenness_vals_weighted(idx);
        Central_InOut.pagerank(i)=pagerank_vals(idx);
    end
end

for i=1:1:size(Optic_InOut,1)
        inStats=Optic_InOut.In{i};
        outStats=Optic_InOut.Out{i};
        if ~isempty(inStats)
            Optic_InOut.InNeuronNumber{i}=sum(cell2mat(inStats(:,3)));
            Optic_InOut.InNeuronTypeNumber{i}=size(inStats,1);
            Optic_InOut.InNeuronRatio{i}=Optic_InOut.InNeuronNumber{i}/Optic_InOut.InNeuronTypeNumber{i};

        end
        if ~isempty(outStats)
            Optic_InOut.OutNeuronNumber{i}=sum(cell2mat(outStats(:,3)));
            Optic_InOut.OutNeuronTypeNumber{i}=size(outStats,1);
            Optic_InOut.OutNeuronRatio{i}=Optic_InOut.OutNeuronNumber{i}/Optic_InOut.OutNeuronTypeNumber{i};

        end
    target_root_id=Optic_InOut.root_id(i);
    idx = find(rootIds == target_root_id);
    if ~isempty(idx)
        Optic_InOut.betweenness_unweighted(i)=betweenness_vals_unweighted(idx);
        Optic_InOut.betweenness_weighted(i)=betweenness_vals_weighted(idx);
        Optic_InOut.pagerank(i)=pagerank_vals(idx);
    end
end
%
[type_RightFFP_InOut,~,ic]=unique(RightFFP_InOut.type);
type_RightFFP_InOut=table(type_RightFFP_InOut,'VariableNames',{'type'});
for i=1:1:size(type_RightFFP_InOut,1)
    idx=ic==i;
    type_RightFFP_InOut.InNeuronNumber{i}=mean(cell2mat(RightFFP_InOut.InNeuronNumber(idx)),'omitmissing');
    type_RightFFP_InOut.InNeuronTypeNumber{i}=mean(cell2mat(RightFFP_InOut.InNeuronTypeNumber(idx)),'omitmissing');
    type_RightFFP_InOut.InNeuronRatio{i}=mean(cell2mat(RightFFP_InOut.InNeuronRatio(idx)),'omitmissing');

    type_RightFFP_InOut.OutNeuronNumber{i}=mean(cell2mat(RightFFP_InOut.OutNeuronNumber(idx)),'omitmissing');
    type_RightFFP_InOut.OutNeuronTypeNumber{i}=mean(cell2mat(RightFFP_InOut.OutNeuronTypeNumber(idx)),'omitmissing');
    type_RightFFP_InOut.OutNeuronRatio{i}=mean(cell2mat(RightFFP_InOut.OutNeuronRatio(idx)),'omitmissing');

    type_RightFFP_InOut.betweenness_unweighted{i}=mean((RightFFP_InOut.betweenness_unweighted(idx)),'omitmissing');
    type_RightFFP_InOut.betweenness_weighted{i}=mean((RightFFP_InOut.betweenness_weighted(idx)),'omitmissing');
    type_RightFFP_InOut.pagerank{i}=mean((RightFFP_InOut.pagerank(idx)),'omitmissing');

end

[type_Central_InOut,~,ic]=unique(Central_InOut.type);
type_Central_InOut=table(type_Central_InOut,'VariableNames',{'type'});
for i=1:1:size(type_Central_InOut,1)
    idx=ic==i;
    type_Central_InOut.InNeuronNumber{i}=mean(cell2mat(Central_InOut.InNeuronNumber(idx)),'omitmissing');
    type_Central_InOut.InNeuronTypeNumber{i}=mean(cell2mat(Central_InOut.InNeuronTypeNumber(idx)),'omitmissing');
    type_Central_InOut.InNeuronRatio{i}=mean(cell2mat(Central_InOut.InNeuronRatio(idx)),'omitmissing');

    type_Central_InOut.OutNeuronNumber{i}=mean(cell2mat(Central_InOut.OutNeuronNumber(idx)),'omitmissing');
    type_Central_InOut.OutNeuronTypeNumber{i}=mean(cell2mat(Central_InOut.OutNeuronTypeNumber(idx)),'omitmissing');
    type_Central_InOut.OutNeuronRatio{i}=mean(cell2mat(Central_InOut.OutNeuronRatio(idx)),'omitmissing');

    type_Central_InOut.betweenness_unweighted{i}=mean((Central_InOut.betweenness_unweighted(idx)),'omitmissing');
    type_Central_InOut.betweenness_weighted{i}=mean((Central_InOut.betweenness_weighted(idx)),'omitmissing');
    type_Central_InOut.pagerank{i}=mean((Central_InOut.pagerank(idx)),'omitmissing');

end

[type_Optic_InOut,~,ic]=unique(Optic_InOut.type);
type_Optic_InOut=table(type_Optic_InOut,'VariableNames',{'type'});
for i=1:1:size(type_Optic_InOut,1)
    idx=ic==i;
    type_Optic_InOut.InNeuronNumber{i}=mean(cell2mat(Optic_InOut.InNeuronNumber(idx)),'omitmissing');
    type_Optic_InOut.InNeuronTypeNumber{i}=mean(cell2mat(Optic_InOut.InNeuronTypeNumber(idx)),'omitmissing');
    type_Optic_InOut.InNeuronRatio{i}=mean(cell2mat(Optic_InOut.InNeuronRatio(idx)),'omitmissing');

    type_Optic_InOut.OutNeuronNumber{i}=mean(cell2mat(Optic_InOut.OutNeuronNumber(idx)),'omitmissing');
    type_Optic_InOut.OutNeuronTypeNumber{i}=mean(cell2mat(Optic_InOut.OutNeuronTypeNumber(idx)),'omitmissing');
    type_Optic_InOut.OutNeuronRatio{i}=mean(cell2mat(Optic_InOut.OutNeuronRatio(idx)),'omitmissing');

    type_Optic_InOut.betweenness_unweighted{i}=mean((Optic_InOut.betweenness_unweighted(idx)),'omitmissing');
    type_Optic_InOut.betweenness_weighted{i}=mean((Optic_InOut.betweenness_weighted(idx)),'omitmissing');
    type_Optic_InOut.pagerank{i}=mean((Optic_InOut.pagerank(idx)),'omitmissing');

end

%% Save the per-type fan-in / fan-out summaries
save(fullfile(baseDir, 'Processed_Data', 'FFP_optic_central_fan_in_out.mat'), ...
    'type_Optic_InOut', 'type_Central_InOut', 'type_RightFFP_InOut');