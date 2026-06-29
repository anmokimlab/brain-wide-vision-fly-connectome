clear all; close all; clc;

baseDir = fileparts(fileparts(mfilename('fullpath')));  % resolved from this script's location (portable; no absolute paths)

% BLP neuron classification (provides BLP_R_NPIs), saved by
% Figures/fig_6B_postPI_prePI_BLP_neurons.m
load(fullfile(baseDir, 'Processed_Data', 'BLP_neurons_thr0.mat'))

%% Load Codex data
opt = detectImportOptions(fullfile(baseDir, 'Codex_Data', 'neurons.csv'));
opt = setvartype(opt,'root_id','int64');
FAFBNeurons = readtable(fullfile(baseDir, 'Codex_Data', 'neurons.csv'),opt);

opt = detectImportOptions(fullfile(baseDir, 'Codex_Data', 'connections_no_threshold.csv'));
opt = setvartype(opt,'pre_root_id','int64');
opt = setvartype(opt,'post_root_id','int64');
FAFBConnections = readtable(fullfile(baseDir, 'Codex_Data', 'connections_no_threshold.csv'),opt);
FAFBConnections(FAFBConnections.syn_count<5,:)=[];   % keep only connections with >= 5 synapses

opt = detectImportOptions(fullfile(baseDir, 'Codex_Data', 'consolidated_cell_types.csv'));
opt = setvartype(opt,'root_id','int64');
FAFBTypes = readtable(fullfile(baseDir, 'Codex_Data', 'consolidated_cell_types.csv'),opt);

NT_Types=unique(FAFBNeurons.nt_type);
NT_Types(1,:)=[];
NT_Types{end+1,1}='Unidentified';

FAFBNeuropils=unique(FAFBConnections.neuropil);
FAFBNeuropil_OpticLobeRight={'LA_R';'AME_R';'ME_R';'LO_R';'LOP_R'};
FAFBNeuropil_OpticLobeLeft={'LA_L';'AME_L';'ME_L';'LO_L';'LOP_L'};
FAFBNeuropil_Central=FAFBNeuropils(~ismember(FAFBNeuropils,FAFBNeuropil_OpticLobeRight));
FAFBNeuropil_Central=FAFBNeuropil_Central(~ismember(FAFBNeuropil_Central,FAFBNeuropil_OpticLobeLeft));
FAFBNeuropil_Central=FAFBNeuropil_Central(~ismember(FAFBNeuropil_Central,'UNASGD'));

%% Assign each BLP_R neuron its neurotransmitter type
for i=1:1:size(BLP_R_NPIs,1)
    idx=FAFBNeurons.root_id==BLP_R_NPIs.root_id(i);
    if ~any(idx)
        BLP_R_NPIs.nt_type{i}='Unidentified';
    elseif isempty(FAFBNeurons.nt_type{idx})
        BLP_R_NPIs.nt_type{i}='Unidentified';
    else
        BLP_R_NPIs.nt_type{i}=FAFBNeurons.nt_type{idx};
    end
end

for i=1:1:size(NT_Types,1)
    idx=strcmp(BLP_R_NPIs.nt_type,NT_Types{i,1});
    NT_Types{i,2}=sum(idx);
end

%% Per-neuropil input -> output connection table (synapse-weighted)
BLP_Connection=cell(length(FAFBNeuropils),length(FAFBNeuropils));

for i=1:1:size(BLP_R_NPIs,1)
    current_root_id=BLP_R_NPIs.root_id(i);
    PostConnections=FAFBConnections(FAFBConnections.post_root_id==current_root_id,:);
    PreConnections=FAFBConnections(FAFBConnections.pre_root_id==current_root_id,:);
    [OutNeuropils, ~, ic]=unique(PreConnections.neuropil);
    for j=1:1:size(OutNeuropils,1)
        idx=ic==j;
        OutNeuropils{j,2}=sum(PreConnections.syn_count(idx));
    end
    [InNeuropils, ~, ic]=unique(PostConnections.neuropil);
    for j=1:1:size(InNeuropils,1)
        idx=ic==j;
        InNeuropils{j,2}=sum(PostConnections.syn_count(idx));
    end

    for j=1:1:size(InNeuropils,1)
        idx_In=find(strcmp(FAFBNeuropils,InNeuropils{j,1}));
        syn_ratio_In=InNeuropils{j,2}/sum(cell2mat(InNeuropils(:,2)));
        for k=1:1:size(OutNeuropils,1)
            idx_Out=find(strcmp(FAFBNeuropils,OutNeuropils{k,1}));

            temp=BLP_Connection{idx_In,idx_Out};
            temp{end+1,1}=current_root_id;
            temp{end,2}=OutNeuropils{k,2}*syn_ratio_In;
            temp{end,3}=BLP_R_NPIs.nt_type{i};
            BLP_Connection{idx_In,idx_Out}=temp;
        end
    end
end

%% Merge AME into ME (right and left), assuming root_ids are unique
aMe_R_idx=find(strcmp(FAFBNeuropils,'AME_R'));
Me_R_idx=find(strcmp(FAFBNeuropils,'ME_R'));
aMe_L_idx=find(strcmp(FAFBNeuropils,'AME_L'));
Me_L_idx=find(strcmp(FAFBNeuropils,'ME_L'));

for i=1:1:size(BLP_Connection,2)
    Me_R_temp=BLP_Connection{Me_R_idx,i};
    AMe_R_Temp=BLP_Connection{aMe_R_idx,i};
    if isempty(Me_R_temp)||isempty(AMe_R_Temp)
        BLP_Connection{Me_R_idx,i}=[Me_R_temp;AMe_R_Temp];
    else
        [MergedAMeIdx,MergedMeIdx]=(ismember(cell2mat(AMe_R_Temp(:,1)),cell2mat(Me_R_temp(:,1))));
        MergedAMeIdx=find(MergedAMeIdx);
        MergedMeIdx(MergedMeIdx==0)=[];
        if ~isempty(MergedMeIdx)
            for j=1:1:size(MergedMeIdx,1)
                Me_R_temp{MergedMeIdx(j),2}= Me_R_temp{MergedMeIdx(j),2}+AMe_R_Temp{MergedAMeIdx(j),2};
            end
        end
        % Append AMe_R-only neurons (not present in Me_R)
        NonMergedAMeIdx = ~ismember(cell2mat(AMe_R_Temp(:,1)), cell2mat(Me_R_temp(:,1)));
        Me_R_temp = [Me_R_temp; AMe_R_Temp(NonMergedAMeIdx, :)];
        BLP_Connection{Me_R_idx,i}=Me_R_temp;
    end
end

for i=1:1:size(BLP_Connection,1)
    Me_L_temp=BLP_Connection{i,Me_L_idx};
    AMe_L_Temp=BLP_Connection{i,aMe_L_idx};
    if isempty(Me_L_temp)||isempty(AMe_L_Temp)
        BLP_Connection{i,Me_L_idx}=[Me_L_temp;AMe_L_Temp];
    else
        [MergedAMeIdx,MergedMeIdx]=(ismember(cell2mat(AMe_L_Temp(:,1)),cell2mat(Me_L_temp(:,1))));
        MergedAMeIdx=find(MergedAMeIdx);
        MergedMeIdx(MergedMeIdx==0)=[];
        if ~isempty(MergedMeIdx)
            for j=1:1:size(MergedMeIdx,1)
                Me_L_temp{MergedMeIdx(j),2}= Me_L_temp{MergedMeIdx(j),2}+AMe_L_Temp{MergedAMeIdx(j),2};
            end
        end
        % Append AMe_L-only neurons (not present in Me_L)
        NonMergedAMeIdx = ~ismember(cell2mat(AMe_L_Temp(:,1)), cell2mat(Me_L_temp(:,1)));
        Me_L_temp = [Me_L_temp; AMe_L_Temp(NonMergedAMeIdx, :)];
        BLP_Connection{i,Me_L_idx}=Me_L_temp;
    end
end

%% Right optic lobe -> left optic lobe connection counts
% Indices into FAFBNeuropils: Me R 55, Lo R 45, LoP R 43 / Me L 54, Lo L 44, LoP L 42
idx_RightOpticLobes=[55 45 43];
idx_LeftOpticLobes=[54 44 42];

BLP_Connection_Matrix=zeros(length(idx_RightOpticLobes),length(idx_LeftOpticLobes));
for i=1:1:length(idx_RightOpticLobes)
    for j=1:1:length(idx_LeftOpticLobes)
        currentConnection=BLP_Connection{idx_RightOpticLobes(i),idx_LeftOpticLobes(j)};
        if isempty(currentConnection)
            BLP_Connection_Matrix(i,j)=0;
        else
            BLP_Connection_Matrix(i,j)=size(currentConnection,1);
        end
    end
end

%% Collect the connecting neuron names / types per neuropil pair
BLP_Connection_NeuronsNames=cell(length(idx_RightOpticLobes),length(idx_LeftOpticLobes));
BLP_Connection_NeuronsNamesType=cell(length(idx_RightOpticLobes),length(idx_LeftOpticLobes));

for i=1:1:length(idx_RightOpticLobes)
    for j=1:1:length(idx_LeftOpticLobes)
        currentConnection=BLP_Connection{idx_RightOpticLobes(i),idx_LeftOpticLobes(j)};
        if isempty(currentConnection)
            continue;
        else
            for k=1:1:size(currentConnection,1)
                idx_Types=FAFBTypes.root_id==currentConnection{k,1};
                currentConnection{k,4}=FAFBTypes.primary_type{idx_Types};
            end
            [uniquecurrentConnectionType,~,ic] = unique(currentConnection(:,4));
            a_counts = accumarray(ic,1);
            uniquecurrentConnectionType(:,2)=num2cell(a_counts);
            for k=1:1:size(uniquecurrentConnectionType,1)
                idx=ic==k;
                uniquecurrentConnectionType{k,3}=sum(cell2mat(currentConnection(idx,2)));
            end
            uniquecurrentConnectionType=sortrows(uniquecurrentConnectionType,3,'descend');
            BLP_Connection_NeuronsNames{i,j}=currentConnection;
            BLP_Connection_NeuronsNamesType{i,j}=uniquecurrentConnectionType;
        end
    end
end

%% Total neuron count per neuropil (added as an extra row/column)
ME_R_NumberNeurons=[BLP_Connection_NeuronsNames{1, 1};BLP_Connection_NeuronsNames{1, 2};BLP_Connection_NeuronsNames{1, 3}];
ME_R_NumberNeurons=size(unique(cell2mat(ME_R_NumberNeurons(:,1))),1);
LO_R_NumberNeurons=[BLP_Connection_NeuronsNames{2, 1};BLP_Connection_NeuronsNames{2, 2};BLP_Connection_NeuronsNames{2, 3}];
LO_R_NumberNeurons=size(unique(cell2mat(LO_R_NumberNeurons(:,1))),1);
LoP_R_NumberNeurons=[BLP_Connection_NeuronsNames{3, 1};BLP_Connection_NeuronsNames{3, 2};BLP_Connection_NeuronsNames{3, 3}];
LoP_R_NumberNeurons=size(unique(cell2mat(LoP_R_NumberNeurons(:,1))),1);

ME_L_NumberNeurons=[BLP_Connection_NeuronsNames{1, 1};BLP_Connection_NeuronsNames{2, 1};BLP_Connection_NeuronsNames{3, 1}];
ME_L_NumberNeurons=size(unique(cell2mat(ME_L_NumberNeurons(:,1))),1);
LO_L_NumberNeurons=[BLP_Connection_NeuronsNames{1, 2};BLP_Connection_NeuronsNames{2, 2};BLP_Connection_NeuronsNames{3, 2}];
LO_L_NumberNeurons=size(unique(cell2mat(LO_L_NumberNeurons(:,1))),1);
LoP_L_NumberNeurons=[BLP_Connection_NeuronsNames{1, 3};BLP_Connection_NeuronsNames{2, 3};BLP_Connection_NeuronsNames{3, 3}];
LoP_L_NumberNeurons=size(unique(cell2mat(LoP_L_NumberNeurons(:,1))),1);

BLP_Connection_Matrix(4,1)=ME_L_NumberNeurons;
BLP_Connection_Matrix(4,2)=LO_L_NumberNeurons;
BLP_Connection_Matrix(4,3)=LoP_L_NumberNeurons;

BLP_Connection_Matrix(1,4)=ME_R_NumberNeurons;
BLP_Connection_Matrix(2,4)=LO_R_NumberNeurons;
BLP_Connection_Matrix(3,4)=LoP_R_NumberNeurons;

%% Figure 1 (panel 6D): R-to-L optic-lobe bubble chart
% The extra (4th) row and column encode the total neuron count per neuropil.
NeuronN_Thr=0;
figure(1);set(gcf,'Color','w')
[YData, XData]=find(BLP_Connection_Matrix>=NeuronN_Thr);
bubbleSizes=BLP_Connection_Matrix(BLP_Connection_Matrix>=NeuronN_Thr);
bubblechart(XData',YData',bubbleSizes');
bubblelegend('Neuron','Location','eastoutside')
set(gca,'XTick',1:1:3,'XTickLabels',{'Me L','Lo L','LoP L'},...
    'YTick',1:1:3,'YTickLabels',{'Me R','Lo R','LoP R'},'Box','off','TickDir','out')
title('Right Bino Neuron')

%% Per-neuropil neurotransmitter (NT) composition
for i=1:1:size(BLP_Connection_NeuronsNames,1)
    for j=1:1:size(BLP_Connection_NeuronsNames,2)
        if size(BLP_Connection_NeuronsNames{i,j},1)<NeuronN_Thr
            BLP_Connection_NeuronsNames{i,j}=[];
        end
    end
end

NT_Me_R={}; NT_Me_L={};
NT_Lo_R={}; NT_Lo_L={};
NT_LoP_R={}; NT_LoP_L={};

for i=1:1:size(BLP_Connection_NeuronsNames,2)
    NT_Me_R=[NT_Me_R;BLP_Connection_NeuronsNames{1,i}];
    NT_Lo_R=[NT_Lo_R;BLP_Connection_NeuronsNames{2,i}];
    NT_LoP_R=[NT_LoP_R;BLP_Connection_NeuronsNames{3,i}];
end

for i=1:1:size(BLP_Connection_NeuronsNames,1)
    NT_Me_L=[NT_Me_L;BLP_Connection_NeuronsNames{i,1}];
    NT_Lo_L=[NT_Lo_L;BLP_Connection_NeuronsNames{i,2}];
    NT_LoP_L=[NT_LoP_L;BLP_Connection_NeuronsNames{i,3}];
end

NT_Me_R_unique=num2cell(unique(cell2mat(NT_Me_R(:,1))));
for i=1:1:size(NT_Me_R_unique,1)
    idx=find(cell2mat(NT_Me_R(:,1))==NT_Me_R_unique{i,1},1);
    NT_Me_R_unique{i,2}=NT_Me_R{idx,3};
end
[NT_Me_R,~,ic] = unique(NT_Me_R_unique(:,2));
NT_Me_R(:,2)=num2cell(accumarray(ic,1));

NT_Lo_R_unique=num2cell(unique(cell2mat(NT_Lo_R(:,1))));
for i=1:1:size(NT_Lo_R_unique,1)
    idx=find(cell2mat(NT_Lo_R(:,1))==NT_Lo_R_unique{i,1},1);
    NT_Lo_R_unique{i,2}=NT_Lo_R{idx,3};
end
[NT_Lo_R,~,ic] = unique(NT_Lo_R_unique(:,2));
NT_Lo_R(:,2)=num2cell(accumarray(ic,1));

NT_LoP_R_unique=num2cell(unique(cell2mat(NT_LoP_R(:,1))));
for i=1:1:size(NT_LoP_R_unique,1)
    idx=find(cell2mat(NT_LoP_R(:,1))==NT_LoP_R_unique{i,1},1);
    NT_LoP_R_unique{i,2}=NT_LoP_R{idx,3};
end
[NT_LoP_R,~,ic] = unique(NT_LoP_R_unique(:,2));
NT_LoP_R(:,2)=num2cell(accumarray(ic,1));

NT_Me_L_unique=num2cell(unique(cell2mat(NT_Me_L(:,1))));
for i=1:1:size(NT_Me_L_unique,1)
    idx=find(cell2mat(NT_Me_L(:,1))==NT_Me_L_unique{i,1},1);
    NT_Me_L_unique{i,2}=NT_Me_L{idx,3};
end
[NT_Me_L,~,ic] = unique(NT_Me_L_unique(:,2));
NT_Me_L(:,2)=num2cell(accumarray(ic,1));

NT_Lo_L_unique=num2cell(unique(cell2mat(NT_Lo_L(:,1))));
for i=1:1:size(NT_Lo_L_unique,1)
    idx=find(cell2mat(NT_Lo_L(:,1))==NT_Lo_L_unique{i,1},1);
    NT_Lo_L_unique{i,2}=NT_Lo_L{idx,3};
end
[NT_Lo_L,~,ic] = unique(NT_Lo_L_unique(:,2));
NT_Lo_L(:,2)=num2cell(accumarray(ic,1));

NT_LoP_L_unique=num2cell(unique(cell2mat(NT_LoP_L(:,1))));
for i=1:1:size(NT_LoP_L_unique,1)
    idx=find(cell2mat(NT_LoP_L(:,1))==NT_LoP_L_unique{i,1},1);
    NT_LoP_L_unique{i,2}=NT_LoP_L{idx,3};
end
[NT_LoP_L,~,ic] = unique(NT_LoP_L_unique(:,2));
NT_LoP_L(:,2)=num2cell(accumarray(ic,1));

%% Figure 2 (panel 6C): BLP_R cell-type distribution (counts entered manually)
figure(2);set(gcf,'Color','w')
x = ["LC14b" "LC14a1" "LC14a2" "MeMe e02" "MTe07" "MeMe e01" "MeMe e08" "MeMe e07" "l-LNv" "Others"];
y = [18 15 12 9 6 6 6 5 4 25];
bar(x,y)
set(gca,'Box','off','TickDir','out')
