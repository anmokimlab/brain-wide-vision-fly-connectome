clear all; close all; clc;

baseDir = fileparts(fileparts(mfilename('fullpath')));  % resolved from this script's location (portable; no absolute paths)

% FFP / FBP / BDP neuron classification (provides RightFBP_NPIs), saved by
% Figures/fig_1D_E_postPI_prePI_right_neurons.m
load(fullfile(baseDir, 'Processed_Data', 'right_neurons_thr0.mat'))

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

%% Map each FBP neuron to its neurotransmitter (NT) type
for i=1:1:size(RightFBP_NPIs,1)
    idx=FAFBNeurons.root_id==RightFBP_NPIs.root_id(i);
    if ~any(idx)
        RightFBP_NPIs.nt_type{i}='Unidentified';
    elseif isempty(FAFBNeurons.nt_type{idx})
        RightFBP_NPIs.nt_type{i}='Unidentified';
    else
        RightFBP_NPIs.nt_type{i}=FAFBNeurons.nt_type{idx};
    end
end

for i=1:1:size(NT_Types,1)
    idx=strcmp(RightFBP_NPIs.nt_type,NT_Types{i,1});
    NT_Types{i,2}=sum(idx);

end
%%  From optic lobe to higher brain region (total and NT types)
FBP_Central_Connection=cell(length(FAFBNeuropils),length(FAFBNeuropils));

for i=1:1:size(RightFBP_NPIs,1)
    current_root_id=RightFBP_NPIs.root_id(i);
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

            temp=FBP_Central_Connection{idx_In,idx_Out};
            temp{end+1,1}=current_root_id;
            temp{end,2}=OutNeuropils{k,2}*syn_ratio_In;
            temp{end,3}=RightFBP_NPIs.nt_type{i};
            FBP_Central_Connection{idx_In,idx_Out}=temp;

        end
    end
end


%% Merge AME into ME
% assumes RightFBP_NPIs.root_id is unique
aMe_R_idx=find(strcmp(FAFBNeuropils,'AME_R'));
Me_R_idx=find(strcmp(FAFBNeuropils,'ME_R'));
aMe_L_idx=find(strcmp(FAFBNeuropils,'AME_L'));
Me_L_idx=find(strcmp(FAFBNeuropils,'ME_L'));

for i=1:1:size(FBP_Central_Connection,1)
    Me_R_temp=FBP_Central_Connection{i,Me_R_idx};
    AMe_R_Temp=FBP_Central_Connection{i,aMe_R_idx};
    if isempty(Me_R_temp)||isempty(AMe_R_Temp)
        FBP_Central_Connection{i,Me_R_idx}=[Me_R_temp;AMe_R_Temp];
    else
        [MergedAMeIdx,MergedMeIdx]=(ismember(cell2mat(AMe_R_Temp(:,1)),cell2mat(Me_R_temp(:,1))));
        MergedAMeIdx=find(MergedAMeIdx);
        MergedMeIdx(MergedMeIdx==0)=[];
        if ~isempty(MergedMeIdx)
            for j=1:1:size(MergedMeIdx,1)
                Me_R_temp{MergedMeIdx(j),2}= Me_R_temp{MergedMeIdx(j),2}+AMe_R_Temp{MergedAMeIdx(j),2}; % ratios, so summing is fine
            end
        end
        % --- merge neurons unique to AME_R ---
        % use ismember again to find neurons NOT present in ME_R
        NonMergedAMeIdx = ~ismember(cell2mat(AMe_R_Temp(:,1)), cell2mat(Me_R_temp(:,1)));

        % append those unique neurons to the bottom of Me_R_temp
        Me_R_temp = [Me_R_temp; AMe_R_Temp(NonMergedAMeIdx, :)];
        % --- end ---

        % Me_R_temp now contains all merged neurons
        FBP_Central_Connection{i,Me_R_idx}=Me_R_temp;
    end
end
%% Bubble-chart connectivity matrix
%%% indices into FAFBNeuropils: Me R = 55, Lo R = 45, LoP R = 43
idx_RightOpticLobes=[55 45 43];
idx_CentralBrains=[1 2 5:37 40 41 46:53 56:74 76:79];

FBP_Central_Connection_Matrix=zeros(length(idx_CentralBrains),length(idx_RightOpticLobes));

for i=1:1:length(idx_CentralBrains)
    for j=1:1:length(idx_RightOpticLobes)
        currentConnection=FBP_Central_Connection{idx_CentralBrains(i),idx_RightOpticLobes(j)};
        if isempty(currentConnection)
            FBP_Central_Connection_Matrix(i,j)=0;
        else
            FBP_Central_Connection_Matrix(i,j)=size(currentConnection,1);
        end
    end
end

NeuronN_Thr=1;   % minimum neuron count required to draw a bubble

%% sending neuron names...
FBP_Central_Optic_Connection=FBP_Central_Connection(idx_CentralBrains,idx_RightOpticLobes);

unique_neuron_size_matrix_optic_lobe=zeros(size(FBP_Central_Optic_Connection,2),1);
for i=1:1:3
    FBP_Names=[];
    for j=1:1:size(FBP_Central_Optic_Connection,1)
        FBP_Names=[FBP_Names; FBP_Central_Optic_Connection{j, i}];
    end
    if ~isempty(FBP_Names)
        unique_neuron_size_matrix_optic_lobe(i,1)=length(unique(cell2mat(FBP_Names(:,1))));
    end
end

unique_neuron_size_matrix_central_brain=zeros(size(FBP_Central_Optic_Connection,1),1);

for j=1:1:size(FBP_Central_Optic_Connection,1)
    FBP_Names=[];
    for i=1:1:3
        FBP_Names=[FBP_Names; FBP_Central_Optic_Connection{j,i}];
    end
    if ~isempty(FBP_Names)
        unique_neuron_size_matrix_central_brain(j,1)=length(unique(cell2mat(FBP_Names(:,1))));
    end
end
[~,idx_SortedCentralBrains]=sort(unique_neuron_size_matrix_central_brain,'descend');
SortedCentralBrains=FAFBNeuropil_Central(idx_SortedCentralBrains);
WantToSeeN=1:1:10;
FBP_Optic_Central_Connection_sorted=FBP_Central_Optic_Connection(idx_SortedCentralBrains(WantToSeeN),:);

FBP_Central_Connection_NeuronsNames=cell(length(WantToSeeN),length(idx_RightOpticLobes));
FBP_Central_Connection_NeuronsNamesType=cell(length(WantToSeeN),length(idx_RightOpticLobes));

for i=1:1:length(WantToSeeN)
    for j=1:1:length(idx_RightOpticLobes)
        currentConnection=FBP_Optic_Central_Connection_sorted{i,j};
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
            FBP_Central_Connection_NeuronsNames{i,j}=currentConnection;
            FBP_Central_Connection_NeuronsNamesType{i,j}=uniquecurrentConnectionType;
        end
    end
end

%% figure 1: sorted bubble chart
figure(1);set(gcf,'Color','w')

temp_matrix=FBP_Central_Connection_Matrix(idx_SortedCentralBrains(WantToSeeN),:);
% extra column: total unique neurons per central-brain region
temp_matrix(WantToSeeN,4)=unique_neuron_size_matrix_central_brain(idx_SortedCentralBrains(WantToSeeN));

unique_neuron_size_matrix_optic_lobe_wanted=zeros(size(FBP_Optic_Central_Connection_sorted,2),1);
for i=1:1:3
    FBP_Names=[];
    for j=1:1:size(FBP_Optic_Central_Connection_sorted,1)
        FBP_Names=[FBP_Names; FBP_Optic_Central_Connection_sorted{j,i}];
    end
    if ~isempty(FBP_Names)
        unique_neuron_size_matrix_optic_lobe_wanted(i,1)=length(unique(cell2mat(FBP_Names(:,1))));
    end
end

% extra row: total unique neurons per right optic-lobe region
temp_matrix(length(WantToSeeN)+1,1:1:3)=unique_neuron_size_matrix_optic_lobe_wanted';

[YData, XData]=find(temp_matrix>=NeuronN_Thr);
FBP_Central_Connection_Matrix_bubble=temp_matrix(temp_matrix>=NeuronN_Thr);
XData=XData';
YData=YData';
FBP_Central_Connection_Matrix_bubble=FBP_Central_Connection_Matrix_bubble';
bubblechart(XData,YData,FBP_Central_Connection_Matrix_bubble,MarkerFaceAlpha=0.8,MarkerFaceColor='#48494B',MarkerEdgeColor='#48494B',MarkerEdgeAlpha=0.8);
bubblelegend('Neuron','Location','eastoutside')

set(gca,'YTick',1:1:length(WantToSeeN),'YTickLabels',FAFBNeuropil_Central(idx_SortedCentralBrains(WantToSeeN)),...
    'XTick',1:1:size(FBP_Central_Connection_Matrix,2),'XTickLabels',{'Me R','Lo R','LoP R'},'Box','off','TickDir','out')
title('Right Feedback Neuron')
grid on



%% NT aggregation (exported to Excel for the NT pie charts)
NT_Me={};
NT_Lo={};
NT_LoP={};
NT_PLP={};
NT_IPS={};
NT_SPS={};
NT_WED={};
NT_PVLP={};
NT_GNG={};
NT_SAD={};
NT_SPSL={};
NT_IB={};
NT_ICL={};


for i=1:1:size(FBP_Central_Connection_NeuronsNames,1)
    for j=1:1:size(FBP_Central_Connection_NeuronsNames,2)
        if size(FBP_Central_Connection_NeuronsNames{i,j},1)<NeuronN_Thr
            FBP_Central_Connection_NeuronsNames{i,j}=[];
        end
    end
end

for i=1:1:size(FBP_Central_Connection_NeuronsNames,1)
    NT_Me=[NT_Me;FBP_Central_Connection_NeuronsNames{i,1}];
    NT_Lo=[NT_Lo;FBP_Central_Connection_NeuronsNames{i,2}];
    NT_LoP=[NT_LoP;FBP_Central_Connection_NeuronsNames{i,3}];
end

for i=1:1:size(FBP_Central_Connection_NeuronsNames,2)

    NT_SPS=[NT_SPS;FBP_Central_Connection_NeuronsNames{1,i}];
    NT_IPS=[NT_IPS;FBP_Central_Connection_NeuronsNames{2,i}];
    NT_PLP=[NT_PLP;FBP_Central_Connection_NeuronsNames{3,i}];
    NT_GNG=[NT_GNG;FBP_Central_Connection_NeuronsNames{4,i}];
    NT_SPSL=[NT_SPSL;FBP_Central_Connection_NeuronsNames{5,i}];
    NT_WED=[NT_WED;FBP_Central_Connection_NeuronsNames{6,i}];
    NT_SAD=[NT_SAD;FBP_Central_Connection_NeuronsNames{7,i}];
    NT_ICL=[NT_ICL;FBP_Central_Connection_NeuronsNames{8,i}];
    NT_IB=[NT_IB;FBP_Central_Connection_NeuronsNames{9,i}];
    NT_PVLP=[NT_PVLP;FBP_Central_Connection_NeuronsNames{10,i}];


end
%%
NT_Me_unique = num2cell(unique(cell2mat(NT_Me(:,1))));
for i = 1:size(NT_Me_unique,1)
    idx = find(cell2mat(NT_Me(:,1)) == NT_Me_unique{i,1}, 1);
    NT_Me_unique{i,2} = NT_Me{idx,3};
end
[NT_Me,~,ic] = unique(NT_Me_unique(:,2));
a_counts = accumarray(ic,1);
NT_Me(:,2) = num2cell(a_counts);

% Lo
NT_Lo_unique = num2cell(unique(cell2mat(NT_Lo(:,1))));
for i = 1:size(NT_Lo_unique,1)
    idx = find(cell2mat(NT_Lo(:,1)) == NT_Lo_unique{i,1}, 1);
    NT_Lo_unique{i,2} = NT_Lo{idx,3};
end
[NT_Lo,~,ic] = unique(NT_Lo_unique(:,2));
a_counts = accumarray(ic,1);
NT_Lo(:,2) = num2cell(a_counts);

% LoP
NT_LoP_unique = num2cell(unique(cell2mat(NT_LoP(:,1))));
for i = 1:size(NT_LoP_unique,1)
    idx = find(cell2mat(NT_LoP(:,1)) == NT_LoP_unique{i,1}, 1);
    NT_LoP_unique{i,2} = NT_LoP{idx,3};
end
[NT_LoP,~,ic] = unique(NT_LoP_unique(:,2));
a_counts = accumarray(ic,1);
NT_LoP(:,2) = num2cell(a_counts);

NT_PLP_unique = num2cell(unique(cell2mat(NT_PLP(:,1))));
for i = 1:size(NT_PLP_unique,1)
    idx = find(cell2mat(NT_PLP(:,1)) == NT_PLP_unique{i,1}, 1);
    NT_PLP_unique{i,2} = NT_PLP{idx,3};
end
[NT_PLP,~,ic] = unique(NT_PLP_unique(:,2));
a_counts = accumarray(ic,1);
NT_PLP(:,2) = num2cell(a_counts);

NT_IPS_unique = num2cell(unique(cell2mat(NT_IPS(:,1))));
for i = 1:size(NT_IPS_unique,1)
    idx = find(cell2mat(NT_IPS(:,1)) == NT_IPS_unique{i,1}, 1);
    NT_IPS_unique{i,2} = NT_IPS{idx,3};
end
[NT_IPS,~,ic] = unique(NT_IPS_unique(:,2));
a_counts = accumarray(ic,1);
NT_IPS(:,2) = num2cell(a_counts);

NT_SPS_unique = num2cell(unique(cell2mat(NT_SPS(:,1))));
for i = 1:size(NT_SPS_unique,1)
    idx = find(cell2mat(NT_SPS(:,1)) == NT_SPS_unique{i,1}, 1);
    NT_SPS_unique{i,2} = NT_SPS{idx,3};
end
[NT_SPS,~,ic] = unique(NT_SPS_unique(:,2));
a_counts = accumarray(ic,1);
NT_SPS(:,2) = num2cell(a_counts);

NT_WED_unique = num2cell(unique(cell2mat(NT_WED(:,1))));
for i = 1:size(NT_WED_unique,1)
    idx = find(cell2mat(NT_WED(:,1)) == NT_WED_unique{i,1}, 1);
    NT_WED_unique{i,2} = NT_WED{idx,3};
end
[NT_WED,~,ic] = unique(NT_WED_unique(:,2));
a_counts = accumarray(ic,1);
NT_WED(:,2) = num2cell(a_counts);

NT_PVLP_unique = num2cell(unique(cell2mat(NT_PVLP(:,1))));
for i = 1:size(NT_PVLP_unique,1)
    idx = find(cell2mat(NT_PVLP(:,1)) == NT_PVLP_unique{i,1}, 1);
    NT_PVLP_unique{i,2} = NT_PVLP{idx,3};
end
[NT_PVLP,~,ic] = unique(NT_PVLP_unique(:,2));
a_counts = accumarray(ic,1);
NT_PVLP(:,2) = num2cell(a_counts);

NT_GNG_unique = num2cell(unique(cell2mat(NT_GNG(:,1))));
for i = 1:size(NT_GNG_unique,1)
    idx = find(cell2mat(NT_GNG(:,1)) == NT_GNG_unique{i,1}, 1);
    NT_GNG_unique{i,2} = NT_GNG{idx,3};
end
[NT_GNG,~,ic] = unique(NT_GNG_unique(:,2));
a_counts = accumarray(ic,1);
NT_GNG(:,2) = num2cell(a_counts);

NT_SAD_unique = num2cell(unique(cell2mat(NT_SAD(:,1))));
for i = 1:size(NT_SAD_unique,1)
    idx = find(cell2mat(NT_SAD(:,1)) == NT_SAD_unique{i,1}, 1);
    NT_SAD_unique{i,2} = NT_SAD{idx,3};
end
[NT_SAD,~,ic] = unique(NT_SAD_unique(:,2));
a_counts = accumarray(ic,1);
NT_SAD(:,2) = num2cell(a_counts);

NT_IB_unique = num2cell(unique(cell2mat(NT_IB(:,1))));
for i = 1:size(NT_IB_unique,1)
    idx = find(cell2mat(NT_IB(:,1)) == NT_IB_unique{i,1}, 1);
    NT_IB_unique{i,2} = NT_IB{idx,3};
end
[NT_IB,~,ic] = unique(NT_IB_unique(:,2));
a_counts = accumarray(ic,1);
NT_IB(:,2) = num2cell(a_counts);

NT_SPSL_unique = num2cell(unique(cell2mat(NT_SPSL(:,1))));
for i = 1:size(NT_SPSL_unique,1)
    idx = find(cell2mat(NT_SPSL(:,1)) == NT_SPSL_unique{i,1}, 1);
    NT_SPSL_unique{i,2} = NT_SPSL{idx,3};
end
[NT_SPSL,~,ic] = unique(NT_SPSL_unique(:,2));
a_counts = accumarray(ic,1);
NT_SPSL(:,2) = num2cell(a_counts);

NT_ICL_unique = num2cell(unique(cell2mat(NT_ICL(:,1))));
for i = 1:size(NT_ICL_unique,1)
    idx = find(cell2mat(NT_ICL(:,1)) == NT_ICL_unique{i,1}, 1);
    NT_ICL_unique{i,2} = NT_ICL{idx,3};
end
[NT_ICL,~,ic] = unique(NT_ICL_unique(:,2));
a_counts = accumarray(ic,1);
NT_ICL(:,2) = num2cell(a_counts);


%% figure 2: FBP neuron type distribution (counts entered manually)
figure(2);set(gcf,'Color','w')
% data
labels = {'cLP','cL','cM','OA','LT','aMe','Others'};
values = [92 46 42 16 14 11 40];


bar(values,'FaceColor',[0.2 0.4 0.6]); % basic bar chart

% X / Y axis labels
set(gca,'XTick',1:numel(labels),'XTickLabel',labels,'FontSize',12,'Box','off','TickDir','out');
xtickangle(45); % rotate X tick labels
ylabel('Count');
title('Neuron Type Distribution');
ylim([0 100])

