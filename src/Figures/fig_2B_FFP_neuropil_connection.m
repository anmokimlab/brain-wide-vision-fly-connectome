clear all; close all; clc;

baseDir = fileparts(fileparts(mfilename('fullpath')));  % resolved from this script's location (portable; no absolute paths)

% FFP / FBP neuron classification (provides RightFFP_NPIs), saved by
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
%%

for i=1:1:size(RightFFP_NPIs,1)
    idx=FAFBNeurons.root_id==RightFFP_NPIs.root_id(i);
    if ~any(idx)
        RightFFP_NPIs.nt_type{i}='Unidentified';
    elseif isempty(FAFBNeurons.nt_type{idx})
        RightFFP_NPIs.nt_type{i}='Unidentified';
    else
        RightFFP_NPIs.nt_type{i}=FAFBNeurons.nt_type{idx};
    end
end

for i=1:1:size(NT_Types,1)
    idx=strcmp(RightFFP_NPIs.nt_type,NT_Types{i,1});
    NT_Types{i,2}=sum(idx);
end
%%  From Optic lobe to higher brain region // total and NT types
FFP_Central_Connection=cell(length(FAFBNeuropils),length(FAFBNeuropils));

for i=1:1:size(RightFFP_NPIs,1)
    current_root_id=RightFFP_NPIs.root_id(i);
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

            temp=FFP_Central_Connection{idx_In,idx_Out};
            temp{end+1,1}=current_root_id;
            temp{end,2}=OutNeuropils{k,2}*syn_ratio_In;
            temp{end,3}=RightFFP_NPIs.nt_type{i};
            FFP_Central_Connection{idx_In,idx_Out}=temp;

        end
    end
end


%% Merge AME into ME
aMe_R_idx=find(strcmp(FAFBNeuropils,'AME_R'));
Me_R_idx=find(strcmp(FAFBNeuropils,'ME_R'));
aMe_L_idx=find(strcmp(FAFBNeuropils,'AME_L'));
Me_L_idx=find(strcmp(FAFBNeuropils,'ME_L'));

for i=1:1:size(FFP_Central_Connection,2)
    Me_R_temp=FFP_Central_Connection{Me_R_idx,i};
    AMe_R_Temp=FFP_Central_Connection{aMe_R_idx,i};
    if isempty(Me_R_temp)||isempty(AMe_R_Temp)
        FFP_Central_Connection{Me_R_idx,i}=[Me_R_temp;AMe_R_Temp];
    else
        [MergedAMeIdx,MergedMeIdx]=(ismember(cell2mat(AMe_R_Temp(:,1)),cell2mat(Me_R_temp(:,1))));
        MergedAMeIdx=find(MergedAMeIdx);
        MergedMeIdx(MergedMeIdx==0)=[];
        if ~isempty(MergedMeIdx)
            for j=1:1:size(MergedMeIdx,1)
                Me_R_temp{MergedMeIdx(j),2}= Me_R_temp{MergedMeIdx(j),2}+AMe_R_Temp{MergedAMeIdx(j),2};
            end
        end
        % --- Append the neurons unique to AMe_R ---
        % Use 'ismember' again to find the neurons NOT already in Me_R.
        NonMergedAMeIdx = ~ismember(cell2mat(AMe_R_Temp(:,1)), cell2mat(Me_R_temp(:,1)));

        % Append those unique neurons (NonMergedAMeIdx) to the bottom of Me_R_temp.
        Me_R_temp = [Me_R_temp; AMe_R_Temp(NonMergedAMeIdx, :)];

        % Now 'Me_R_temp' holds the merged Me_R + AMe_R neuron list.
        FFP_Central_Connection{Me_R_idx,i}=Me_R_temp;
    end
end


%% Build the optic-lobe -> central-brain connection-count matrix
%%% Me R 55, Lo R 45, LoP R 43
idx_RightOpticLobes=[55 45 43];
idx_CentralBrains=[1 2 5:37 40 41 46:53 56:74 76:79];

FFP_Central_Connection_Matrix=zeros(length(idx_RightOpticLobes),length(idx_CentralBrains));

for i=1:1:length(idx_RightOpticLobes)
    for j=1:1:length(idx_CentralBrains)
        currentConnection=FFP_Central_Connection{idx_RightOpticLobes(i),idx_CentralBrains(j)};
        if isempty(currentConnection)
            FFP_Central_Connection_Matrix(i,j)=0;
        else
            % number of distinct FFP neurons linking this optic lobe -> central region
            FFP_Central_Connection_Matrix(i,j)=size(currentConnection,1);
        end
    end
end

NeuronN_Thr=1;   % minimum neuron count required to draw a bubble (used by Figure 2)


%% sending neuron names...

FFP_Optic_Central_Connection=FFP_Central_Connection(idx_RightOpticLobes,idx_CentralBrains);

unique_neuron_size_matrix_optic_lobe=zeros(size(FFP_Optic_Central_Connection,1),1);
for i=1:1:3
    FFP_Names=[];
    for j=1:1:size(FFP_Optic_Central_Connection,2)
        FFP_Names=[FFP_Names; FFP_Optic_Central_Connection{i, j}];
    end
    if ~isempty(FFP_Names)
        unique_neuron_size_matrix_optic_lobe(i,1)=length(unique(cell2mat(FFP_Names(:,1))));
    end
end

unique_neuron_size_matrix_central_brain=zeros(size(FFP_Optic_Central_Connection,2),1);

for j=1:1:size(FFP_Optic_Central_Connection,2)
    FFP_Names=[];
    for i=1:1:3
        FFP_Names=[FFP_Names; FFP_Optic_Central_Connection{i, j}];
    end
    if ~isempty(FFP_Names)
        unique_neuron_size_matrix_central_brain(j,1)=length(unique(cell2mat(FFP_Names(:,1))));
    end
end
[~,idx_SortedCentralBrains]=sort(unique_neuron_size_matrix_central_brain,'descend');
SortedCentralBrains=FAFBNeuropil_Central(idx_SortedCentralBrains);
WantToSeeN=1:1:10;
FFP_Optic_Central_Connection_sorted=FFP_Optic_Central_Connection(:,idx_SortedCentralBrains(WantToSeeN));

FFP_Central_Connection_NeuronsNames=cell(length(idx_RightOpticLobes),length(WantToSeeN));
FFP_Central_Connection_NeuronsNamesType=cell(length(idx_RightOpticLobes),length(WantToSeeN));

for i=1:1:length(idx_RightOpticLobes)
    for j=1:1:length(WantToSeeN)
        currentConnection=FFP_Optic_Central_Connection_sorted{i,j};
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
            FFP_Central_Connection_NeuronsNames{i,j}=currentConnection;
            FFP_Central_Connection_NeuronsNamesType{i,j}=uniquecurrentConnectionType;
        end
    end
end
%% Figure 1: bubble chart of FFP optic-lobe -> central-brain connections
% Each bubble's size = number of FFP neurons linking an optic-lobe region (rows:
% Me R / Lo R / LoP R) to a central-brain region (columns: top 10 by neuron count).
figure(1);set(gcf,'Color','w')

temp_matrix=FFP_Central_Connection_Matrix(:,idx_SortedCentralBrains(WantToSeeN));

% Intentionally add an extra row (4): the total neuron count of each central
% region, so that bubble's size encodes the number of neurons in that neuropil.
temp_matrix(4,WantToSeeN)=unique_neuron_size_matrix_central_brain(idx_SortedCentralBrains(WantToSeeN))';

unique_neuron_size_matrix_optic_lobe_wanted=zeros(size(FFP_Optic_Central_Connection_sorted,1),1);
for i=1:1:3
    FFP_Names=[];
    for j=1:1:size(FFP_Optic_Central_Connection_sorted,2)
        FFP_Names=[FFP_Names; FFP_Optic_Central_Connection_sorted{i, j}];
    end
    if ~isempty(FFP_Names)
        unique_neuron_size_matrix_optic_lobe_wanted(i,1)=length(unique(cell2mat(FFP_Names(:,1))));
    end
end


% Intentionally add an extra column: the total neuron count of each optic-lobe
% region, so that bubble's size also encodes the number of neurons in that neuropil.
temp_matrix(1:1:3,length(WantToSeeN)+1)=unique_neuron_size_matrix_optic_lobe_wanted;

[YData, XData]=find(temp_matrix>=NeuronN_Thr);
FFP_Central_Connection_Matrix_bubble=temp_matrix(temp_matrix>=NeuronN_Thr);
XData=XData';
YData=YData';
FFP_Central_Connection_Matrix_bubble=FFP_Central_Connection_Matrix_bubble';
bubblechart(XData,YData,FFP_Central_Connection_Matrix_bubble,MarkerFaceAlpha=0.8,MarkerFaceColor='#48494B',MarkerEdgeColor='#48494B',MarkerEdgeAlpha=0.8);
bubblelegend('Neuron','Location','eastoutside')

set(gca,'XTick',1:1:length(WantToSeeN),'XTickLabels',FAFBNeuropil_Central(idx_SortedCentralBrains(WantToSeeN)),...
    'YTick',1:1:size(FFP_Central_Connection_Matrix,1),'YTickLabels',{'Me R','Lo R','LoP R'},'Box','off','TickDir','out')
title('Right FeedForward Neuron')
grid on

%% Per-neuropil neurotransmitter (NT) composition (source data for the NT pie charts)
% Each NT_<neuropil> variable (NT_Me, NT_Lo, NT_LoP, NT_PLP, NT_PVLP, NT_AVLP,
% NT_SPS, NT_AOTU, NT_IPS, NT_SLP, NT_IB, NT_SCL, NT_ICL) is built below as a
% {NT_type, neuron_count} list for that neuropil. These counts were copied into
% Excel to draw the NT pie charts, saved at:
%   Processed_Data\FFP_neuropil_NT.xlsx

NT_Me={};
NT_Lo={};
NT_LoP={};
NT_PLP={};
NT_PVLP={};
NT_AVLP={};
NT_SPS={};
NT_AOTU={};
NT_IPS={};
NT_SLP={};
NT_IB={};
NT_SCL={};
NT_ICL={};

for i=1:1:size(FFP_Central_Connection_NeuronsNames,2)
    NT_Me=[NT_Me;FFP_Central_Connection_NeuronsNames{1,i}];
    NT_Lo=[NT_Lo;FFP_Central_Connection_NeuronsNames{2,i}];
    NT_LoP=[NT_LoP;FFP_Central_Connection_NeuronsNames{3,i}];
end

for i=1:1:size(FFP_Central_Connection_NeuronsNames,1)
    NT_PLP=[NT_PLP;FFP_Central_Connection_NeuronsNames{i,1}];
    NT_PVLP=[NT_PVLP;FFP_Central_Connection_NeuronsNames{i,2}];
    NT_SPS=[NT_SPS;FFP_Central_Connection_NeuronsNames{i,3}];
    NT_AOTU=[NT_AOTU;FFP_Central_Connection_NeuronsNames{i,4}];
    NT_AVLP=[NT_AVLP;FFP_Central_Connection_NeuronsNames{i,5}];
    NT_SCL=[NT_SCL;FFP_Central_Connection_NeuronsNames{i,6}];
    NT_ICL=[NT_ICL;FFP_Central_Connection_NeuronsNames{i,7}];
    NT_IPS=[NT_IPS;FFP_Central_Connection_NeuronsNames{i,8}];
    NT_SLP=[NT_SLP;FFP_Central_Connection_NeuronsNames{i,9}];
    NT_IB=[NT_IB;FFP_Central_Connection_NeuronsNames{i,10}];
end

%% Me
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

% PLP
NT_PLP_unique = num2cell(unique(cell2mat(NT_PLP(:,1))));
for i = 1:size(NT_PLP_unique,1)
    idx = find(cell2mat(NT_PLP(:,1)) == NT_PLP_unique{i,1}, 1);
    NT_PLP_unique{i,2} = NT_PLP{idx,3};
end
[NT_PLP,~,ic] = unique(NT_PLP_unique(:,2));
a_counts = accumarray(ic,1);
NT_PLP(:,2) = num2cell(a_counts);

% PVLP
NT_PVLP_unique = num2cell(unique(cell2mat(NT_PVLP(:,1))));
for i = 1:size(NT_PVLP_unique,1)
    idx = find(cell2mat(NT_PVLP(:,1)) == NT_PVLP_unique{i,1}, 1);
    NT_PVLP_unique{i,2} = NT_PVLP{idx,3};
end
[NT_PVLP,~,ic] = unique(NT_PVLP_unique(:,2));
a_counts = accumarray(ic,1);
NT_PVLP(:,2) = num2cell(a_counts);

% AVLP
NT_AVLP_unique = num2cell(unique(cell2mat(NT_AVLP(:,1))));
for i = 1:size(NT_AVLP_unique,1)
    idx = find(cell2mat(NT_AVLP(:,1)) == NT_AVLP_unique{i,1}, 1);
    NT_AVLP_unique{i,2} = NT_AVLP{idx,3};
end
[NT_AVLP,~,ic] = unique(NT_AVLP_unique(:,2));
a_counts = accumarray(ic,1);
NT_AVLP(:,2) = num2cell(a_counts);

% SPS
NT_SPS_unique = num2cell(unique(cell2mat(NT_SPS(:,1))));
for i = 1:size(NT_SPS_unique,1)
    idx = find(cell2mat(NT_SPS(:,1)) == NT_SPS_unique{i,1}, 1);
    NT_SPS_unique{i,2} = NT_SPS{idx,3};
end
[NT_SPS,~,ic] = unique(NT_SPS_unique(:,2));
a_counts = accumarray(ic,1);
NT_SPS(:,2) = num2cell(a_counts);

% AOTU
NT_AOTU_unique = num2cell(unique(cell2mat(NT_AOTU(:,1))));
for i = 1:size(NT_AOTU_unique,1)
    idx = find(cell2mat(NT_AOTU(:,1)) == NT_AOTU_unique{i,1}, 1);
    NT_AOTU_unique{i,2} = NT_AOTU{idx,3};
end
[NT_AOTU,~,ic] = unique(NT_AOTU_unique(:,2));
a_counts = accumarray(ic,1);
NT_AOTU(:,2) = num2cell(a_counts);

% IPS
NT_IPS_unique = num2cell(unique(cell2mat(NT_IPS(:,1))));
for i = 1:size(NT_IPS_unique,1)
    idx = find(cell2mat(NT_IPS(:,1)) == NT_IPS_unique{i,1}, 1);
    NT_IPS_unique{i,2} = NT_IPS{idx,3};
end
[NT_IPS,~,ic] = unique(NT_IPS_unique(:,2));
a_counts = accumarray(ic,1);
NT_IPS(:,2) = num2cell(a_counts);

% SLP
NT_SLP_unique = num2cell(unique(cell2mat(NT_SLP(:,1))));
for i = 1:size(NT_SLP_unique,1)
    idx = find(cell2mat(NT_SLP(:,1)) == NT_SLP_unique{i,1}, 1);
    NT_SLP_unique{i,2} = NT_SLP{idx,3};
end
[NT_SLP,~,ic] = unique(NT_SLP_unique(:,2));
a_counts = accumarray(ic,1);
NT_SLP(:,2) = num2cell(a_counts);

% LH
NT_IB_unique = num2cell(unique(cell2mat(NT_IB(:,1))));
for i = 1:size(NT_IB_unique,1)
    idx = find(cell2mat(NT_IB(:,1)) == NT_IB_unique{i,1}, 1);
    NT_IB_unique{i,2} = NT_IB{idx,3};
end
[NT_IB,~,ic] = unique(NT_IB_unique(:,2));
a_counts = accumarray(ic,1);
NT_IB(:,2) = num2cell(a_counts);

% SCL
NT_SCL_unique = num2cell(unique(cell2mat(NT_SCL(:,1))));
for i = 1:size(NT_SCL_unique,1)
    idx = find(cell2mat(NT_SCL(:,1)) == NT_SCL_unique{i,1}, 1);
    NT_SCL_unique{i,2} = NT_SCL{idx,3};
end
[NT_SCL,~,ic] = unique(NT_SCL_unique(:,2));
a_counts = accumarray(ic,1);
NT_SCL(:,2) = num2cell(a_counts);

% ICL
NT_ICL_unique = num2cell(unique(cell2mat(NT_ICL(:,1))));
for i = 1:size(NT_ICL_unique,1)
    idx = find(cell2mat(NT_ICL(:,1)) == NT_ICL_unique{i,1}, 1);
    NT_ICL_unique{i,2} = NT_ICL{idx,3};
end
[NT_ICL,~,ic] = unique(NT_ICL_unique(:,2));
a_counts = accumarray(ic,1);
NT_ICL(:,2) = num2cell(a_counts);


%% Figure 2: FFP neuron type distribution (counts entered manually)
figure(2);set(gcf,'Color','w')
% Data definition
labels = {'LC','MTe','LLPC','MeTu','LPLC','LT','LPC','LPT','aMe','Others'};
values = [1756,352,333,317,228,210,82,47,43,68];


bar(values,'FaceColor',[0.2 0.4 0.6]); % basic bar chart

% X / Y axis labels
set(gca,'XTick',1:numel(labels),'XTickLabel',labels,'FontSize',12,'Box','off','TickDir','out');
xtickangle(45); % tilt the x-axis labels
ylabel('Count');
title('Neuron Type Distribution');

% Show the value above each bar
for i = 1:length(values)
    text(i,values(i)+30,num2str(values(i)),'HorizontalAlignment','center','FontSize',10);
end
ylim([0 2000])
