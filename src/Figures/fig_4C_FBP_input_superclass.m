%% Load data
clear all; close all; clc

baseDir = fileparts(fileparts(mfilename('fullpath')));  % resolved from this script's location (portable; no absolute paths)

% FFP / FBP / BDP neuron classification (provides RightFBP_NPIs), saved by
% Figures/fig_1D_E_postPI_prePI_right_neurons.m
load(fullfile(baseDir, 'Processed_Data', 'right_neurons_thr0.mat'))

% Group the FBP neurons by cell type
[RightFBP_type,~,ic]=unique(RightFBP_NPIs.type);
RightFBP_type=table(RightFBP_type,'VariableNames',{'type'});
for i=1:1:size(RightFBP_type,1)
    idx=ic==i;
    RightFBP_type.root_id{i}=RightFBP_NPIs.root_id(idx);
end

%% Load Codex data
opt = detectImportOptions(fullfile(baseDir, 'Codex_Data', 'connections_no_threshold.csv'));
opt = setvartype(opt,'pre_root_id','int64');
opt = setvartype(opt,'post_root_id','int64');
FAFBConnections = readtable(fullfile(baseDir, 'Codex_Data', 'connections_no_threshold.csv'),opt);

opt = detectImportOptions(fullfile(baseDir, 'Codex_Data', 'classification.csv'));
opt = setvartype(opt,'root_id','int64');
FAFBClassification = readtable(fullfile(baseDir, 'Codex_Data', 'classification.csv'),opt);

opt = detectImportOptions(fullfile(baseDir, 'Codex_Data', 'consolidated_cell_types.csv'));
opt = setvartype(opt,'root_id','int64');
FAFBConsolidated_type = readtable(fullfile(baseDir, 'Codex_Data', 'consolidated_cell_types.csv'),opt);

FAFB_Superclass=unique(FAFBClassification.super_class);

%% Collect the central-brain input neurons of each FBP type, grouped by superclass
for i=1:1:size(RightFBP_type,1)
    current_FBP_root_id=RightFBP_type.root_id{i};
    InConnections=FAFBConnections(ismember(FAFBConnections.post_root_id,current_FBP_root_id),:);
    % Keep only inputs from the central brain (drop optic-lobe neuropils, both sides)
    InConnections(ismember(InConnections.neuropil,{'LA_R';'AME_R';'ME_R';'LO_R';'LOP_R'}),:)=[];
    InConnections(ismember(InConnections.neuropil,{'LA_L';'AME_L';'ME_L';'LO_L';'LOP_L'}),:)=[];

    [Unique_inputs,~,ic]=unique(InConnections.pre_root_id);
    Unique_inputs=num2cell(Unique_inputs);
    for j=1:1:size(Unique_inputs,1)
        idx_syn=ic==j;
        Unique_inputs{j,2}=sum(InConnections.syn_count(idx_syn));
        idx_superclass=FAFBClassification.root_id==Unique_inputs{j,1};
        Unique_inputs{j,3}=FAFBClassification.super_class{idx_superclass};
    end
    In_SuperClass=cell2table(FAFB_Superclass,'VariableNames',{'superclass'});
    for j=1:1:size(In_SuperClass,1)
        idx_superclass=strcmp(Unique_inputs(:,3),In_SuperClass.superclass{j});
        In_SuperClass.syn_count{j}=sum(cell2mat(Unique_inputs(idx_superclass,2)));
        In_SuperClass.root_ids{j}=Unique_inputs(idx_superclass,[1 2]);
    end
    RightFBP_type.input_superclass{i}=In_SuperClass;
end

%% Colors (one per superclass, in FAFB_Superclass / alphabetical order)
colors=[ 0.6350, 0.5090, 0.2540;... % ascending
    0.4940, 0.1840, 0.5560;...      % central
    0.8500, 0.1500, 0.2000;...      % descending
    1.00,0.07,0.65;...              % endocrine
    0.8500, 0.3250, 0.0980;         % motor
    0.9290, 0.6940, 0.1250;...      % optic
    0.3010, 0.7450, 0.9330;...      % sensory
    0.4660, 0.6740, 0.1880;...      % visual_centrifugal
    0.0000, 0.4470, 0.7410          % visual_projection
    ];

%% Group the FBP types by their target optic lobe: Me / Lo / Lop / Multi
% (index lists into RightFBP_type, from the FBP output-neuropil classification, s10)
Me_FBP_idx=[4 22 25 26 29 30 31 32 33 68 69 70 71 72 73 74 75 76 78 81 82 83 84 85 86 87 93];
Lo_FBP_idx=[2 10 11 12 13 14 15 16 17 18 19 21 34 35 36 37 38 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 95];
Lop_FBP_idx=[3 7 8 9 28 62 63 64 65 66];
Multi_FBP_idx=[1 5 6 23 24 27 39 60 61 67 77 79 80 88 89 90 91 92 94 96];

% Group-averaged input superclass composition (normalized per neuron, then averaged)
Me_FBP_input_superclass=zeros(size(FAFB_Superclass));
for i=1:1:size(Me_FBP_idx,2)
    temp=cell2mat(RightFBP_type.input_superclass{Me_FBP_idx(i)}.syn_count);
    temp=temp/sum(temp);
    Me_FBP_input_superclass=Me_FBP_input_superclass +temp;
end
Me_FBP_input_superclass=Me_FBP_input_superclass./i;
Me_FBP_input_superclass=[FAFB_Superclass num2cell(Me_FBP_input_superclass)];

Lo_FBP_input_superclass=zeros(size(FAFB_Superclass));
for i=1:1:size(Lo_FBP_idx,2)
    temp=cell2mat(RightFBP_type.input_superclass{Lo_FBP_idx(i)}.syn_count);
    temp=temp/sum(temp);
    Lo_FBP_input_superclass=Lo_FBP_input_superclass +temp;
end
Lo_FBP_input_superclass=Lo_FBP_input_superclass./i;
Lo_FBP_input_superclass=[FAFB_Superclass num2cell(Lo_FBP_input_superclass)];

Lop_FBP_input_superclass=zeros(size(FAFB_Superclass));
for i=1:1:size(Lop_FBP_idx,2)
    temp=cell2mat(RightFBP_type.input_superclass{Lop_FBP_idx(i)}.syn_count);
    temp=temp/sum(temp);
    Lop_FBP_input_superclass=Lop_FBP_input_superclass +temp;
end
Lop_FBP_input_superclass=Lop_FBP_input_superclass./i;
Lop_FBP_input_superclass=[FAFB_Superclass num2cell(Lop_FBP_input_superclass)];

Multi_FBP_input_superclass=zeros(size(FAFB_Superclass));
for i=1:1:size(Multi_FBP_idx,2)
    temp=cell2mat(RightFBP_type.input_superclass{Multi_FBP_idx(i)}.syn_count);
    temp=temp/sum(temp);
    Multi_FBP_input_superclass=Multi_FBP_input_superclass +temp;
end
Multi_FBP_input_superclass=Multi_FBP_input_superclass./i;
Multi_FBP_input_superclass=[FAFB_Superclass num2cell(Multi_FBP_input_superclass)];

All_idx=[Me_FBP_idx Lo_FBP_idx Lop_FBP_idx Multi_FBP_idx];

%% Figure 1 (paper S2A): per-type input superclass composition
figure(1);set(gcf,'Color','w'); hold on;
for i=1:1:size(All_idx,2)
    values=cell2mat(RightFBP_type.input_superclass{All_idx(i)}.syn_count);
    values=values/sum(values);
    b= bar(i,values,'stacked','EdgeColor','flat');

    for j = 1:length(b)
        b(j).FaceColor = 'flat';   % set to 'flat'
        b(j).CData = colors(j, :); % color of each stack
    end
end
hold off;
set(gca,'TickDir','out','XTick',1:1:size(RightFBP_type,1),'XTickLabel',RightFBP_type.type(All_idx))
title('Input')
ylim([0 1])

%% Figure 2 (paper 4C): per-group input superclass composition
figure(2);set(gcf,'Color','w'); hold on;

b= bar(1,cell2mat(Me_FBP_input_superclass(:,2)),'stacked','EdgeColor','flat');
for j = 1:length(b)
    b(j).FaceColor = 'flat';
    b(j).CData = colors(j, :);
end

b= bar(2,cell2mat(Lo_FBP_input_superclass(:,2)),'stacked','EdgeColor','flat');
for j = 1:length(b)
    b(j).FaceColor = 'flat';
    b(j).CData = colors(j, :);
end

b= bar(3,cell2mat(Lop_FBP_input_superclass(:,2)),'stacked','EdgeColor','flat');
for j = 1:length(b)
    b(j).FaceColor = 'flat';
    b(j).CData = colors(j, :);
end

b= bar(4,cell2mat(Multi_FBP_input_superclass(:,2)),'stacked','EdgeColor','flat');
for j = 1:length(b)
    b(j).FaceColor = 'flat';
    b(j).CData = colors(j, :);
end
hold off;

set(gca,'TickDir','out','XTick',1:1:4,'XTickLabel',{'Me FBP','Lo FBP','Lop FBP','Multi FBP'})
ylim([0 1])
title('Input')

%% Reassign each type's 'central' input by the superclass of those central neurons' own top-5 inputs
RightFBP_type.input_superclass2=RightFBP_type.input_superclass;
TopN=5;
for i=1:1:size(RightFBP_type,1)
    centralsToFBP=RightFBP_type.input_superclass{i}.root_ids{2};   % column 2 = 'central' superclass
    for j=size(centralsToFBP,1):-1:1
        inConnections_centralsToFBP=FAFBConnections(FAFBConnections.post_root_id==centralsToFBP{j,1},:);
        for k=1:1:size(inConnections_centralsToFBP,1)
            idx_type=FAFBConsolidated_type.root_id==inConnections_centralsToFBP.pre_root_id(k);
            if any(idx_type)
                inConnections_centralsToFBP.pre_type{k}=FAFBConsolidated_type.primary_type{idx_type};
            else
                inConnections_centralsToFBP.pre_type{k}='NoLabel';
            end
            idx_superclass=FAFBClassification.root_id==inConnections_centralsToFBP.pre_root_id(k);
            inConnections_centralsToFBP.pre_superclass{k}=FAFBClassification.super_class{idx_superclass};
        end
        if isempty(inConnections_centralsToFBP)
            continue;
        end
        noLabel_idx= strcmp(inConnections_centralsToFBP.pre_type,'NoLabel');
        inConnections_centralsToFBP.pre_superclass(noLabel_idx)={'NoLabel'};

        [unique_inputType_centralsToFBP,~,ic]=unique(inConnections_centralsToFBP.pre_type);
        for k=1:1:size(unique_inputType_centralsToFBP,1)
            idx=ic==k;
            unique_inputType_centralsToFBP{k,2}=sum(inConnections_centralsToFBP.syn_count(idx));
            [unique_inputType_centralsToFBP_superclass,~,icc]=unique(inConnections_centralsToFBP.pre_superclass(idx));
            if length(unique_inputType_centralsToFBP_superclass)~=1
                temp=inConnections_centralsToFBP(idx,:);
                for l=1:1:size(unique_inputType_centralsToFBP_superclass,1)
                    idxx=icc==l;
                    unique_inputType_centralsToFBP_superclass{l,2}=sum(temp.syn_count(idxx));
                end
                [~,idx_max]=max(cell2mat(unique_inputType_centralsToFBP_superclass(:,2)));
                unique_inputType_centralsToFBP{k,3}=unique_inputType_centralsToFBP_superclass{idx_max,1};
            else
                unique_inputType_centralsToFBP(k,3)=unique_inputType_centralsToFBP_superclass;
            end
        end
        unique_inputType_centralsToFBP=sortrows(unique_inputType_centralsToFBP,2,'descend');
        TopN_inputType_centralsToFBP=unique_inputType_centralsToFBP(1:1:min(TopN,size(unique_inputType_centralsToFBP,1)),:);
        superClass_Input_centralsToFBP=cell2table(FAFB_Superclass,'VariableNames',{'superclass'});
        for k=1:1:size(superClass_Input_centralsToFBP,1)
            idx_superclass=strcmp(TopN_inputType_centralsToFBP(:,3),superClass_Input_centralsToFBP.superclass{k});
            superClass_Input_centralsToFBP.syn_count{k}=sum(cell2mat(TopN_inputType_centralsToFBP(idx_superclass,2)));
        end
        % Redistribute by ratio, excluding 'central' itself
        superClass_Input_centralsToFBP.syn_count{2}=0;
        if sum(cell2mat(superClass_Input_centralsToFBP.syn_count))==0
            continue;
        else
            RatioSuperClass=cell2mat(superClass_Input_centralsToFBP.syn_count)/sum(cell2mat(superClass_Input_centralsToFBP.syn_count));
            for k=1:1:size(RatioSuperClass,1)
                if RatioSuperClass(k)==0
                    continue;
                else
                    RightFBP_type.input_superclass2{i}.syn_count{k}=...
                        RightFBP_type.input_superclass2{i}.syn_count{k}+...
                        RightFBP_type.input_superclass2{i}.root_ids{2,1}{j,2}*RatioSuperClass(k);
                    RightFBP_type.input_superclass2{i}.syn_count{2}=...
                        RightFBP_type.input_superclass2{i}.syn_count{2}-...
                        RightFBP_type.input_superclass2{i}.root_ids{2,1}{j,2}*RatioSuperClass(k);
                end
            end

            RightFBP_type.input_superclass2{i}.root_ids{2,1}(j, :)=[];

        end
    end

end

%% Figure 3 (paper S2A): per-type input superclass after central reassignment
figure(3);set(gcf,'Color','w'); hold on;
for i=1:1:size(All_idx,2)
    values=cell2mat(RightFBP_type.input_superclass2{All_idx(i)}.syn_count);
    values=values/sum(values);
    b= bar(i,values,'stacked','EdgeColor','flat');

    for j = 1:length(b)
        b(j).FaceColor = 'flat';
        b(j).CData = colors(j, :);
    end
end
hold off;
ylim([0 1])
set(gca,'TickDir','out','XTick',1:1:size(RightFBP_type,1),'XTickLabel',RightFBP_type.type(All_idx))
title('center neuron replaced input')

%% Group-averaged input superclass composition after central reassignment
Me_FBP_input_superclass2=zeros(size(FAFB_Superclass));
for i=1:1:size(Me_FBP_idx,2)
    temp=cell2mat(RightFBP_type.input_superclass2{Me_FBP_idx(i)}.syn_count);
    temp=temp/sum(temp);
    Me_FBP_input_superclass2=Me_FBP_input_superclass2 +temp;
end
Me_FBP_input_superclass2=Me_FBP_input_superclass2./i;
Me_FBP_input_superclass2=[FAFB_Superclass num2cell(Me_FBP_input_superclass2)];

Lo_FBP_input_superclass2=zeros(size(FAFB_Superclass));
for i=1:1:size(Lo_FBP_idx,2)
    temp=cell2mat(RightFBP_type.input_superclass2{Lo_FBP_idx(i)}.syn_count);
    temp=temp/sum(temp);
    Lo_FBP_input_superclass2=Lo_FBP_input_superclass2 +temp;
end
Lo_FBP_input_superclass2=Lo_FBP_input_superclass2./i;
Lo_FBP_input_superclass2=[FAFB_Superclass num2cell(Lo_FBP_input_superclass2)];

Lop_FBP_input_superclass2=zeros(size(FAFB_Superclass));
for i=1:1:size(Lop_FBP_idx,2)
    temp=cell2mat(RightFBP_type.input_superclass2{Lop_FBP_idx(i)}.syn_count);
    temp=temp/sum(temp);
    Lop_FBP_input_superclass2=Lop_FBP_input_superclass2 +temp;
end
Lop_FBP_input_superclass2=Lop_FBP_input_superclass2./i;
Lop_FBP_input_superclass2=[FAFB_Superclass num2cell(Lop_FBP_input_superclass2)];

Multi_FBP_input_superclass2=zeros(size(FAFB_Superclass));
for i=1:1:size(Multi_FBP_idx,2)
    temp=cell2mat(RightFBP_type.input_superclass2{Multi_FBP_idx(i)}.syn_count);
    temp=temp/sum(temp);
    Multi_FBP_input_superclass2=Multi_FBP_input_superclass2 +temp;
end
Multi_FBP_input_superclass2=Multi_FBP_input_superclass2./i;
Multi_FBP_input_superclass2=[FAFB_Superclass num2cell(Multi_FBP_input_superclass2)];

FBP_input_superclass2=zeros(size(FAFB_Superclass));
for i=1:1:size(RightFBP_type,1)
    temp=cell2mat(RightFBP_type.input_superclass2{i}.syn_count);
    temp=temp/sum(temp);
    FBP_input_superclass2=FBP_input_superclass2 +temp;
end
FBP_input_superclass2=FBP_input_superclass2./i;
FBP_input_superclass2=[FAFB_Superclass num2cell(FBP_input_superclass2)];

%% Figure 4 (paper 4C): per-group input superclass after central reassignment
figure(4);set(gcf,'Color','w'); hold on;

b= bar(1,cell2mat(Me_FBP_input_superclass2(:,2)),'stacked','EdgeColor','flat');
for j = 1:length(b)
    b(j).FaceColor = 'flat';
    b(j).CData = colors(j, :);
end

b= bar(2,cell2mat(Lo_FBP_input_superclass2(:,2)),'stacked','EdgeColor','flat');
for j = 1:length(b)
    b(j).FaceColor = 'flat';
    b(j).CData = colors(j, :);
end

b= bar(3,cell2mat(Lop_FBP_input_superclass2(:,2)),'stacked','EdgeColor','flat');
for j = 1:length(b)
    b(j).FaceColor = 'flat';
    b(j).CData = colors(j, :);
end

b= bar(4,cell2mat(Multi_FBP_input_superclass2(:,2)),'stacked','EdgeColor','flat');
for j = 1:length(b)
    b(j).FaceColor = 'flat';
    b(j).CData = colors(j, :);
end
hold off;

set(gca,'TickDir','out','XTick',1:1:4,'XTickLabel',{'Me FBP','Lo FBP','Lop FBP','Multi FBP'})
ylim([0 1])
title('center neuron replaced input')

%% Selected-superclass stacked bars: per-group input (original)
data1 = containers.Map(FAFB_Superclass, cell2mat(Me_FBP_input_superclass(:,2)));
data2 = containers.Map(FAFB_Superclass, cell2mat(Lo_FBP_input_superclass(:,2)));
data3 = containers.Map(FAFB_Superclass, cell2mat(Lop_FBP_input_superclass(:,2)));
data4 = containers.Map(FAFB_Superclass, cell2mat(Multi_FBP_input_superclass(:,2)));
dataList = {data1, data2, data3, data4};
selected = {'ascending' 'central' 'descending' 'visual_centrifugal' 'visual_projection'};

colorCell = mat2cell(colors, ones(1, size(colors, 1)), 3);
colorMap = containers.Map(FAFB_Superclass, colorCell);

% Extract the colors for the selected labels only
SelectedColors = zeros(numel(selected), 3);
for i = 1:numel(selected)
    SelectedColors(i, :) = colorMap(selected{i});
end

plotStackedSelectedBar(dataList, selected, SelectedColors);

%% Selected-superclass stacked bars: per-group input (central reassigned)
data1 = containers.Map(FAFB_Superclass, cell2mat(Me_FBP_input_superclass2(:,2)));
data2 = containers.Map(FAFB_Superclass, cell2mat(Lo_FBP_input_superclass2(:,2)));
data3 = containers.Map(FAFB_Superclass, cell2mat(Lop_FBP_input_superclass2(:,2)));
data4 = containers.Map(FAFB_Superclass, cell2mat(Multi_FBP_input_superclass2(:,2)));
dataList = {data1, data2, data3, data4};
selected = {'ascending' 'central' 'descending' 'visual_centrifugal' 'visual_projection'};

colorCell = mat2cell(colors, ones(1, size(colors, 1)), 3);
colorMap = containers.Map(FAFB_Superclass, colorCell);

% Extract the colors for the selected labels only
SelectedColors = zeros(numel(selected), 3);
for i = 1:numel(selected)
    SelectedColors(i, :) = colorMap(selected{i});
end

plotStackedSelectedBar(dataList, selected, SelectedColors);

%% Local function
function plotStackedSelectedBar(dataList, selectedLabels, colors)
    % dataList:       cell array of containers.Map objects (or a struct array)
    % selectedLabels: cell array of the labels of interest
    % colors:         RGB color per selected label (n x 3 matrix)

    numSamples = numel(dataList);
    numSelected = numel(selectedLabels);
    stackedData = zeros(numSamples, numSelected + 1);  % +1 for "Others"

    allKeys = dataList{1}.keys;  % assume all maps share the same keys

    % Fill the data
    for i = 1:numSamples
        if isstruct(dataList)
            dataMap = dataList(i).data;
        else
            dataMap = dataList{i};
        end
        otherSum = 0;
        for k = 1:numel(allKeys)
            label = allKeys{k};
            value = dataMap(label);
            selIdx = find(strcmp(selectedLabels, label));
            if ~isempty(selIdx)
                stackedData(i, selIdx) = value;
            else
                otherSum = otherSum + value;
            end
        end
        stackedData(i, end) = otherSum;  % store "Others" in the last column
    end

    % Colors (selected colors + an "Others" color)
    othersColor = [137 137 137]/255;
    allColors = [colors; othersColor];

    % Plot
    figure; set(gcf,'Color','w');
    b = bar(stackedData, 'stacked', 'FaceColor', 'flat','EdgeColor','flat');
    for j = 1:(numSelected + 1)
        b(j).FaceColor = allColors(j, :);
    end
    set(gca,'TickDir','out','Box','off');
    xlabel('Sample Index');
    ylabel('Value');
    legend([selectedLabels, {'Others'}], 'Location', 'northeastoutside');
    ylim([0 1])
end
