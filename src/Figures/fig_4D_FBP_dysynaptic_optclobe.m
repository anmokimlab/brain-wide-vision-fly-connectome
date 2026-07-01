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

opt = detectImportOptions(fullfile(baseDir, 'Codex_Data', 'consolidated_cell_types.csv'));
opt = setvartype(opt,'root_id','int64');
FAFBConsolidated_type = readtable(fullfile(baseDir, 'Codex_Data', 'consolidated_cell_types.csv'),opt);

%% Group FBP types by their target optic lobe (Me / Lo / Lop / Multi)
% (index lists into RightFBP_type, from the FBP output-neuropil classification, s10)
Me_FBP_idx=[4 22 25 26 29 30 31 32 33 68 69 70 71 72 73 74 75 76 78 81 82 83 84 85 86 87 93];
Lo_FBP_idx=[2 10 11 12 13 14 15 16 17 18 19 21 34 35 36 37 38 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 95];
Lop_FBP_idx=[3 7 8 9 28 62 63 64 65 66];
Multi_FBP_idx=[1 5 6 23 24 27 39 60 61 67 77 79 80 88 89 90 91 92 94 96];

opticlobes={'AME_R','ME_R','LO_R','LOP_R','AME_L','ME_L','LO_L','LOP_L'};

FBP_Me=RightFBP_type(Me_FBP_idx,:);
FBP_Lo=RightFBP_type(Lo_FBP_idx,:);
FBP_Lop=RightFBP_type(Lop_FBP_idx,:);
FBP_Multi=RightFBP_type(Multi_FBP_idx,:);

% Disynaptic optic-lobe input:
%   optic lobe -> upstream (central, non-optic) neuron -> FBP neuron.
% For each FBP type, find its upstream non-optic partners, then measure how much
% optic-lobe input those upstream neurons receive, weighted by their synapse
% contribution to the FBP neurons. Columns: Me R, Me L, Lo R, Lo L, LoP R, LoP L.

%% Lop-targeting FBP types
for i=1:1:size(FBP_Lop,1)
    Wantrootids=FBP_Lop.root_id{i};

    [upstreamNeurons]=seeConnection_root_id_NoOptic(Wantrootids,FAFBConnections,FAFBConsolidated_type);
    opticlobesSynapses=zeros(size(upstreamNeurons,1),size(opticlobes,2)-2);

    for j=1:1:size(upstreamNeurons,1)
        currentUpstream_root_id=upstreamNeurons(j,1);
        idx=ismember(FAFBConnections.post_root_id,currentUpstream_root_id);

        InConnections_upstream=FAFBConnections(idx,:);
        Total_InConnections_upstream=sum(InConnections_upstream.syn_count);

        opticlobesSynapses(j,1)=(sum(InConnections_upstream.syn_count(strcmp(InConnections_upstream.neuropil,'AME_R')))+sum(InConnections_upstream.syn_count(strcmp(InConnections_upstream.neuropil,'ME_R'))))/Total_InConnections_upstream*upstreamNeurons(j,2);
        opticlobesSynapses(j,2)=(sum(InConnections_upstream.syn_count(strcmp(InConnections_upstream.neuropil,'AME_L')))+sum(InConnections_upstream.syn_count(strcmp(InConnections_upstream.neuropil,'ME_L'))))/Total_InConnections_upstream*upstreamNeurons(j,2);

        opticlobesSynapses(j,3)=sum(InConnections_upstream.syn_count(strcmp(InConnections_upstream.neuropil,'LO_R')))/Total_InConnections_upstream*upstreamNeurons(j,2);
        opticlobesSynapses(j,4)=sum(InConnections_upstream.syn_count(strcmp(InConnections_upstream.neuropil,'LO_L')))/Total_InConnections_upstream*upstreamNeurons(j,2);

        opticlobesSynapses(j,5)=sum(InConnections_upstream.syn_count(strcmp(InConnections_upstream.neuropil,'LOP_R')))/Total_InConnections_upstream*upstreamNeurons(j,2);
        opticlobesSynapses(j,6)=sum(InConnections_upstream.syn_count(strcmp(InConnections_upstream.neuropil,'LOP_L')))/Total_InConnections_upstream*upstreamNeurons(j,2);

    end

    FBP_Lop.inputsFromOpticlobes{i}=sum(opticlobesSynapses)/sum(upstreamNeurons(:,2))*100;

end

% Per-type normalized LoP optic-lobe input fractions (%)
inputsFromOpticlobes_LOP=[];
for i=1:1:size(FBP_Lop,1)
    inputsFromOpticlobes_LOP=[inputsFromOpticlobes_LOP;FBP_Lop.inputsFromOpticlobes{i}/sum(FBP_Lop.inputsFromOpticlobes{i})*100];
end

%% Lo-targeting FBP types
for i=1:1:size(FBP_Lo,1)
    Wantrootids=FBP_Lo.root_id{i};

    [upstreamNeurons]=seeConnection_root_id_NoOptic(Wantrootids,FAFBConnections,FAFBConsolidated_type);
    opticlobesSynapses=zeros(size(upstreamNeurons,1),size(opticlobes,2)-2);

    for j=1:1:size(upstreamNeurons,1)
        currentUpstream_root_id=upstreamNeurons(j,1);
        idx=ismember(FAFBConnections.post_root_id,currentUpstream_root_id);

        InConnections_upstream=FAFBConnections(idx,:);
        Total_InConnections_upstream=sum(InConnections_upstream.syn_count);
        opticlobesSynapses(j,1)=(sum(InConnections_upstream.syn_count(strcmp(InConnections_upstream.neuropil,'AME_R')))+sum(InConnections_upstream.syn_count(strcmp(InConnections_upstream.neuropil,'ME_R'))))/Total_InConnections_upstream*upstreamNeurons(j,2);
        opticlobesSynapses(j,2)=(sum(InConnections_upstream.syn_count(strcmp(InConnections_upstream.neuropil,'AME_L')))+sum(InConnections_upstream.syn_count(strcmp(InConnections_upstream.neuropil,'ME_L'))))/Total_InConnections_upstream*upstreamNeurons(j,2);

        opticlobesSynapses(j,3)=sum(InConnections_upstream.syn_count(strcmp(InConnections_upstream.neuropil,'LO_R')))/Total_InConnections_upstream*upstreamNeurons(j,2);
        opticlobesSynapses(j,4)=sum(InConnections_upstream.syn_count(strcmp(InConnections_upstream.neuropil,'LO_L')))/Total_InConnections_upstream*upstreamNeurons(j,2);

        opticlobesSynapses(j,5)=sum(InConnections_upstream.syn_count(strcmp(InConnections_upstream.neuropil,'LOP_R')))/Total_InConnections_upstream*upstreamNeurons(j,2);
        opticlobesSynapses(j,6)=sum(InConnections_upstream.syn_count(strcmp(InConnections_upstream.neuropil,'LOP_L')))/Total_InConnections_upstream*upstreamNeurons(j,2);

    end
    FBP_Lo.inputsFromOpticlobes{i}=sum(opticlobesSynapses)/sum(upstreamNeurons(:,2))*100;
end

% Per-type normalized Lo optic-lobe input fractions (%)
inputsFromOpticlobes_LO=[];
for i=1:1:size(FBP_Lo,1)
    inputsFromOpticlobes_LO=[inputsFromOpticlobes_LO;FBP_Lo.inputsFromOpticlobes{i}/sum(FBP_Lo.inputsFromOpticlobes{i})*100];
end

%% Me-targeting FBP types
for i=1:1:size(FBP_Me,1)
    Wantrootids=FBP_Me.root_id{i};

    [upstreamNeurons]=seeConnection_root_id_NoOptic(Wantrootids,FAFBConnections,FAFBConsolidated_type);
    opticlobesSynapses=zeros(size(upstreamNeurons,1),size(opticlobes,2)-2);

    for j=1:1:size(upstreamNeurons,1)
        currentUpstream_root_id=upstreamNeurons(j,1);
        idx=ismember(FAFBConnections.post_root_id,currentUpstream_root_id);

        InConnections_upstream=FAFBConnections(idx,:);
        Total_InConnections_upstream=sum(InConnections_upstream.syn_count);
        opticlobesSynapses(j,1)=(sum(InConnections_upstream.syn_count(strcmp(InConnections_upstream.neuropil,'AME_R')))+sum(InConnections_upstream.syn_count(strcmp(InConnections_upstream.neuropil,'ME_R'))))/Total_InConnections_upstream*upstreamNeurons(j,2);
        opticlobesSynapses(j,2)=(sum(InConnections_upstream.syn_count(strcmp(InConnections_upstream.neuropil,'AME_L')))+sum(InConnections_upstream.syn_count(strcmp(InConnections_upstream.neuropil,'ME_L'))))/Total_InConnections_upstream*upstreamNeurons(j,2);

        opticlobesSynapses(j,3)=sum(InConnections_upstream.syn_count(strcmp(InConnections_upstream.neuropil,'LO_R')))/Total_InConnections_upstream*upstreamNeurons(j,2);
        opticlobesSynapses(j,4)=sum(InConnections_upstream.syn_count(strcmp(InConnections_upstream.neuropil,'LO_L')))/Total_InConnections_upstream*upstreamNeurons(j,2);

        opticlobesSynapses(j,5)=sum(InConnections_upstream.syn_count(strcmp(InConnections_upstream.neuropil,'LOP_R')))/Total_InConnections_upstream*upstreamNeurons(j,2);
        opticlobesSynapses(j,6)=sum(InConnections_upstream.syn_count(strcmp(InConnections_upstream.neuropil,'LOP_L')))/Total_InConnections_upstream*upstreamNeurons(j,2);

    end
    FBP_Me.inputsFromOpticlobes{i}=sum(opticlobesSynapses)/sum(upstreamNeurons(:,2))*100;
end

% Per-type normalized Me optic-lobe input fractions (%)
inputsFromOpticlobes_Me=[];
for i=1:1:size(FBP_Me,1)
    inputsFromOpticlobes_Me=[inputsFromOpticlobes_Me;FBP_Me.inputsFromOpticlobes{i}/sum(FBP_Me.inputsFromOpticlobes{i})*100];
end

%% Multi-targeting FBP types
for i=1:1:size(FBP_Multi,1)
    Wantrootids=FBP_Multi.root_id{i};

    [upstreamNeurons]=seeConnection_root_id_NoOptic(Wantrootids,FAFBConnections,FAFBConsolidated_type);
    opticlobesSynapses=zeros(size(upstreamNeurons,1),size(opticlobes,2)-2);

    for j=1:1:size(upstreamNeurons,1)
        currentUpstream_root_id=upstreamNeurons(j,1);
        idx=ismember(FAFBConnections.post_root_id,currentUpstream_root_id);

        InConnections_upstream=FAFBConnections(idx,:);
        Total_InConnections_upstream=sum(InConnections_upstream.syn_count);

        opticlobesSynapses(j,1)=(sum(InConnections_upstream.syn_count(strcmp(InConnections_upstream.neuropil,'AME_R')))+sum(InConnections_upstream.syn_count(strcmp(InConnections_upstream.neuropil,'ME_R'))))/Total_InConnections_upstream*upstreamNeurons(j,2);
        opticlobesSynapses(j,2)=(sum(InConnections_upstream.syn_count(strcmp(InConnections_upstream.neuropil,'AME_L')))+sum(InConnections_upstream.syn_count(strcmp(InConnections_upstream.neuropil,'ME_L'))))/Total_InConnections_upstream*upstreamNeurons(j,2);

        opticlobesSynapses(j,3)=sum(InConnections_upstream.syn_count(strcmp(InConnections_upstream.neuropil,'LO_R')))/Total_InConnections_upstream*upstreamNeurons(j,2);
        opticlobesSynapses(j,4)=sum(InConnections_upstream.syn_count(strcmp(InConnections_upstream.neuropil,'LO_L')))/Total_InConnections_upstream*upstreamNeurons(j,2);

        opticlobesSynapses(j,5)=sum(InConnections_upstream.syn_count(strcmp(InConnections_upstream.neuropil,'LOP_R')))/Total_InConnections_upstream*upstreamNeurons(j,2);
        opticlobesSynapses(j,6)=sum(InConnections_upstream.syn_count(strcmp(InConnections_upstream.neuropil,'LOP_L')))/Total_InConnections_upstream*upstreamNeurons(j,2);

    end
    FBP_Multi.inputsFromOpticlobes{i}=sum(opticlobesSynapses)/sum(upstreamNeurons(:,2))*100;

end

% Per-type normalized Multi optic-lobe input fractions (%)
inputsFromOpticlobes_Multi=[];
for i=1:1:size(FBP_Multi,1)
    inputsFromOpticlobes_Multi=[inputsFromOpticlobes_Multi;FBP_Multi.inputsFromOpticlobes{i}/sum(FBP_Multi.inputsFromOpticlobes{i})*100];
end

%% Figure 1 (paper 4D): mean +/- std disynaptic optic-lobe input per neuropil
matrices = {inputsFromOpticlobes_Me, inputsFromOpticlobes_LO, ...
            inputsFromOpticlobes_LOP, inputsFromOpticlobes_Multi};
conditionNames = {'CB-Me','CB-LO', 'CB-LOP', 'CB-Multi'};

% Convert NaN to 0
for i = 1:length(matrices)
    matrices{i}(isnan(matrices{i})) = 0;
end

num_neuropils = size(matrices{1}, 2);
neuropilLabels = {'Me R','Me L','Lo R','Lo L','LoP R','LoP L'};

% Mean and standard deviation per neuropil (rows) and condition (columns)
meanInputs = zeros(num_neuropils, length(matrices));
stdInputs = zeros(num_neuropils, length(matrices));
for i = 1:length(matrices)
    meanInputs(:, i) = mean(matrices{i}, 1)';     % column mean
    stdInputs(:, i) = std(matrices{i}, 0, 1)';    % column std
end

% Visualize (imagesc + per-cell text)
figure(1); set(gcf,'Color','w');
imagesc(meanInputs);
colormap(gray(256));
colorbar;
title('Mean ± STD of Inputs per Neuropil');
xlabel('Condition');
ylabel('Neuropils');
set(gca, 'XTick', 1:length(conditionNames), 'XTickLabel', conditionNames);
set(gca, 'YTick', 1:num_neuropils, 'YTickLabel', neuropilLabels);

% Add a text label to each cell
for i = 1:num_neuropils
    for j = 1:length(matrices)
        text(j, i, sprintf('%.2f±%.2f', meanInputs(i,j), stdInputs(i,j)), ...
            'HorizontalAlignment', 'center', 'FontSize', 9, 'Color', 'black');
    end
end

%% Figure 2 (paper S2B): Me-targeting FBP types, Me input (contra L vs ipsi R)
figure(2); set(gcf,'Color','w')
for i=1:1:size(FBP_Me,1)
    inputFromME_R(i)= FBP_Me.inputsFromOpticlobes{i}(1);
    inputFromME_L(i)= FBP_Me.inputsFromOpticlobes{i}(2);
end
y = [inputFromME_L' inputFromME_R'];

b = bar(FBP_Me.type, y);
b(1).FaceColor = [0.4660 0.6740 0.1880];  % left (contralateral)
b(2).FaceColor = [0 0.4470 0.7410];       % right (ipsilateral)

set(gca, 'TickDir', 'out', 'Box', 'off');
ylim([0 55]);
ylabel('Input from Me (%)');
title('Input');

%% Figure 3 (paper S2B): Lo-targeting FBP types, Lo input (contra L vs ipsi R)
figure(3); set(gcf,'Color','w')
for i=1:1:size(FBP_Lo,1)
    inputFromLO_R(i)= FBP_Lo.inputsFromOpticlobes{i}(3);
    inputFromLO_L(i)= FBP_Lo.inputsFromOpticlobes{i}(4);
end
y = [inputFromLO_L' inputFromLO_R'];

b = bar(FBP_Lo.type, y);
b(1).FaceColor = [0.4660 0.6740 0.1880];  % left (contralateral)
b(2).FaceColor = [0 0.4470 0.7410];       % right (ipsilateral)

set(gca, 'TickDir', 'out', 'Box', 'off');
ylim([0 55]);
ylabel('Input from Lo (%)');
title('Input');

%% Figure 4 (paper S2B): Lop-targeting FBP types, LoP input (contra L vs ipsi R)
figure(4); set(gcf,'Color','w')
for i=1:1:size(FBP_Lop,1)
    inputFromLOP_R(i)= FBP_Lop.inputsFromOpticlobes{i}(5);
    inputFromLOP_L(i)= FBP_Lop.inputsFromOpticlobes{i}(6);
end
y = [inputFromLOP_L' inputFromLOP_R'];

b = bar(FBP_Lop.type, y);
b(1).FaceColor = [0.4660 0.6740 0.1880];  % left (contralateral)
b(2).FaceColor = [0 0.4470 0.7410];       % right (ipsilateral)

set(gca, 'TickDir', 'out', 'Box', 'off');
ylim([0 55]);
ylabel('Input from LoP (%)');
title('Input');


%% Local function
function [upstreamNeurons] = seeConnection_root_id_NoOptic(Want_root_ids,FAFBConnections,FAFB_consolidated_cell_types)
% Return the upstream (presynaptic) partners of Want_root_ids, excluding optic-lobe
% inputs. upstreamNeurons: column 1 = pre_root_id, column 2 = total synapse count
% onto Want_root_ids, sorted by synapse count (descending).

idx=ismember(FAFBConnections.post_root_id,Want_root_ids);
upstreamConnection=FAFBConnections(idx,:);
opticlobes={'LA_R','AME_R','ME_R','LO_R','LOP_R','LA_L','AME_L','ME_L','LO_L','LOP_L'};
upstreamConnection(ismember(upstreamConnection.neuropil,opticlobes),:)=[];

[upstreamNeurons,~,ic]=unique(upstreamConnection.pre_root_id);
for i=1:1:size(upstreamNeurons,1)
    idx=ic==i;
    upstreamNeurons(i,2)=sum(upstreamConnection.syn_count(idx));
end
    upstreamNeurons=sortrows(upstreamNeurons,2,'descend');

end
