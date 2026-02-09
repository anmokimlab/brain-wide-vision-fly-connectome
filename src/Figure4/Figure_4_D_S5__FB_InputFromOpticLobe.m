%% 데이터 가져오기
clear all; close all; clc
load('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\Figure1_Overall\Right_Neurons_Thr0.mat')
load('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\Figure1_Overall\FAFB_NPI_Thr0.mat')


[RightFB_type,~,ic]=unique(RightFB_NPIs.type);

RightFB_type=table(RightFB_type,'VariableNames',{'type'});

for i=1:1:size(RightFB_type,1)
    idx=ic==i;
    RightFB_type.root_id{i}=RightFB_NPIs.root_id(idx);
end
%%
% opt = detectImportOptions('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_v3\CodexData\connections_no_threshold.csv');
% opt = setvartype(opt,'pre_root_id','int64');
% opt = setvartype(opt,'post_root_id','int64');
% FAFBConnections = readtable('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_v3\CodexData\connections_no_threshold.csv',opt);

opt = detectImportOptions('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\connections_no_threshold.csv');
opt = setvartype(opt,'pre_root_id','int64');
opt = setvartype(opt,'post_root_id','int64');
FAFBConnections = readtable('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\connections_no_threshold.csv',opt);

opt = detectImportOptions('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\classification.csv');
opt = setvartype(opt,'root_id','int64');
FAFBClassification = readtable('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\classification.csv',opt);

opt = detectImportOptions('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\consolidated_cell_types.csv');
opt = setvartype(opt,'root_id','int64');
FAFBConsolidated_type = readtable('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\consolidated_cell_types.csv',opt);

FAFB_Superclass=unique(FAFBClassification.super_class);
%% ME, LO, LOP MULTI 로 나눠서 해보자
Me_FB_idx=[4 22 25 26 29 30 31 32 33 68 69 70 71 72 73 74 75 76 78 81 82 83 84 85 86 87 93];
Lo_FB_idx=[2 10 11 12 13 14 15 16 17 18 19 21 34 35 36 37 38 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 95];
Lop_FB_idx=[3 7 8 9 28 62 63 64 65 66];
Multi_FB_idx=[1 5 6 23 24 27 39 60 61 67 77 79 80 88 89 90 91 92 94 96];


opticlobes={'AME_R','ME_R','LO_R','LOP_R','AME_L','ME_L','LO_L','LOP_L'};


FB_Me=RightFB_type(Me_FB_idx,:);
FB_Lo=RightFB_type(Lo_FB_idx,:);
FB_Lop=RightFB_type(Lop_FB_idx,:);
FB_Multi=RightFB_type(Multi_FB_idx,:);


%%% right neuron 만...
for i=1:1:size(FB_Lop,1)
    Wantrootids=FB_Lop.root_id{i};

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

    FB_Lop.inputsFromOpticlobes{i}=sum(opticlobesSynapses)/sum(upstreamNeurons(:,2))*100;

end

%% 최종적으로LOP 정리
inputsFromOpticlobes_LOP=[];

for i=1:1:size(FB_Lop,1)
    inputsFromOpticlobes_LOP=[inputsFromOpticlobes_LOP;FB_Lop.inputsFromOpticlobes{i}];

end
inputsFromOpticlobes_LoP_mean=mean(inputsFromOpticlobes_LOP,'omitnan');
inputsFromOpticlobes_LoP_std=std(inputsFromOpticlobes_LOP,'omitnan');

%% LO

%%% right neuron 만...
for i=1:1:size(FB_Lo,1)
    Wantrootids=FB_Lo.root_id{i};

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
    % if sum(opticlobesSynapses(:))==0
    %     Target_Table.inputsFromOpticlobes{i}=zeros(size(opticlobes));
    % else
    FB_Lo.inputsFromOpticlobes{i}=sum(opticlobesSynapses)/sum(upstreamNeurons(:,2))*100;
    % end
end

%% 최종적으로  LO
inputsFromOpticlobes_LO=[];

for i=1:1:size(FB_Lo,1)
    inputsFromOpticlobes_LO=[inputsFromOpticlobes_LO;FB_Lo.inputsFromOpticlobes{i}];

end
inputsFromOpticlobes_Lo_mean=mean(inputsFromOpticlobes_LO,'omitnan');
inputsFromOpticlobes_Lo_std=std(inputsFromOpticlobes_LO,'omitnan');

%% Me

%%% right neuron 만...
for i=1:1:size(FB_Me,1)

     Wantrootids=FB_Me.root_id{i};

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
    % if sum(opticlobesSynapses(:))==0
    % Target_Table.inputsFromOpticlobes{i}=zeros(size(opticlobes));
    % else
    FB_Me.inputsFromOpticlobes{i}=sum(opticlobesSynapses)/sum(upstreamNeurons(:,2))*100;
    % end
end

%% 최종적으로 Me
inputsFromOpticlobes_Me=[];

for i=1:1:size(FB_Me,1)
    inputsFromOpticlobes_Me=[inputsFromOpticlobes_Me;FB_Me.inputsFromOpticlobes{i}];

end
inputsFromOpticlobes_Me_mean=mean(inputsFromOpticlobes_Me,'omitnan');
inputsFromOpticlobes_Me_std=std(inputsFromOpticlobes_Me,'omitnan');

%% Multi

%%% right neuron 만...
for i=1:1:size(FB_Multi,1)
     Wantrootids=FB_Multi.root_id{i};


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
    % if sum(opticlobesSynapses(:))==0
    %     Target_Table.inputsFromOpticlobes{i}=zeros(size(opticlobes));
    % else
    % Target_Table.inputsFromOpticlobes{i}=sum(opticlobesSynapses)/sum(opticlobesSynapses(:))*100;
    FB_Multi.inputsFromOpticlobes{i}=sum(opticlobesSynapses)/sum(upstreamNeurons(:,2))*100;

    % end
end

%% 최종적으로 Multi
inputsFromOpticlobes_Multi=[];

for i=1:1:size(FB_Multi,1)
    inputsFromOpticlobes_Multi=[inputsFromOpticlobes_Multi;FB_Multi.inputsFromOpticlobes{i}];

end
inputsFromOpticlobes_Multi_mean=mean(inputsFromOpticlobes_Multi,'omitnan');
inputsFromOpticlobes_Multi_std=std(inputsFromOpticlobes_Multi,'omitnan');
%%
% 데이터 불러오기
% load('InputsFromOpticlobes.mat'); % 파일 이름 맞춰서

% 네 개의 매트릭스 로드
matrices = {inputsFromOpticlobes_Me, inputsFromOpticlobes_LO, inputsFromOpticlobes_LOP, inputsFromOpticlobes_Multi};
conditionNames = {'CB-Me','CB-LO', 'CB-LOP', 'CB-Multi'};

% NaN을 0으로 변환
for i = 1:length(matrices)
    matrices{i}(isnan(matrices{i})) = 0;
end

% Neuropil 개수 가정 (8개)
num_neuropils = size(matrices{1}, 2);

% 각 조건별 뉴로필 평균값 계산 (행: neuropil, 열: condition)
meanInputs = zeros(num_neuropils, length(matrices));
for i = 1:length(matrices)
    meanInputs(:, i) = mean(matrices{i}, 1)';  % Transpose해서 세로로 넣음
end

% Heatmap 그리기
figure;set(gcf,'Color','w')
h = heatmap(conditionNames, 1:num_neuropils, meanInputs, ...
    'Colormap', bone, ...
    'ColorbarVisible', 'on');

ylabel('Neuropils');
xlabel('Condition');
title('Mean Input Fractions per Neuropil (rotated view)');
h.ColorLimits = [0 max(meanInputs(:))];

% (선택) y축 라벨 이름 지정
h.YDisplayLabels = {'Me R','Me L','Lo R','Lo L','LoP R','LoP L'};
%%
% 네 개의 매트릭스 불러오기
matrices = {inputsFromOpticlobes_Me, inputsFromOpticlobes_LO, ...
            inputsFromOpticlobes_LOP, inputsFromOpticlobes_Multi};
conditionNames = {'CB-Me','CB-LO', 'CB-LOP', 'CB-Multi'};

% NaN을 0으로 변환
for i = 1:length(matrices)
    matrices{i}(isnan(matrices{i})) = 0;
end

% Neuropil 개수
num_neuropils = size(matrices{1}, 2);
neuropilLabels = {'Me R','Me L','Lo R','Lo L','LoP R','LoP L'};  % 필요시 수정

% 평균과 표준편차 계산
meanInputs = zeros(num_neuropils, length(matrices));
stdInputs = zeros(num_neuropils, length(matrices));
for i = 1:length(matrices)
    meanInputs(:, i) = mean(matrices{i}, 1)';     % 열 평균
    stdInputs(:, i) = std(matrices{i}, 0, 1)';    % 열 표준편차
end

% 시각화 (imagesc + text)
figure;
imagesc(meanInputs);
fullMap = gray(256);         % 기본 bone colormap
trimmedMap = fullMap(70:end, :);  % 앞 29단계 제거 (0~29번 색상 생략)

% trimmedMap을 colormap으로 설정
colormap(fullMap);
colorbar;
title('Mean ± STD of Inputs per Neuropil');
xlabel('Condition');
ylabel('Neuropils');

% 축 라벨 설정
set(gca, 'XTick', 1:length(conditionNames), 'XTickLabel', conditionNames);
set(gca, 'YTick', 1:num_neuropils, 'YTickLabel', neuropilLabels);

% 셀마다 텍스트 추가
for i = 1:num_neuropils
    for j = 1:length(matrices)
        text(j, i, sprintf('%.2f±%.2f', meanInputs(i,j), stdInputs(i,j)), ...
            'HorizontalAlignment', 'center', 'FontSize', 9, 'Color', 'black');
    end
end
print(gcf,'-depsc2','-vector','OpticLobeInputFBType.eps')
%% median
% 네 개의 매트릭스 불러오기
matrices = {inputsFromOpticlobes_Me, inputsFromOpticlobes_LO, ...
            inputsFromOpticlobes_LOP, inputsFromOpticlobes_Multi};
conditionNames = {'CB-Me','CB-LO', 'CB-LOP', 'CB-Multi'};

% NaN을 0으로 변환
for i = 1:length(matrices)
    matrices{i}(isnan(matrices{i})) = 0;
end

% Neuropil 개수
num_neuropils = size(matrices{1}, 2);
neuropilLabels = {'Me R','Me L','Lo R','Lo L','LoP R','LoP L'};  % 필요 시 수정

% 중앙값, Q1, Q3 계산
medInputs = zeros(num_neuropils, length(matrices));
q1Inputs = zeros(num_neuropils, length(matrices));
q3Inputs = zeros(num_neuropils, length(matrices));

for i = 1:length(matrices)
    mat = matrices{i};
    medInputs(:, i) = median(mat, 1)';
    q1Inputs(:, i) = prctile(mat, 25, 1)';
    q3Inputs(:, i) = prctile(mat, 75, 1)';
end

% 시각화 (imagesc + text)
figure;
imagesc(medInputs);

% bone colormap trimming
fullMap = bone(256);
trimmedMap = fullMap(100:end, :);  % 어두운 색 제거
colormap(trimmedMap);
colorbar;

title('Median [Q1–Q3] of Inputs per Neuropil');
xlabel('Condition');
ylabel('Neuropils');
set(gca, 'XTick', 1:length(conditionNames), 'XTickLabel', conditionNames);
set(gca, 'YTick', 1:num_neuropils, 'YTickLabel', neuropilLabels);

% 셀마다 텍스트 추가 ("2.00 [0.50–3.50]" 형식)
for i = 1:num_neuropils
    for j = 1:length(matrices)
        med = medInputs(i,j);
        q1 = q1Inputs(i,j);
        q3 = q3Inputs(i,j);
        text(j, i, sprintf('%.2f [%.2f–%.2f]', med, q1, q3), ...
            'HorizontalAlignment', 'center', 'FontSize', 9, 'Color', 'black');
    end
end

%%
% 데이터 불러오기
% load('InputsFromOpticlobes.mat'); % 파일 이름 맞춰서

% 네 개의 매트릭스 로드
matrices = {inputsFromOpticlobes_Me, inputsFromOpticlobes_LO, inputsFromOpticlobes_LOP, inputsFromOpticlobes_Multi};
conditionNames = {'CB-Me','CB-LO', 'CB-LOP', 'CB-Multi'};

% NaN을 0으로 변환
for i = 1:length(matrices)
    matrices{i}(isnan(matrices{i})) = 0;
end

% Neuropil 개수 가정 (8개)
num_neuropils = size(matrices{1}, 2);

% 각 조건별 뉴로필 평균값 계산 (행: neuropil, 열: condition)
medianInputs = zeros(num_neuropils, length(matrices));
for i = 1:length(matrices)
    medianInputs(:, i) = median(matrices{i}, 1)';  % Transpose해서 세로로 넣음
end

% Heatmap 그리기
figure;set(gcf,'Color','w')
h = heatmap(conditionNames, 1:num_neuropils, medianInputs, ...
    'Colormap', bone, ...
    'ColorbarVisible', 'on');

ylabel('Neuropils');
xlabel('Condition');
title('Mean Input Fractions per Neuropil (rotated view)');
h.ColorLimits = [0 max(medianInputs(:))];

% (선택) y축 라벨 이름 지정
h.YDisplayLabels = {'Me R','Me L','Lo R','Lo L','LoP R','LoP L'};
% print(gcf,'-depsc2','-vector','OpticLobeInputFBType.eps')
%%
figure;set(gcf,'Color','w')
for i=1:1:size(FB_Me,1)
    inputFromME_R(i)= FB_Me.inputsFromOpticlobes{i}(1);
    inputFromME_L(i)= FB_Me.inputsFromOpticlobes{i}(2);
    
end
y = [inputFromME_L' inputFromME_R'];

b = bar(FB_Me.type, y);
b(1).FaceColor = [0.4660 0.6740 0.1880];  % 첫째 열 색 (RGB 0~1)
b(2).FaceColor = [0 0.4470 0.7410];  % 둘째 열 색

set(gca, 'TickDir', 'out', 'Box', 'off');
ylim([0 55]);
ylabel('Overlap coefficient (Gaussian smoothed)');
title('input');
print(gcf, '-depsc2', '-vector', 'Input_From_Me_FB_ME_bar.eps');

%%
figure;set(gcf,'Color','w')
for i=1:1:size(FB_Lo,1)
    inputFromLO_R(i)= FB_Lo.inputsFromOpticlobes{i}(3);
    inputFromLO_L(i)= FB_Lo.inputsFromOpticlobes{i}(4);
    
end
y = [inputFromLO_L' inputFromLO_R'];

b = bar(FB_Lo.type, y);
b(1).FaceColor = [0.4660 0.6740 0.1880];  % 첫째 열 색 (RGB 0~1)
b(2).FaceColor = [0 0.4470 0.7410];  % 둘째 열 색

set(gca, 'TickDir', 'out', 'Box', 'off');
ylim([0 55]);
ylabel('Overlap coefficient (Gaussian smoothed)');
title('input');
print(gcf, '-depsc2', '-vector', 'Input_From_Lo_FB_LO_bar.eps');
%%
figure;set(gcf,'Color','w')
for i=1:1:size(FB_Lop,1)
    inputFromLOP_R(i)= FB_Lop.inputsFromOpticlobes{i}(5);
    inputFromLOP_L(i)= FB_Lop.inputsFromOpticlobes{i}(6);
    
end
y = [inputFromLOP_L' inputFromLOP_R'];

b = bar(FB_Lop.type, y);
b(1).FaceColor = [0.4660 0.6740 0.1880];  % 첫째 열 색 (RGB 0~1)
b(2).FaceColor = [0 0.4470 0.7410];  % 둘째 열 색

set(gca, 'TickDir', 'out', 'Box', 'off');
ylim([0 55]);
ylabel('Overlap coefficient (Gaussian smoothed)');
title('input ');
print(gcf, '-depsc2', '-vector', 'Input_From_Lo_FB_LOP_bar.eps');
