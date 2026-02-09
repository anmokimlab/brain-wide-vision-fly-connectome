%% LC9 PVLP004 연결성 조사... 밑 자체 연결
clear all; close all; clc
opt = detectImportOptions('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\consolidated_cell_types.csv');
opt = setvartype(opt,'root_id','int64');
FAFBConsolidated_type = readtable('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\consolidated_cell_types.csv',opt);


% opt = detectImportOptions('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_v3\CodexData\connections.csv');
% opt = setvartype(opt,'pre_root_id','int64');
% opt = setvartype(opt,'post_root_id','int64');
% FAFBConnections = readtable('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_v3\CodexData\connections.csv',opt);

opt = detectImportOptions('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\connections_no_threshold.csv');
opt = setvartype(opt,'pre_root_id','int64');
opt = setvartype(opt,'post_root_id','int64');
FAFBConnections = readtable('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\connections_no_threshold.csv',opt);

opt = detectImportOptions('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\classification.csv');
opt = setvartype(opt,'root_id','int64');
FAFBClassification = readtable('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\classification.csv',opt);

% === [수정됨] 뉴로필 정의 ===
% Central/Optic을 명확히 구분하기 위해 정의를 추가합니다.
FAFBNeuropils = unique(FAFBConnections.neuropil);
FAFBNeuropil_OpticLobeRight = {'LA_R';'AME_R';'ME_R';'LO_R';'LOP_R'};
FAFBNeuropil_OpticLobeLeft  = {'LA_L';'AME_L';'ME_L';'LO_L';'LOP_L'};
% Central은 Optic Lobe와 UNASGD를 제외한 영역으로 정의
FAFBNeuropil_Central = FAFBNeuropils(~ismember(FAFBNeuropils,FAFBNeuropil_OpticLobeRight));
FAFBNeuropil_Central = FAFBNeuropil_Central(~ismember(FAFBNeuropil_Central,FAFBNeuropil_OpticLobeLeft));
FAFBNeuropil_Central = FAFBNeuropil_Central(~ismember(FAFBNeuropil_Central,'UNASGD'));
    
load('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\Figure1_Overall\FAFB_NPI_Thr0.mat')

LC9_root_ids=FAFBNPIs.root_id(strcmp(FAFBNPIs.type,'LC9'));

LC9_R_root_ids=[];
LC9_L_root_ids=[];
PVLP004_R_root_ids=[];
PVLP004_L_root_ids=[];

for i=1:1:size(LC9_root_ids,1)
    current_root_id=LC9_root_ids(i);
    idx_classification=FAFBClassification.root_id==current_root_id;
    if strcmp(FAFBClassification.side(idx_classification),'right')
        LC9_R_root_ids=[LC9_R_root_ids; current_root_id];
    elseif strcmp(FAFBClassification.side(idx_classification),'left')
        LC9_L_root_ids=[LC9_L_root_ids; current_root_id];
    end
end

%%
LC9_interConnection_OpticLobe=zeros(size(LC9_R_root_ids,1));
LC9_interConnection_Central=zeros(size(LC9_R_root_ids));


for i=1:1:size(LC9_R_root_ids,1)
    current_root_id=LC9_R_root_ids(i);

    OutConections=FAFBConnections(ismember(FAFBConnections.pre_root_id,current_root_id),:);
    OutConections_LC9=OutConections(ismember(OutConections.post_root_id,LC9_R_root_ids),:);

    OutConections_LC9 = sortrows(OutConections_LC9,"post_root_id","ascend");
    OutConections_LC9_Optic=OutConections_LC9(ismember(OutConections_LC9.neuropil,FAFBNeuropil_OpticLobeRight),:);
    OutConections_LC9_Central=OutConections_LC9(ismember(OutConections_LC9.neuropil,FAFBNeuropil_Central),:);
    
    % 1열과 2열의 조합 찾기
    [uniquePairs, ~, idx] = unique(OutConections_LC9_Optic(:, 1:2), 'rows');
    
    % 각 그룹의 4열 합산
    summedValues = accumarray(idx, table2array(OutConections_LC9_Optic(:,4)));
    
    for j=1:1:size(uniquePairs,1)
        idx_out=find(LC9_R_root_ids==uniquePairs.post_root_id(j));
        LC9_interConnection_OpticLobe(i,idx_out)=summedValues(j);
    end

       % 1열과 2열의 조합 찾기
    [uniquePairs, ~, idx] = unique(OutConections_LC9_Central(:, 1:2), 'rows');
    
    % 각 그룹의 4열 합산
    summedValues = accumarray(idx, table2array(OutConections_LC9_Central(:,4)));
    
    for j=1:1:size(uniquePairs,1)
        idx_out=find(LC9_R_root_ids==uniquePairs.post_root_id(j));
        LC9_interConnection_Central(i,idx_out)=summedValues(j);
    end
end
LC9_interConnection=LC9_interConnection_Central+LC9_interConnection_OpticLobe;
%%
figure(1);set(gcf,'Color','w')
imagesc(LC9_interConnection)
grid on;
set(gca,'Box','off','TickDir','out');

axis('square')
clim([0 20])
%%
% 1) 특징 행렬 (예: in/out 연결 프로파일을 붙인 X)
X = [LC9_interConnection LC9_interConnection'];   % n x 2n

% 2) 샘플 간 유클리드 거리
Y = pdist(X, 'euclidean');                        % == pdist(X)와 동일

% 3) 계층 군집 (average linkage; 유클리드 기반)
Z = linkage(Y, 'average');

% 4) 잎 순서 최적화도 동일한 거리로
order = optimalleaforder(Z, Y);

%%
figure(2);set(gcf,'Color','w')

imagesc(LC9_interConnection(order, order))
set(gca,'Box','off','TickDir','out');

axis('square')
title('LC9 interConnection All')
print(gcf, '-depsc2', '-vector', 'figure3_LC9 interConnection All.eps');
clim([0 20])

figure(3);set(gcf,'Color','w')

imagesc(LC9_interConnection_OpticLobe(order, order))
set(gca,'Box','off','TickDir','out');

axis('square')
title('LC9 interConnection optic')
print(gcf, '-depsc2', '-vector', 'figure3_LC9 interConnection Optic.eps');
clim([0 20])

figure(4);set(gcf,'Color','w')

imagesc(LC9_interConnection_Central(order, order))
set(gca,'Box','off','TickDir','out');

axis('square')
title('LC9 interConnection central')
print(gcf, '-depsc2', '-vector', 'figure3_LC9 interConnection Central.eps');
clim([0 20])


%%
% Z: linkage 결과
% order: 정렬 순서
% figure와 dendrogram
figure(); set(gcf,'Color','w')
[H, T, perm] = dendrogram(Z, 0, 'Reorder', order, colorThreshold=41);
set(gca,'TickDir','out', 'Box','off')
% print(gcf, '-depsc2', '-vector', 'figure3_LC9 interConnection Central_dendrogram.eps');


%%
% linkage 결과 Z와 reorder 순서 order가 있다고 가정
colorThreshold =41; % default와 동일
T = cluster(Z, 'cutoff', colorThreshold, 'criterion', 'distance');  % 색상에 따라 클러스터 나누기

% 클러스터 번호 확인 및 색상별 인덱스 저장
unique_clusters = unique(T);
cluster_indices = cell(length(unique_clusters), 1);
for i = 1:length(unique_clusters)
    cluster_indices{i,1} = find(T == unique_clusters(i));
    cluster_indices{i,2}=LC9_R_root_ids(cluster_indices{i,1});
end
