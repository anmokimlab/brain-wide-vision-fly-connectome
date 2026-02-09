%% 데이터 가져오기
clear all; close all; clc;

opt = detectImportOptions('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\Figure2_FF\Clustering_Result\Leiden_RightFF_OUTPUT_Central_No_VPN_Thr0_Resol_3.csv');
opt = setvartype(opt,'root_id','int64');
ClusteringResult=readtable('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\Figure2_FF\Clustering_Result\Leiden_RightFF_OUTPUT_Central_No_VPN_Thr0_Resol_3.csv',opt);

opt = detectImportOptions('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\connections_no_threshold.csv');
opt = setvartype(opt,'pre_root_id','int64');
opt = setvartype(opt,'post_root_id','int64');
FAFBConnections = readtable('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\connections_no_threshold.csv',opt);

FAFBNeuropils=unique(FAFBConnections.neuropil);
FAFBNeuropil_OpticLobeRight={'AME_R','ME_R','LO_R','LOP_R','LA_R'};
FAFBNeuropil_OpticLobeLeft={'AME_L','ME_L','LO_L','LOP_L','LA_L'};
FAFBNeuropil_Central=FAFBNeuropils(~ismember(FAFBNeuropils,FAFBNeuropil_OpticLobeRight));
FAFBNeuropil_Central=FAFBNeuropil_Central(~ismember(FAFBNeuropil_Central,FAFBNeuropil_OpticLobeLeft));
FAFBNeuropil_Central=FAFBNeuropil_Central(~ismember(FAFBNeuropil_Central,'UNASGD'));

opt = detectImportOptions('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\consolidated_cell_types.csv');
opt = setvartype(opt,'root_id','int64');
FAFB_consolidated_cell_types = readtable('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\consolidated_cell_types.csv',opt);

opt = detectImportOptions('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\classification.csv');
opt = setvartype(opt,'root_id','int64');
FAFBClassification = readtable('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\classification.csv',opt);

opt = detectImportOptions('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\Figure2_FF\Clustering_Result\Leiden_RightFF_OUTPUT_Intercluster_Connectivity_Central_No_VPN_Thr0_Resol_3.csv');
BipartiteConnectivity = readmatrix('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\Figure2_FF\Clustering_Result\Leiden_RightFF_OUTPUT_Intercluster_Connectivity_Central_No_VPN_Thr0_Resol_3.csv',opt);
BipartiteConnectivity = reshape(BipartiteConnectivity,[max(ClusteringResult.Cluster)+1,max(ClusteringResult.Cluster)+1]);

opt = detectImportOptions('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\Figure2_FF\Clustering_Result\dendrogram_raw_Central_No_VPN_Thr0_Resol_3.csv');
RightFF_Dendrogram = readmatrix('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\Figure2_FF\Clustering_Result\dendrogram_raw_Central_No_VPN_Thr0_Resol_3.csv',opt);
RightFF_Dendrogram(:,4)=[];
RightFF_Dendrogram(:,1)=RightFF_Dendrogram(:,1)+1;
RightFF_Dendrogram(:,2)=RightFF_Dendrogram(:,2)+1;

load('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\Figure2_FF\Clustering_Result\Data\Pre_Post_neurons_FF_Opticlobe_Central_No_VPN_Thr0.mat')
ClusteringResult_PostNeurons=readtable('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\Figure2_FF\Clustering_Result\Leiden_Post_RightFF_Central_No_VPN_Thr0.csv');
ClusteringResult_PostNeurons(1,:) = [];
ClusteringResult_PostNeurons = removevars(ClusteringResult_PostNeurons, "Var1");
ClusteringResult_PostNeurons.Properties.VariableNames(1) = "root_id";
ClusteringResult_PostNeurons.Properties.VariableNames(2) = "Cluster";
ClusteringResult_PostNeurons.root_id=Post_Want;
load('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\Figure2_FF\Clustering_Result/Data\RightFF_OUTPUT_matrix_Central_No_VPN_Thr0.mat')
% WantMatrix_In=WantMatrix_In';
load('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\Figure1_Overall\Right_Neurons_Thr0.mat')

%% RightFF 정리
ClusteringValues={};
for i=min(ClusteringResult.Cluster):1:max(ClusteringResult.Cluster)
    ClusteringValues{i+1,1}=num2str(i);

end

ClusteringValueCounts={};
ClusteringValueCounts(:,1)=ClusteringValues;
for i=1:1:size(ClusteringValueCounts,1)
    ClusteringValueCounts{i,2}=sum(ClusteringResult.Cluster==str2double(ClusteringValueCounts{i,1}));
    ClusteringValueCounts{i,3}=(ClusteringResult.Cluster==str2double(ClusteringValueCounts{i,1}));
    ClusteringValueCounts{i,5}=ClusteringResult.root_id(ClusteringValueCounts{i,3});
end
% ClusteringValueCounts=[{'root',size(ClusteringResult,1),ones(size(ClusteringResult,1),1),[],ClusteringResult.root_id};ClusteringValueCounts];

for i=1:1:size(ClusteringValueCounts,1)
    AllType=ClusteringResult.type(ClusteringValueCounts{i,3});
    UniqueType=unique(AllType);
    for j=1:1:size(UniqueType,1)
        UniqueType{j,2}=sum(strcmp(AllType,UniqueType{j,1}));
        UniqueType{j,3}=ClusteringValueCounts{i,5}(strcmp(AllType,UniqueType{j,1}));
        UniqueType{j,4}=sum(strcmp(ClusteringResult.type,UniqueType{j,1}));
    end

    if ~isempty(UniqueType)
        UniqueType=sortrows(UniqueType,2,'descend');
        if isempty(UniqueType{1,1})
            if size(UniqueType,1)>=3
                ClusteringValueCounts{i,4}=[UniqueType{2,1} ', ' UniqueType{3,1}];
            elseif size(UniqueType,1)==1
                ClusteringValueCounts{i,4}='UnLabeled';
            else
                ClusteringValueCounts{i,4}=[UniqueType{2,1}];
            end
        else
            if size(UniqueType,1)>=2
                ClusteringValueCounts{i,4}=[UniqueType{1,1} ', ' UniqueType{2,1}];
            else
                ClusteringValueCounts{i,4}=[UniqueType{1,1}];
            end

        end
    else
        ClusteringValueCounts{i,4}='UnLabeled';
    end
    ClusteringValueCounts{i,6}=UniqueType;
end

%% PostNeuron 정리 정리
ClusteringValues={};
for i=min(ClusteringResult_PostNeurons.Cluster):1:max(ClusteringResult_PostNeurons.Cluster)
    ClusteringValues{i+1,1}=num2str(i);

end

ClusteringValueCounts_PostNeurons={};
ClusteringValueCounts_PostNeurons(:,1)=ClusteringValues;
for i=1:1:size(ClusteringValueCounts_PostNeurons,1)
    ClusteringValueCounts_PostNeurons{i,2}=sum(ClusteringResult_PostNeurons.Cluster==str2double(ClusteringValueCounts_PostNeurons{i,1}));
    ClusteringValueCounts_PostNeurons{i,3}=(ClusteringResult_PostNeurons.Cluster==str2double(ClusteringValueCounts_PostNeurons{i,1}));
    ClusteringValueCounts_PostNeurons{i,5}=ClusteringResult_PostNeurons.root_id(ClusteringValueCounts_PostNeurons{i,3});
end

for i=1:1:size(ClusteringResult_PostNeurons,1)
    idx_Consol=find(FAFB_consolidated_cell_types.root_id==ClusteringResult_PostNeurons.root_id(i));
    if ~isempty(idx_Consol)
        ClusteringResult_PostNeurons.type{i}=FAFB_consolidated_cell_types.primary_type{idx_Consol};
    else
        ClusteringResult_PostNeurons.type{i}='';
    end
end

for i=1:1:size(ClusteringValueCounts_PostNeurons,1)
    AllType=ClusteringResult_PostNeurons.type(ClusteringValueCounts_PostNeurons{i,3});
    UniqueType=unique(AllType);
    current_FF_Cluster=ClusteringValueCounts{i,6}  ;
    for j=1:1:size(UniqueType,1)
        UniqueType{j,2}=sum(strcmp(AllType,UniqueType{j,1}));
        UniqueType{j,3}=ClusteringValueCounts_PostNeurons{i,5}(strcmp(AllType,UniqueType{j,1}));
        UniqueType{j,4}=sum(FAFBConnections.syn_count(ismember(FAFBConnections.post_root_id,UniqueType{j,3})&ismember(FAFBConnections.pre_root_id,ClusteringValueCounts{i, 5})));
        % mergingNumber=0;
        % mergedNeurons={};
        % for k=1:1:size(current_FF_Cluster,1)
        %
        %     current_FF_Cluster_root_id=current_FF_Cluster{k, 3};
        %     idx_pre=ismember(FAFBConnections.pre_root_id,current_FF_Cluster_root_id);
        %     idx_post=ismember(FAFBConnections.post_root_id,UniqueType{j,3});
        %     % if sum(idx_pre&idx_post)>0
        %     %     mergingNumber=mergingNumber+1;
        %     %     mergedNeurons{mergingNumber,1}= current_FF_Cluster{k,1};
        %     % end
        % end
        % UniqueType{j,5}=mergingNumber;
        % UniqueType{j,6}=mergedNeurons;
    end
    if ~isempty(UniqueType)
        UniqueType=sortrows(UniqueType,4,'descend');
        if isempty(UniqueType{1,1})
            if size(UniqueType,1)>=3
                ClusteringValueCounts_PostNeurons{i,4}=[UniqueType{2,1} ', ' UniqueType{3,1}];
            elseif size(UniqueType,1)==1
                ClusteringValueCounts_PostNeurons{i,4}='UnLabeled';
            else
                ClusteringValueCounts_PostNeurons{i,4}=[UniqueType{2,1}];
            end
        else
            if size(UniqueType,1)>=2
                ClusteringValueCounts_PostNeurons{i,4}=[UniqueType{1,1} ', ' UniqueType{2,1}];
            else
                ClusteringValueCounts_PostNeurons{i,4}=[UniqueType{1,1}];
            end

        end
    else
        ClusteringValueCounts_PostNeurons{i,4}='UnLabeled';
    end
    ClusteringValueCounts_PostNeurons{i,6}=UniqueType;
end

%% 직접 label 매기기  5 % 까지
for i=1:1:size(ClusteringValueCounts,1)
    ClusteringValueCounts{i,7}=size(ClusteringValueCounts{i,5},1)*0.05;
end
for i=1:1:size(ClusteringValueCounts_PostNeurons,1)
    ClusteringValueCounts_PostNeurons{i,7}=size(ClusteringValueCounts_PostNeurons{i,5},1)*0.05;
end

RightFF_labels={};
Post_RightFF_labels={};
for i=1:1:size(ClusteringValueCounts,1)
    RightFF_labels_temp=[];
    for j=1:1:size(ClusteringValueCounts{i, 6},1)
        if ClusteringValueCounts{i, 6}{j, 2}>ClusteringValueCounts{i,7}
            RightFF_labels_temp=[RightFF_labels_temp ', ' ClusteringValueCounts{i, 6}{j, 1}];
        end
    end
    RightFF_labels_temp(1:2)=[];
    RightFF_labels{i,1}=RightFF_labels_temp;

    if i<=size(ClusteringValueCounts_PostNeurons,1)
        Post_RightFF_labels{i,1}=ClusteringValueCounts_PostNeurons{i,4};
    else
        Post_RightFF_labels{i,1}='';
    end
end

%% figure 3 bipartite graph connectivity
figure(3);
set(gcf,'Color','w')
imagesc(log10(1+BipartiteConnectivity))
set(gca,'XTick',1:1:size(Post_RightFF_labels,1),'XTickLabel',Post_RightFF_labels,'YTick',1:1:size(RightFF_labels,1),'YTickLabel',RightFF_labels)
colors=(colormap('abyss'));
colormap(colors)
grid on
%%
% temp_RightFF_labels=RightFF_labels(1:1:29);
temp_RightFF_labels=RightFF_labels;

figure(4);
set(gcf,'Color','w')
[H, ~, DendrogramOrder]=dendrogram(RightFF_Dendrogram,0,'ColorThreshold',1,'Labels',temp_RightFF_labels);
set(H,'LineWidth',2)
set(gca,'TickDir','out')
ylim([0 2])
print(gcf,'-depsc2','-vector','figure2_dendrogram.eps')
%%
% RightFF_Dendrogram 은 linkage matrix
leafOrder = get_leaf_order_from_linkage(RightFF_Dendrogram);

% 덴드로그램 시각화
figure(5);
set(gcf,'Color','w')
[H,~,order] = dendrogram(RightFF_Dendrogram, 0, ...
    'Reorder', leafOrder, 'ColorThreshold', 1, 'Labels', temp_RightFF_labels);
set(H,'LineWidth',2);
set(gca,'TickDir','out');
% ylim([0 10]);
% print(gcf,'-depsc2','-vector','figure2_dendrogram.eps')


%%
figure(6);
set(gcf,'Color','w')
imagesc(log10(1+BipartiteConnectivity(leafOrder,leafOrder)))
% imagesc(BipartiteConnectivity(DendrogramOrder,DendrogramOrder))
% xlim([0.5 29.5])
% ylim([0.5 29.5])
% set(gca,'XTick',1:1:size(Post_RightFF_labels,1),'YTick',1:1:size(RightFF_labels,1));
set(gca,'XTick',1:1:size(Post_RightFF_labels,1),'XTickLabel',Post_RightFF_labels(leafOrder),'YTick',1:1:size(RightFF_labels,1),'YTickLabel',RightFF_labels(leafOrder));
set(gca,'TickDir','out','Box','off')
colors=(colormap('sky'));
colormap(colors)
axis square
grid on
% print(gcf,'-depsc2','-vector','figure2_ClusteringResult.eps')

%%
% RightFF_Dendrogram 은 linkage(Z) 결과라고 가정
threshold = 1;

% 특정 거리(threshold) 기준으로 클러스터 번호 할당
T = cluster(RightFF_Dendrogram, 'cutoff', threshold, 'criterion', 'distance');

% 고유 클러스터들 찾기
unique_clusters = unique(T);

% 각 클러스터에 포함된 리프 인덱스 구하기
cluster_groups = cell(length(unique_clusters), 1);
for i = 1:length(unique_clusters)
    cluster_id = unique_clusters(i);
    cluster_groups{i} = find(T == cluster_id);
end

% 결과 출력
for i = 1:length(cluster_groups)
    % fprintf('Cluster %d: ', i);
    temp=ClusteringValueCounts(cluster_groups{i},5);
    cluster_groups{i,2}=RightFF_labels(cluster_groups{i}) ;
    cluster_groups{i,3}=vertcat(temp{:}) ;
end

% cluster_groups(11:end,:)=[];
%%
opt = detectImportOptions('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\classification.csv');
opt = setvartype(opt,'root_id','int64');
FAFBClassification = readtable('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\classification.csv',opt);

Central_root_id=FAFBClassification.root_id(strcmp(FAFBClassification.super_class,'central'));
Optic_root_id=FAFBClassification.root_id(strcmp(FAFBClassification.super_class,'optic'));
VPN_root_id=FAFBClassification.root_id(strcmp(FAFBClassification.super_class,'visual_projection'));
VCN_root_id=FAFBClassification.root_id(strcmp(FAFBClassification.super_class,'visual_centrifugal'));

%%
for i=1:1:size(cluster_groups,1)
    %%% 타입 계산
    cluster_groups_cell_types={};
    current_root_ids=cluster_groups{i,3};
    for j=1:1:size(current_root_ids,1)
        idx=FAFB_consolidated_cell_types.root_id==current_root_ids(j);
        cluster_groups_cell_types{j,1}=FAFB_consolidated_cell_types.primary_type{idx};
    end
    [unique_cluster_groups_cell_types,~,ic]=unique(cluster_groups_cell_types);
    for j=1:1:size(unique_cluster_groups_cell_types,1)
        idx=ic==j;
        unique_cluster_groups_cell_types{j,2}=sum(idx);
        unique_cluster_groups_cell_types{j,3}=current_root_ids(idx);
    end
    unique_cluster_groups_cell_types=sortrows(unique_cluster_groups_cell_types,2,'descend');
    cluster_groups{i,4}=unique_cluster_groups_cell_types;


    %%% 출력 뭐가 가장 많이 받는지 계산해야지 시냅스로?
    post_idx=ismember(FAFBConnections.pre_root_id,current_root_ids);
    post_Connections=FAFBConnections(post_idx,:);
    % 같은기준 매트릭스 계산이랑
    OpticR=ismember(post_Connections.neuropil,{'LA_R', 'ME_R','AME_R', 'LO_R','LOP_R','LA_L', 'ME_L','AME_L', 'LO_L','LOP_L'});
    post_Connections(OpticR,:)=[];
    post_Connections(ismember(post_Connections.post_root_id,VPN_root_id),:)=[];
    %%%%%%%%%%%%%%%%%
    % post_Connections(ismember(post_Connections.post_root_id,Optic_root_id),:)=[];
    %%%%%%%%%%%%%%%%% 여기 봐야함....
    post_Neurons=unique(post_Connections.post_root_id);
    post_Neurons_cell_types={};
    for j=1:1:size(post_Neurons,1)
        idx=FAFB_consolidated_cell_types.root_id==post_Neurons(j);
        if any(idx)
            post_Neurons_cell_types{j,1}=FAFB_consolidated_cell_types.primary_type{idx};
        else
            post_Neurons_cell_types{j,1}='Others';
        end
    end
    [unique_post_Neurons_cell_types,~,ic]=unique(post_Neurons_cell_types);
    for j=1:1:size(unique_post_Neurons_cell_types,1)
        idx=ic==j;
        unique_post_Neurons_cell_types{j,2}=sum(idx);
        unique_post_Neurons_cell_types{j,3}=post_Neurons(idx);
        unique_post_Neurons_cell_types{j,4}=sum(post_Connections.syn_count(ismember(post_Connections.post_root_id,unique_post_Neurons_cell_types{j,3})));
    end
    unique_post_Neurons_cell_types=sortrows(unique_post_Neurons_cell_types,4,'descend');
    cluster_groups{i,6}=unique_post_Neurons_cell_types;
    %%% 출력 시냅스...
    [Post_neuropils,~,ic]=unique(post_Connections.neuropil);

    for j=1:1:size(Post_neuropils,1)
        idx=ic==j;
        Post_neuropils{j,2}=sum(post_Connections.syn_count(idx));
    end
    Post_neuropils=sortrows(Post_neuropils,2,'descend');
    cluster_groups{i,5}=Post_neuropils;
end
%% Dendrogram order 순서대로
ClusterSharednessMatrix=zeros(size(leafOrder,1));

for i=1:1:size(ClusterSharednessMatrix,1)
    Cluster1_root_id=ClusteringValueCounts{leafOrder(i),5};
    idx_1=ismember(ClusteringResult.root_id,Cluster1_root_id);
    % Matrix_1=WantMatrix_Out(idx_1,:)>0;
    Matrix_1=WantMatrix_Out(idx_1,:);

    for j = i+1:size(ClusterSharednessMatrix,2)

        Cluster2_root_id=ClusteringValueCounts{leafOrder(j),5};
        idx_2=ismember(ClusteringResult.root_id,Cluster2_root_id);
        % Matrix_2=WantMatrix_Out(idx_2,:)>0;
        Matrix_2=WantMatrix_Out(idx_2,:);
        % 
        % intersection = Matrix_1 * Matrix_2';                 % (200×10), 교집합 크기
        % union_size = bsxfun(@plus, sum(Matrix_1,2), sum(Matrix_2,2)') - intersection;  % (200×10), 합집합 크기
        % jaccard_mat = intersection ./ union_size;
        % jaccard_mat(union_size == 0) = NaN;  % union이 0인 경우 NaN 처리
        % ClusterSharednessMatrix(i,j)=mean(jaccard_mat(:),'omitnan');

        % 초기화
        numerator = zeros(size(Matrix_1,1), size(Matrix_2,1));
        denominator = zeros(size(Matrix_1,1), size(Matrix_2,1));

        % 각 column (postsynaptic neuron)에 대해 min/max 계산 후 누적
        for k = 1:size(Matrix_1,2)  % 각 postsynaptic 뉴런마다
            a_col = Matrix_1(:,k);   % A 클러스터: 200x1
            b_col = Matrix_2(:,k)';  % B 클러스터: 1x10

            min_mat = min(a_col, b_col);  % (200x10) – 각 (i,j)에 대해 min(a_ik, b_jk)
            max_mat = max(a_col, b_col);  % (200x10)

            numerator = numerator + min_mat;
            denominator = denominator + max_mat;
        end

        % Weighted Jaccard matrix
        weighted_jaccard_mat = numerator ./ denominator;
        weighted_jaccard_mat(denominator == 0) = NaN;  % 정의 안 되는 경우

        % 클러스터 간 공유도 = 평균 유사도
        ClusterSharednessMatrix(i,j) = mean(weighted_jaccard_mat(:), 'omitnan');
        ClusterSharednessMatrix(j,i) = mean(weighted_jaccard_mat(:), 'omitnan');

    end
end
%% Dendrogram order 순서대로

for i=1:1:size(ClusterSharednessMatrix,1)
    Cluster1_root_id=ClusteringValueCounts{leafOrder(i),5};
    idx_1=ismember(ClusteringResult.root_id,Cluster1_root_id);
    % Matrix_1=WantMatrix_Out(idx_1,:)>0;
    Matrix_1=WantMatrix_Out(idx_1,:);

 
        % 초기화
        numerator = zeros(size(Matrix_1,1), size(Matrix_1,1));
        denominator = zeros(size(Matrix_1,1), size(Matrix_1,1));

        % 각 column (postsynaptic neuron)에 대해 min/max 계산 후 누적
        for k = 1:size(Matrix_1,2)  % 각 postsynaptic 뉴런마다
            a_col = Matrix_1(:,k);   % A 클러스터: 200x1
            b_col = Matrix_1(:,k)';  % B 클러스터: 1x10

            min_mat = min(a_col, b_col);  % (200x10) – 각 (i,j)에 대해 min(a_ik, b_jk)
            max_mat = max(a_col, b_col);  % (200x10)

            numerator = numerator + min_mat;
            denominator = denominator + max_mat;
        end
        % Weighted Jaccard matrix
        weighted_jaccard_mat = numerator ./ denominator;
        weighted_jaccard_mat(denominator == 0) = NaN;  % 정의 안 되는 경우

        n_wj = size(weighted_jaccard_mat,1);
        weighted_jaccard_mat(1:n_wj+1:end) = NaN;  % diag를 NaN으로

        % 클러스터 간 공유도 = 평균 유사도
        ClusterSharednessMatrix(i,i) = mean(weighted_jaccard_mat(:), 'omitnan');
end

%%
RightFF_labels_Shared=RightFF_labels(leafOrder);

figure(7);
set(gcf,'Color','w')
imagesc(ClusterSharednessMatrix)
% set(gca,'XTick',1:1:size(Post_RightFF_labels,1),'YTick',1:1:size(RightFF_labels,1));
set(gca,'YTick',1:1:size(RightFF_labels,1),'YTickLabel',RightFF_labels(leafOrder));
% set(gca,'XTick',1:1:size(Post_RightFF_labels,1),'XTickLabel',Post_RightFF_labels(DendrogramOrder),'YTick',1:1:size(RightFF_labels,1),'YTickLabel',RightFF_labels(DendrogramOrder));

set(gca,'TickDir','out','Box','off')
colors=(colormap('sky'));
colormap(colors)
axis square
grid on
clim([0 0.05]);
title('ClusterSharedness')


%% MetaClusterDendrogram order 순서대로
MetaClusterOrder{1}=[21;24;28];
MetaClusterOrder{2}=[15;17;25;29];
MetaClusterOrder{3}=[9;18;19;26;27;30];
MetaClusterOrder{4}=[3;5;6;10;23];
MetaClusterOrder{5}=[4;7;12;13;14;22];
MetaClusterOrder{6}=[1;2;8;20];
MetaClusterOrder{7}=[11;16;31;32];

MetaClusterSharednessMatrix=zeros(size(MetaClusterOrder,2));

for i=1:1:size(MetaClusterSharednessMatrix,1)
    Cluster1_root_id=ClusteringValueCounts(MetaClusterOrder{i},5);
    Cluster1_root_id = vertcat(Cluster1_root_id{:});

    idx_1=ismember(ClusteringResult.root_id,Cluster1_root_id);
    % Matrix_1=WantMatrix_Out(idx_1,:)>0;
    Matrix_1=WantMatrix_Out(idx_1,:);

    for j=i+1:1:size(MetaClusterSharednessMatrix,2)
        Cluster2_root_id=ClusteringValueCounts(MetaClusterOrder{j},5);
        Cluster2_root_id = vertcat(Cluster2_root_id{:});

        idx_2=ismember(ClusteringResult.root_id,Cluster2_root_id);
        % Matrix_2=WantMatrix_Out(idx_2,:)>0;
        Matrix_2=WantMatrix_Out(idx_2,:);

        % intersection = Matrix_1 * Matrix_2';                 % (200×10), 교집합 크기
        % union_size = bsxfun(@plus, sum(Matrix_1,2), sum(Matrix_2,2)') - intersection;  % (200×10), 합집합 크기
        % jaccard_mat = intersection ./ union_size;
        % % jaccard_mat(union_size == 0) = NaN;  % union이 0인 경우 NaN 처리
        % MetaClusterSharednessMatrix(i,j)=mean(jaccard_mat(:));

        % 초기화
        numerator = zeros(size(Matrix_1,1), size(Matrix_2,1));
        denominator = zeros(size(Matrix_1,1), size(Matrix_2,1));

        % 각 column (postsynaptic neuron)에 대해 min/max 계산 후 누적
        for k = 1:size(Matrix_1,2)  % 각 postsynaptic 뉴런마다
            a_col = Matrix_1(:,k);   % A 클러스터: 200x1
            b_col = Matrix_2(:,k)';  % B 클러스터: 1x10

            min_mat = min(a_col, b_col);  % (200x10) – 각 (i,j)에 대해 min(a_ik, b_jk)
            max_mat = max(a_col, b_col);  % (200x10)

            numerator = numerator + min_mat;
            denominator = denominator + max_mat;
        end

        % Weighted Jaccard matrix
        weighted_jaccard_mat = numerator ./ denominator;
        weighted_jaccard_mat(denominator == 0) = NaN;  % 정의 안 되는 경우

        % 클러스터 간 공유도 = 평균 유사도
        MetaClusterSharednessMatrix(i,j) = mean(weighted_jaccard_mat(:), 'omitnan');
        MetaClusterSharednessMatrix(j,i) = mean(weighted_jaccard_mat(:), 'omitnan');

    end
end
%%
for i=1:1:size(MetaClusterSharednessMatrix,1)
    Cluster1_root_id=ClusteringValueCounts(MetaClusterOrder{i},5);
    Cluster1_root_id = vertcat(Cluster1_root_id{:});

    idx_1=ismember(ClusteringResult.root_id,Cluster1_root_id);
    % Matrix_1=WantMatrix_Out(idx_1,:)>0;
    Matrix_1=WantMatrix_Out(idx_1,:);

    % 초기화
    numerator = zeros(size(Matrix_1,1), size(Matrix_1,1));
    denominator = zeros(size(Matrix_1,1), size(Matrix_1,1));

    % 각 column (postsynaptic neuron)에 대해 min/max 계산 후 누적
    for k = 1:size(Matrix_1,2)  % 각 postsynaptic 뉴런마다
        a_col = Matrix_1(:,k);   % A 클러스터: 200x1
        b_col = Matrix_1(:,k)';  % B 클러스터: 1x10

        min_mat = min(a_col, b_col);  % (200x10) – 각 (i,j)에 대해 min(a_ik, b_jk)
        max_mat = max(a_col, b_col);  % (200x10)

        numerator = numerator + min_mat;
        denominator = denominator + max_mat;
    end

    % Weighted Jaccard matrix
    weighted_jaccard_mat = numerator ./ denominator;
    weighted_jaccard_mat(denominator == 0) = NaN;  % 정의 안 되는 경우

    n_wj = size(weighted_jaccard_mat,1);
    weighted_jaccard_mat(1:n_wj+1:end) = NaN;  % diag를 NaN으로

    % 클러스터 간 공유도 = 평균 유사도
    MetaClusterSharednessMatrix(i,i) = mean(weighted_jaccard_mat(:), 'omitnan');
end

%%
figure(10);
set(gcf,'Color','w')
imagesc(MetaClusterSharednessMatrix)
% imagesc(BipartiteConnectivity(DendrogramOrder,DendrogramOrder))
xlim([0.5 7.5])
ylim([0.5 7.5])
% set(gca,'XTick',1:1:size(Post_RightFF_labels,1),'YTick',1:1:size(RightFF_labels,1));
% set(gca,'XTick',1:1:size(RightFF_labels,1),'XTickLabel',RightFF_labels(DendrogramOrder),'YTick',1:1:size(RightFF_labels,1),'YTickLabel',RightFF_labels(DendrogramOrder));
% set(gca,'XTick',1:1:size(Post_RightFF_labels,1),'XTickLabel',Post_RightFF_labels(DendrogramOrder),'YTick',1:1:size(RightFF_labels,1),'YTickLabel',RightFF_labels(DendrogramOrder));

set(gca,'TickDir','out','Box','off')
colors=(colormap('sky'));
colormap(colors)
axis square
grid on
clim([0 0.05]);
title('MetaClusterSharedNess')
% print(gcf,'-depsc2','-vector','figure2_ClusteringResult.eps')
%%
% ClusterSharednessMatrix

% 예시 대칭 행렬
M = ClusterSharednessMatrix * 100;

n = size(M,1);

% 대각 요소
diag_vals = diag(M);
mean_diag = mean(diag_vals,'omitnan');
std_diag = std(diag_vals,'omitnan');

% 비대각 요소 (중복 제거: 상삼각 행렬에서 대각 제외)
off_diag_vals = M(triu(true(n),1));
mean_off_diag = mean(off_diag_vals,'omitnan');
std_off_diag = std(off_diag_vals,'omitnan');

% 데이터 정리
means = [mean_diag, mean_off_diag];
stds = [std_diag, std_off_diag];

% 가로 막대그래프 + errorbar
figure(11); clf; set(gcf, 'Color', 'w');
hold on;
bar_handle = barh(means, 'FaceColor', [0.0000, 0.4470, 0.7410], 'EdgeColor', 'flat');
errorbar(means, 1:2 , stds, 'horizontal', '.k', 'LineWidth', 1.5);

% 원소 값 점 찍기 (가로로 점 위치 조정)
% scatter(diag_vals, ones(size(diag_vals)) * 1 + 0.05, 10, 'black', 'filled', 'DisplayName', 'Diagonal values');
% scatter(off_diag_vals, ones(size(off_diag_vals)) * 2 + 0.05, 10, 'black', 'filled', 'DisplayName', 'Off-diagonal values');

% 라벨 및 설정
ylim([0.5 2.5]);
xlim([0 50]);

yticks([1 2]);
yticklabels({'Within-cluster', 'Between-cluster'});
xlabel('Projection similarity (%)');
title('Cluster-level Projection Similarity');

set(gca, 'FontSize', 12, 'TickDir', 'out', 'Box', 'off','YDir','reverse');
set(gcf, 'Color', 'w');
% print(gcf, '-depsc2', '-vector', 'figure2_ClusterIndex_Horizontal.eps');

%%
% 예시 대칭 행렬
M =MetaClusterSharednessMatrix*100;


n = size(M,1);

% 대각 요소
diag_vals = diag(M);
mean_diag = mean(diag_vals,'omitnan');
std_diag = std(diag_vals,'omitnan');

% 비대각 요소 (중복 제거: 상삼각 행렬에서 대각 제외)
off_diag_vals = M(triu(true(n),1));
mean_off_diag = mean(off_diag_vals,'omitnan');
std_off_diag = std(off_diag_vals,'omitnan');

% 데이터 정리
means = [mean_diag, mean_off_diag];
stds = [std_diag, std_off_diag];

% 가로 막대그래프 + errorbar
figure(12); clf; set(gcf, 'Color', 'w');
hold on;
bar_handle = barh(means, 'FaceColor', [0.0000, 0.4470, 0.7410], 'EdgeColor', 'flat');
errorbar(means, 1 :2, stds, 'horizontal', '.k', 'LineWidth', 1.5);

% 원소 값 점 찍기 (가로로 점 위치 조정)
% scatter(diag_vals, ones(size(diag_vals)) * 1 + 0.05, 10, 'black', 'filled', 'DisplayName', 'Diagonal values');
% scatter(off_diag_vals, ones(size(off_diag_vals)) * 2 + 0.05, 10, 'black', 'filled', 'DisplayName', 'Off-diagonal values');

% 라벨 및 설정
ylim([0.5 2.5]);
xlim([0 50]);

yticks([1 2]);
yticklabels({'Within-cluster', 'Between-cluster'});
xlabel('Projection similarity (%)');
title('Cluster-level Projection Similarity');

set(gca, 'FontSize', 12, 'TickDir', 'out', 'Box', 'off','YDir','reverse');
set(gcf, 'Color', 'w');
print(gcf,'-depsc2','-vector','figure2_MetaClusterIndex.eps')

%% 직관적 클러스터/ 슈퍼클러스터간 뉴런별 synapse 간 평균 거리는??

opt = detectImportOptions('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\fafb_v783_princeton_synapse_table.csv');
opt = setvartype(opt,'pre_root_id_720575940','int64');
opt = setvartype(opt,'post_root_id_720575940','int64');
FAFB_synapse_coordinates = readtable('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\fafb_v783_princeton_synapse_table.csv',opt);
FAFB_synapse_coordinates.Properties.VariableNames(11) = "pre_root_id";
FAFB_synapse_coordinates.Properties.VariableNames(12) = "post_root_id";
FAFB_synapse_coordinates.pre_root_id=FAFB_synapse_coordinates.pre_root_id+int64(720575940000000000);
FAFB_synapse_coordinates.post_root_id=FAFB_synapse_coordinates.post_root_id+int64(720575940000000000);
FAFB_synapse_coordinates.pre_x=FAFB_synapse_coordinates.pre_x*1e-9;
FAFB_synapse_coordinates.pre_y=FAFB_synapse_coordinates.pre_y*1e-9;
FAFB_synapse_coordinates.pre_z=FAFB_synapse_coordinates.pre_z*1e-9;
Central_syn_idx=ismember(FAFB_synapse_coordinates.neuropil,FAFBNeuropil_Central);
%% Cluster
ClusterClosenessMatrix=zeros(size(leafOrder,1));

for i=1:1:size(ClusterClosenessMatrix,1)
    Cluster1_root_id=ClusteringValueCounts{leafOrder(i),5};
    n = numel(Cluster1_root_id);
    Cluster1_outSynapses = cell(n, 1);    
    for k=1:1:size(Cluster1_root_id,1)
        current_root_id=Cluster1_root_id(k);
        out_syn_idx=FAFB_synapse_coordinates.pre_root_id==current_root_id&Central_syn_idx;
        % Cluster1_outSynapses{k,1}= FAFB_synapse_coordinates{out_syn_idx,1:3};    
         Cluster1_outSynapses{k} = table2array(FAFB_synapse_coordinates(out_syn_idx, 1:3));
    end

    for j = i+1:size(ClusterClosenessMatrix,2)

        Cluster2_root_id=ClusteringValueCounts{leafOrder(j),5};
  
        m = numel(Cluster2_root_id);
        Cluster2_outSynapses = cell(m, 1);
        for k=1:1:size(Cluster2_root_id,1)
            current_root_id=Cluster2_root_id(k);
            out_syn_idx=FAFB_synapse_coordinates.pre_root_id==current_root_id&Central_syn_idx;
            Cluster2_outSynapses{k} = table2array(FAFB_synapse_coordinates(out_syn_idx, 1:3));
        end
        MeanPairDist=zeros(n,m);
        for k = 1:n
            A = Cluster1_outSynapses{k};   % K_a x 3

            for l = 1:m
                B = Cluster2_outSynapses{l};  % K_b x 3

                % 모든 쌍 거리
                D = pdist2(A, B, 'euclidean');   % K_a x K_b

                % 1) 모든 쌍의 평균거리
                MeanPairDist(k,l) = mean(D(:),'omitnan');
            end
        end


        % 클러스터 간 공유도 = 평균 유사도
        ClusterClosenessMatrix(i,j) = mean(MeanPairDist(:),'omitnan');
        ClusterClosenessMatrix(j,i) = mean(MeanPairDist(:),'omitnan');

    end
end
%%
for i=1:1:size(ClusterClosenessMatrix,1)
    Cluster1_root_id=ClusteringValueCounts{leafOrder(i),5};
    n = numel(Cluster1_root_id);
    Cluster1_outSynapses = cell(n, 1);    
    for k=1:1:size(Cluster1_root_id,1)
        current_root_id=Cluster1_root_id(k);
        out_syn_idx=FAFB_synapse_coordinates.pre_root_id==current_root_id&Central_syn_idx;
        % Cluster1_outSynapses{k,1}= FAFB_synapse_coordinates{out_syn_idx,1:3};    
        Cluster1_outSynapses{k} = table2array(FAFB_synapse_coordinates(out_syn_idx, 1:3));
    end
        MeanPairDist=zeros(n,n);
        for k = 1:n
            A = Cluster1_outSynapses{k};   % K_a x 3

            for l = 1:n
                B = Cluster1_outSynapses{l};  % K_b x 3

                % 모든 쌍 거리
                D = pdist2(A, B, 'euclidean');   % K_a x K_b

                % 1) 모든 쌍의 평균거리
                MeanPairDist(k,l) = mean(D(:),'omitnan');
            end
        end
        n_md = size(MeanPairDist,1);
        MeanPairDist(1:n_md+1:end) = NaN;  % diag를 NaN으로

    % 클러스터 간 공유도 = 평균 유사도
    ClusterClosenessMatrix(i,i) = mean(MeanPairDist(:),'omitnan');    
end
%% meta cluster 평균 거리.. 다른것

MetaClusterClosenessMatrix=zeros(size(MetaClusterOrder,2));

for i=1:1:size(MetaClusterClosenessMatrix,1)
    Cluster1_root_id=ClusteringValueCounts(MetaClusterOrder{i},5);
    Cluster1_root_id = vertcat(Cluster1_root_id{:});
    n = numel(Cluster1_root_id);
    Cluster1_outSynapses = cell(n, 1);    

    for k=1:1:size(Cluster1_root_id,1)
        current_root_id=Cluster1_root_id(k);
        out_syn_idx=FAFB_synapse_coordinates.pre_root_id==current_root_id&Central_syn_idx;
        % Cluster1_outSynapses{k,1}= FAFB_synapse_coordinates{out_syn_idx,1:3};    
         Cluster1_outSynapses{k} = table2array(FAFB_synapse_coordinates(out_syn_idx, 1:3));
    end

    for j = i+1:size(MetaClusterClosenessMatrix,2)

        Cluster2_root_id=ClusteringValueCounts(MetaClusterOrder{j},5);
        Cluster2_root_id = vertcat(Cluster2_root_id{:});  
        
        m = numel(Cluster2_root_id);
        Cluster2_outSynapses = cell(m, 1);
        for k=1:1:size(Cluster2_root_id,1)
            current_root_id=Cluster2_root_id(k);
            out_syn_idx=FAFB_synapse_coordinates.pre_root_id==current_root_id&Central_syn_idx;
            Cluster2_outSynapses{k} = table2array(FAFB_synapse_coordinates(out_syn_idx, 1:3));
        end
        MeanPairDist=zeros(n,m);
        for k = 1:n
            A = Cluster1_outSynapses{k};   % K_a x 3

            for l = 1:m
                B = Cluster2_outSynapses{l};  % K_b x 3

                % 모든 쌍 거리
                D = pdist2(A, B, 'euclidean');   % K_a x K_b

                % 1) 모든 쌍의 평균거리
                MeanPairDist(k,l) = mean(D(:),'omitnan');
            end
        end


        % 클러스터 간 공유도 = 평균 유사도
        MetaClusterClosenessMatrix(i,j) = mean(MeanPairDist(:),'omitnan');
        MetaClusterClosenessMatrix(j,i) = mean(MeanPairDist(:),'omitnan');

    end
end
%% meta cluster 평균 거리.. 같은 것
for i=1:1:size(MetaClusterClosenessMatrix,1)
    Cluster1_root_id=ClusteringValueCounts(MetaClusterOrder{i},5);
    Cluster1_root_id = vertcat(Cluster1_root_id{:});
    n = numel(Cluster1_root_id);
    Cluster1_outSynapses = cell(n, 1);    

    for k=1:1:size(Cluster1_root_id,1)
        current_root_id=Cluster1_root_id(k);
        out_syn_idx=FAFB_synapse_coordinates.pre_root_id==current_root_id & Central_syn_idx;
        % Cluster1_outSynapses{k,1}= FAFB_synapse_coordinates{out_syn_idx,1:3};    
        Cluster1_outSynapses{k} = table2array(FAFB_synapse_coordinates(out_syn_idx, 1:3));
    end

   
    MeanPairDist=zeros(n,n);
    for k = 1:n
        A = Cluster1_outSynapses{k};   % K_a x 3

        for l = 1:n
            B = Cluster1_outSynapses{l};  % K_b x 3

            % 모든 쌍 거리
            D = pdist2(A, B, 'euclidean');   % K_a x K_b

            % 1) 모든 쌍의 평균거리
            MeanPairDist(k,l) = mean(D(:),'omitnan');
        end
    end

    n_md = size(MeanPairDist,1);
    MeanPairDist(1:n_md+1:end) = NaN;  % diag를 NaN으로
    
    MetaClusterClosenessMatrix(i,i) = mean(MeanPairDist(:),'omitnan');

end
%%
% ClusterSharednessMatrix

% 예시 대칭 행렬
M = ClusterClosenessMatrix ;

n = size(M,1);

% 대각 요소
diag_vals = diag(M);
mean_diag = mean(diag_vals,'omitnan');
std_diag = std(diag_vals,'omitnan');

% 비대각 요소 (중복 제거: 상삼각 행렬에서 대각 제외)
off_diag_vals = M(triu(true(n),1));
mean_off_diag = mean(off_diag_vals,'omitnan');
std_off_diag = std(off_diag_vals,'omitnan');

% 데이터 정리
means = [mean_diag, mean_off_diag];
stds = [std_diag, std_off_diag];

% 가로 막대그래프 + errorbar
figure(13); clf; set(gcf, 'Color', 'w');
hold on;
bar_handle = barh(means, 'FaceColor', [0.0000, 0.4470, 0.7410], 'EdgeColor', 'flat');
errorbar(means, 1:2 , stds, 'horizontal', '.k', 'LineWidth', 1.5);

% 원소 값 점 찍기 (가로로 점 위치 조정)
% scatter(diag_vals, ones(size(diag_vals)) * 1 + 0.05, 10, 'black', 'filled', 'DisplayName', 'Diagonal values');
% scatter(off_diag_vals, ones(size(off_diag_vals)) * 2 + 0.05, 10, 'black', 'filled', 'DisplayName', 'Off-diagonal values');

% 라벨 및 설정
ylim([0.5 2.5]);
% xlim([0 50]);

yticks([1 2]);
yticklabels({'Within-cluster', 'Between-cluster'});
xlabel('distance um');
title('Cluster-level outsynapse Distance');

set(gca, 'FontSize', 12, 'TickDir', 'out', 'Box', 'off','YDir','reverse');
set(gcf, 'Color', 'w');
print(gcf, '-depsc2', '-vector', 'figure2_ClusterDistance_Horizontal.eps');


%%
% ClusterSharednessMatrix

% 예시 대칭 행렬
M = MetaClusterClosenessMatrix ;

n = size(M,1);

% 대각 요소
diag_vals = diag(M);
mean_diag = mean(diag_vals,'omitnan');
std_diag = std(diag_vals,'omitnan');

% 비대각 요소 (중복 제거: 상삼각 행렬에서 대각 제외)
off_diag_vals = M(triu(true(n),1));
mean_off_diag = mean(off_diag_vals,'omitnan');
std_off_diag = std(off_diag_vals,'omitnan');

% 데이터 정리
means = [mean_diag, mean_off_diag];
stds = [std_diag, std_off_diag];

% 가로 막대그래프 + errorbar
figure(14); clf; set(gcf, 'Color', 'w');
hold on;
bar_handle = barh(means, 'FaceColor', [0.0000, 0.4470, 0.7410], 'EdgeColor', 'flat');
errorbar(means, 1:2 , stds, 'horizontal', '.k', 'LineWidth', 1.5);

% 원소 값 점 찍기 (가로로 점 위치 조정)
% scatter(diag_vals, ones(size(diag_vals)) * 1 + 0.05, 10, 'black', 'filled', 'DisplayName', 'Diagonal values');
% scatter(off_diag_vals, ones(size(off_diag_vals)) * 2 + 0.05, 10, 'black', 'filled', 'DisplayName', 'Off-diagonal values');

% 라벨 및 설정
ylim([0.5 2.5]);
xlim([0 2e-4]);

yticks([1 2]);
yticklabels({'Within-cluster', 'Between-cluster'});
xlabel('distance um');
title('Meta-Cluster-level outsynapse Distance');

set(gca, 'FontSize', 12, 'TickDir', 'out', 'Box', 'off','YDir','reverse');
set(gcf, 'Color', 'w');
print(gcf, '-depsc2', '-vector', 'figure2_MetaClusterDistance_Horizontal.eps');
%%
figure(15);
imagesc(ClusterClosenessMatrix)
set(gca,'TickDir','out','Box','off')
colors=(colormap('sky'));
colormap(colors)
axis square
grid on
% clim([0 0.05]);
%%
figure(16);
imagesc(MetaClusterClosenessMatrix)
set(gca,'TickDir','out','Box','off')
colors=(colormap('sky'));
colormap(colors)
axis square
grid on
% clim([0 0.05]);
%% 직관적 하나의 뉴런이 같은 클러스터로 분류되어 있는가?

% ------------------------------
% [1] 데이터 준비
% ------------------------------
types = ClusteringResult.type;         % 예: {'LC9','LC9','LC10',...}
clusters = ClusteringResult.Cluster;   % 예: [1,1,2,...]
n_neurons = length(types);

% ------------------------------
% [2] 타입별 클러스터 일치도 계산
% ------------------------------
unique_types = unique(types);
n_types = length(unique_types);

type_mean_ratios = nan(n_types, 1);  % 타입별 평균 일치도 저장

for t = 1:n_types
    type_name = unique_types{t};
    idx = strcmp(types, type_name);      % 해당 타입의 뉴런 인덱스
    type_cluster = clusters(idx);        % 해당 타입 뉴런들의 클러스터 할당

    % 뉴런 수가 1개면 비교 불가
    if sum(idx) <= 1
        continue;
    end

    % 타입 내 모든 쌍에서 같은 클러스터에 속한 경우 계산
    n = sum(idx);
    count_same_cluster = 0;
    total_pairs = 0;

    for i = 1:n
        for j = i+1:n
            total_pairs = total_pairs + 1;
            if type_cluster(i) == type_cluster(j)
                count_same_cluster = count_same_cluster + 1;
            end
        end
    end

    % 타입별 평균 일치도 저장
    type_mean_ratios(t) = count_same_cluster / total_pairs;
end

% ------------------------------
% [3] 시각화
% ------------------------------
figure; set(gcf, 'Color', 'w');
bar(type_mean_ratios, 'FaceColor', [0.2 0.6 0.8]);
xticks(1:n_types);
xticklabels(unique_types);
xtickangle(45);
ylabel('Intra-type clustering consistency');
title('Neuron type clustering consistency (same cluster rate)');
ylim([0 1]);
box off;
%%
% ------------------------------
% [3-NEW] scatter 시각화: 뉴런 수 vs. 일치도 + 타입 이름 표시
% ------------------------------
type_counts = cellfun(@(x) sum(strcmp(types, x)), unique_types);  % 각 타입의 뉴런 수

figure; set(gcf, 'Color', 'w');
scatter(type_counts, type_mean_ratios, 50, 'filled', 'MarkerFaceColor', [0 0.4470 0.7410]);
xlabel('Number of neurons per type');
ylabel('Intra-type clustering consistency');
title('Consistency vs. Type Size');
grid on;
xlim([0, max(type_counts)*1.05]);
ylim([0, 1]);

% --- 타입 이름 텍스트 추가 ---
% hold on;
% for i = 1:length(unique_types)
%     text(type_counts(i), type_mean_ratios(i), unique_types{i}, ...
%         'FontSize', 8, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
% end

print(gcf,'-depsc2','-vector','TypeConsistency.eps')

%%
% 히스토그램 bin 계산
[bin_counts, bin_edges] = histcounts(type_mean_ratios, 10);
bin_centers = bin_edges(1:end-1) + diff(bin_edges)/2;

% 가로 막대 그래프 그리기
figure; set(gcf, 'Color', 'w');
barh(bin_centers, bin_counts, 'FaceColor', [0.2 0.6 0.8], 'EdgeColor', 'k');
xlabel('Count');
ylabel('Intra-type clustering consistency');
title('Horizontal Histogram of Intra-type Clustering Consistency');
xlim([0, max(bin_counts)*1.1]);  % 여유 공간
ylim([0 1])
grid on;
set(gca, 'TickDir', 'out', 'Box', 'off');
print(gcf,'-depsc2','-vector','TypeConsistencyHisto_Horizontal.eps');

%%

%%
% -----------------------------------
% [1] 가중 평균 및 가중 표준편차 계산
% -----------------------------------
type_counts = cellfun(@(x) sum(strcmp(types, x)), unique_types);  % 각 타입의 뉴런 수
x = type_mean_ratios;
w = type_counts;

% 가중 평균
overall_weighted_mean = sum(w .* x, 'omitnan') / sum(w, 'omitnan');

% 가중 표준편차
numerator = sum(w .* (x - overall_weighted_mean).^2, 'omitnan');
denominator = sum(w, 'omitnan');
weighted_std = sqrt(numerator / denominator);

% -----------------------------------
% [2] bar + errorbar 시각화
% -----------------------------------
figure; set(gcf, 'Color', 'w');
bar(1, overall_weighted_mean, 0.5, 'FaceColor', [0.2 0.6 0.8]); hold on;
errorbar(1, overall_weighted_mean, weighted_std, 'k.', 'LineWidth', 1.5);

xlim([0.5 1.5]);
ylim([0 1]);
xticks(1);
xticklabels({'All types'});
ylabel('Weighted mean intra-type clustering consistency');
title('Overall clustering consistency (weighted mean ± std)');
box off;

%% Merging neurons 중 의미찾기
threshold = 0.05; %% 갯수가 아닌 시냅스 강도로 %퍼센트
WantTarget={'LC10a'
'LC10c';
'LC10d';
'LC10e';
'MeTu3c';
'MeTu3b';
'MeTu3a';
'MeTu4a';
'MeTu4c';
'MeTu4d';
'MeTu4b'
'MeTu2a';
'MeTu2b';};
TargetLength=length(WantTarget);

WantTargetN=[];
for i=1:1:TargetLength
    WantTargetN(i)=sum(ismember(FAFB_consolidated_cell_types.primary_type,WantTarget{i,1}));
    WantTarget_root_ids{i,1}=FAFB_consolidated_cell_types.root_id(ismember(FAFB_consolidated_cell_types.primary_type,WantTarget{i,1}));
end
% % ---------- 2. 전체 post_root_id의 총 입력 시냅스 수 계산 ----------
AllPostIDs = unique(FAFBConnections.post_root_id);
% total_syn_map = containers.Map('KeyType','int64','ValueType','double');
% 
% for i = 1:size(FAFBConnections,1)
%     pid = FAFBConnections.post_root_id(i);
%     syn = FAFBConnections.syn_count(i);
%     if isKey(total_syn_map, pid)
%         total_syn_map(pid) = total_syn_map(pid) + syn;
%     else
%         total_syn_map(pid) = syn;
%     end
% end
load('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\total_syn_map.mat')

%%% % ---------- 3. 각 뉴런 군이 보낸 시냅스 수 저장 ----------
N_group = size(WantTarget_root_ids,1); % 군 개수 (예: 4개)
InputCountByGroup = cell(N_group,1); % 각 군에서 입력받은 post ID 및 시냅스 수

for g = 1:N_group
    current_group = WantTarget_root_ids{g,1};
    group_rows = ismember(FAFBConnections.pre_root_id, current_group);

    FAFBConnections_group = FAFBConnections(group_rows, :); % 이 군에서 보낸 연결만 추출

    % post_id 별로 시냅스 수 집계
    post_ids = unique(FAFBConnections_group.post_root_id);
    syn_count_map = containers.Map('KeyType','int64','ValueType','double');

    for i = 1:height(FAFBConnections_group)
        pid = FAFBConnections_group.post_root_id(i);
        syn = FAFBConnections_group.syn_count(i);
        if isKey(syn_count_map, pid)
            syn_count_map(pid) = syn_count_map(pid) + syn;
        else
            syn_count_map(pid) = syn;
        end
    end
    InputCountByGroup{g} = syn_count_map;
end
%%% ---------- 4. 각 post_root_id에 대해 5% 이상 기여한 군 수 세기 ----------
qualified_post_ids = int64([]);

for i = 1:length(AllPostIDs)
    pid = AllPostIDs(i);
    if ~isKey(total_syn_map, pid)
        continue;
    end
    total_input = total_syn_map(pid);
    count_above_thresh = 0;

    for g = 1:N_group
        syn_map = InputCountByGroup{g};
        if isKey(syn_map, pid)
            prop = syn_map(pid) / total_input;
            if prop >= threshold
                count_above_thresh = count_above_thresh + 1;
            end
        end
    end

    if count_above_thresh >= 2
        qualified_post_ids(end+1,1) = pid;
    end
end

%%% 군 개수 자동 설정
N_group = size(WantTarget_root_ids, 1);

% post_id와 각 군에서 받은 비율 저장
n_qualified = length(qualified_post_ids);
group_input_percent = zeros(n_qualified, N_group);  % [n_post_ids x n_groups]

for i = 1:n_qualified
    pid = qualified_post_ids(i);
    total_input = total_syn_map(pid);

    for g = 1:N_group
        syn_map = InputCountByGroup{g};
        if isKey(syn_map, pid)
            group_input_percent(i, g) = syn_map(pid) / total_input;
        else
            group_input_percent(i, g) = 0;
        end
    end
end

% 동적으로 변수명 생성 (첫 번째는 post_root_id)
varNames = {'post_root_id'};
for g = 1:N_group
    % WantTarget_root_ids{i,2}에 저장된 군 이름 사용
    group_name = WantTarget{g,1};
    % 공백 등 제거하고 변수명으로 적절히 가공
    group_varname = matlab.lang.makeValidName(sprintf('%s_percent', group_name));
    varNames{end+1} = group_varname;
end

% 테이블 생성
% 각각 따로 테이블 생성
T_post = table(qualified_post_ids, 'VariableNames', {'post_root_id'});
T_group = array2table(group_input_percent, 'VariableNames', varNames(2:end));  % group 컬럼만

% 수평 결합
Merging_result = [T_post, T_group];
for i=1:1:size(Merging_result,1)
    idx=FAFB_consolidated_cell_types.root_id==Merging_result.post_root_id(i);
    if ~any(idx)
       Merging_result.type{i}='';
       continue;
    end
    Merging_result.type(i)=FAFB_consolidated_cell_types.primary_type(idx);
    
end
Merging_result = sortrows(Merging_result,"type","descend");

%%
sortingOrder=[21 24 28 15 17 25 29 9 18 19 26 27 30 3 5 6 10 23 4 7 12 13 14 22 1 2 8 20 11 16 31 32]-1;
WantMatrix_Out_Sorted=[];

for i=1:1:length(sortingOrder)
    RightFFIdx=ClusteringResult.Cluster==sortingOrder(i);
    WantMatrix_Out_Sorted_Temp=[];
    for j=1:1:length(sortingOrder)
        PostIdx=ClusteringResult_PostNeurons.Cluster==sortingOrder(j);
        WantMatrix_Out_Sorted_Temp=[WantMatrix_Out_Sorted_Temp WantMatrix_Out(RightFFIdx,PostIdx)];
    end
    WantMatrix_Out_Sorted=[WantMatrix_Out_Sorted;WantMatrix_Out_Sorted_Temp];
end

%%
figure;set(gcf,'Color','w')
imagesc(WantMatrix_Out_Sorted>0)
colors=(colormap('gray'));
colormap((colors))
YTick=cumsum(cell2mat(ClusteringValueCounts(sortingOrder+1,2)));
YTick(end)=[];
XTick=cumsum(cell2mat(ClusteringValueCounts_PostNeurons(sortingOrder+1,2)));
XTick(end)=[];
axis image
set(gca,'Box','off','XTick',XTick','YTick',YTick,'XTickLabel',[],'YTickLabel',[],'TickDir','none')
grid on
% 축 속성을 가져옵니다.
ax = gca;

% 그리드 색상을 변경합니다.
ax.GridColor = [0.8 0.8 0.8]; % 빨간색
ax.GridAlpha = 0.2; % 그리드 투명도 설정 (0: 투명, 1: 불투명)

% 원하는 경우 그리드 스타일도 변경할 수 있습니다.
ax.GridLineStyle = '--'; % 점선 스타일

function leafOrder = get_leaf_order_from_linkage(Z)
    % Z: linkage matrix (size: [n-1, 3])
    % return: leafOrder - reordered list of leaf indices

    n = size(Z, 1) + 1;
    total_nodes = 2 * n - 1;

    % 각 클러스터 인덱스마다 포함된 리프 노드 기록
    cluster_leaves = cell(total_nodes, 1);
    
    % 초기 리프 노드: 자기 자신만 포함
    for i = 1:n
        cluster_leaves{i} = i;
    end

    % linkage 행을 따라 병합하며 leaf 집합 생성
    for i = 1:n-1
        c1 = Z(i,1);
        c2 = Z(i,2);
        new_cluster = n + i;
        cluster_leaves{new_cluster} = [cluster_leaves{c1}; cluster_leaves{c2}];
    end

    % 최종 트리 루트 노드에서 리프 순서 반환
    leafOrder = cluster_leaves{2*n - 1};
end
