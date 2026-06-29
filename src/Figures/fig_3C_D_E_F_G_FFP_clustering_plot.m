%% Load data
clear all; close all; clc;

baseDir = fileparts(fileparts(mfilename('fullpath')));  % resolved from this script's location (portable; no absolute paths)

% Leiden cluster assignment of the FFP neurons (rows), produced by
% Data_Processing/s09_FFP_leiden_clustering.ipynb
opt = detectImportOptions(fullfile(baseDir,'Processed_Data','leiden_right_FFP_output_CB_no_VPN_thr0_resol3_100000.csv'));
opt = setvartype(opt,'root_id','int64');
ClusteringResult=readtable(fullfile(baseDir,'Processed_Data','leiden_right_FFP_output_CB_no_VPN_thr0_resol3_100000.csv'),opt);

% Codex connectivity (no synapse threshold)
opt = detectImportOptions(fullfile(baseDir,'Codex_Data','connections_no_threshold.csv'));
opt = setvartype(opt,'pre_root_id','int64');
opt = setvartype(opt,'post_root_id','int64');
FAFBConnections = readtable(fullfile(baseDir,'Codex_Data','connections_no_threshold.csv'),opt);

FAFBNeuropils=unique(FAFBConnections.neuropil);
FAFBNeuropil_OpticLobeRight={'AME_R','ME_R','LO_R','LOP_R','LA_R'};
FAFBNeuropil_OpticLobeLeft={'AME_L','ME_L','LO_L','LOP_L','LA_L'};
FAFBNeuropil_Central=FAFBNeuropils(~ismember(FAFBNeuropils,FAFBNeuropil_OpticLobeRight));
FAFBNeuropil_Central=FAFBNeuropil_Central(~ismember(FAFBNeuropil_Central,FAFBNeuropil_OpticLobeLeft));
FAFBNeuropil_Central=FAFBNeuropil_Central(~ismember(FAFBNeuropil_Central,'UNASGD'));

opt = detectImportOptions(fullfile(baseDir,'Codex_Data','consolidated_cell_types.csv'));
opt = setvartype(opt,'root_id','int64');
FAFB_consolidated_cell_types = readtable(fullfile(baseDir,'Codex_Data','consolidated_cell_types.csv'),opt);

% Inter-cluster (bipartite) connectivity matrix, produced by s09
opt = detectImportOptions(fullfile(baseDir,'Processed_Data','leiden_right_FFP_intercluster_connectivity_CB_no_VPN_thr0_resol3_100000.csv'));
BipartiteConnectivity = readmatrix(fullfile(baseDir,'Processed_Data','leiden_right_FFP_intercluster_connectivity_CB_no_VPN_thr0_resol3_100000.csv'),opt);
BipartiteConnectivity = reshape(BipartiteConnectivity,[max(ClusteringResult.Cluster)+1,max(ClusteringResult.Cluster)+1]);

% Raw-similarity dendrogram (linkage matrix), produced by s09
opt = detectImportOptions(fullfile(baseDir,'Processed_Data','dendrogram_raw_CB_no_VPN_thr0_resol3_100000.csv'));
RightFFP_Dendrogram = readmatrix(fullfile(baseDir,'Processed_Data','dendrogram_raw_CB_no_VPN_thr0_resol3_100000.csv'),opt);
RightFFP_Dendrogram(:,4)=[];
RightFFP_Dendrogram(:,1)=RightFFP_Dendrogram(:,1)+1;
RightFFP_Dendrogram(:,2)=RightFFP_Dendrogram(:,2)+1;

% Post-neuron root_ids (provides Post_Want), produced by s08
load(fullfile(baseDir,'Processed_Data','post_neurons_FFP_opticlobe_central_no_VPN_thr0.mat'))
% Leiden cluster assignment of the post neurons (columns), produced by s09
ClusteringResult_PostNeurons=readtable(fullfile(baseDir,'Processed_Data','leiden_post_right_FFP_CB_no_VPN_thr0_100000.csv'));
ClusteringResult_PostNeurons(1,:) = [];
ClusteringResult_PostNeurons = removevars(ClusteringResult_PostNeurons, "Var1");
ClusteringResult_PostNeurons.Properties.VariableNames(1) = "root_id";
ClusteringResult_PostNeurons.Properties.VariableNames(2) = "Cluster";
ClusteringResult_PostNeurons.root_id=Post_Want;
% FFP -> central-brain output adjacency matrix (provides WantMatrix_Out), produced by s08
load(fullfile(baseDir,'Processed_Data','right_FFP_output_matrix_CB_no_VPN_thr0.mat'))

%% Summarize FFP clusters (rows)
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

%% Summarize post-neuron clusters (columns)
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

%% Manual labelling (down to 5% of cluster size)
for i=1:1:size(ClusteringValueCounts,1)
    ClusteringValueCounts{i,7}=size(ClusteringValueCounts{i,5},1)*0.05;
end
for i=1:1:size(ClusteringValueCounts_PostNeurons,1)
    ClusteringValueCounts_PostNeurons{i,7}=size(ClusteringValueCounts_PostNeurons{i,5},1)*0.05;
end

RightFFP_labels={};
Post_RightFFP_labels={};
for i=1:1:size(ClusteringValueCounts,1)
    RightFFP_labels_temp=[];
    for j=1:1:size(ClusteringValueCounts{i, 6},1)
        if ClusteringValueCounts{i, 6}{j, 2}>ClusteringValueCounts{i,7}
            RightFFP_labels_temp=[RightFFP_labels_temp ', ' ClusteringValueCounts{i, 6}{j, 1}];
        end
    end
    RightFFP_labels_temp(1:2)=[];
    RightFFP_labels{i,1}=RightFFP_labels_temp;

    if i<=size(ClusteringValueCounts_PostNeurons,1)
        Post_RightFFP_labels{i,1}=ClusteringValueCounts_PostNeurons{i,4};
    else
        Post_RightFFP_labels{i,1}='';
    end
end

%% figure 1: dendrogram reordered by leaf order
% RightFFP_Dendrogram is a linkage matrix
leafOrder = get_leaf_order_from_linkage(RightFFP_Dendrogram);
temp_RightFFP_labels=RightFFP_labels;

figure(1);
set(gcf,'Color','w')
[H,~,order] = dendrogram(RightFFP_Dendrogram, 0, ...
    'Reorder', leafOrder, 'ColorThreshold', 1, 'Labels', temp_RightFFP_labels);
set(H,'LineWidth',2);
set(gca,'TickDir','out');

%% figure 2: inter-cluster connectivity reordered by leaf order
figure(2);
set(gcf,'Color','w')
imagesc(log10(1+BipartiteConnectivity(leafOrder,leafOrder)))
set(gca,'XTick',1:1:size(Post_RightFFP_labels,1),'XTickLabel',Post_RightFFP_labels(leafOrder),'YTick',1:1:size(RightFFP_labels,1),'YTickLabel',RightFFP_labels(leafOrder));
set(gca,'TickDir','out','Box','off')
colors=(colormap('sky'));
colormap(colors)
axis square
grid on

%% Cluster projection similarity (weighted Jaccard), in dendrogram (leaf) order
ClusterSharednessMatrix=zeros(size(leafOrder,1));

for i=1:1:size(ClusterSharednessMatrix,1)
    Cluster1_root_id=ClusteringValueCounts{leafOrder(i),5};
    idx_1=ismember(ClusteringResult.root_id,Cluster1_root_id);
    Matrix_1=WantMatrix_Out(idx_1,:);

    for j = i+1:size(ClusterSharednessMatrix,2)

        Cluster2_root_id=ClusteringValueCounts{leafOrder(j),5};
        idx_2=ismember(ClusteringResult.root_id,Cluster2_root_id);
        Matrix_2=WantMatrix_Out(idx_2,:);

        % Initialize
        numerator = zeros(size(Matrix_1,1), size(Matrix_2,1));
        denominator = zeros(size(Matrix_1,1), size(Matrix_2,1));

        % Accumulate min/max over each column (postsynaptic neuron)
        for k = 1:size(Matrix_1,2)  % for each postsynaptic neuron
            a_col = Matrix_1(:,k);   % cluster A: 200x1
            b_col = Matrix_2(:,k)';  % cluster B: 1x10

            min_mat = min(a_col, b_col);  % (200x10) - min(a_ik, b_jk) for each (i,j)
            max_mat = max(a_col, b_col);  % (200x10)

            numerator = numerator + min_mat;
            denominator = denominator + max_mat;
        end

        % Weighted Jaccard matrix
        weighted_jaccard_mat = numerator ./ denominator;
        weighted_jaccard_mat(denominator == 0) = NaN;  % undefined case

        % cluster-to-cluster sharedness = mean similarity
        ClusterSharednessMatrix(i,j) = mean(weighted_jaccard_mat(:), 'omitnan');
        ClusterSharednessMatrix(j,i) = mean(weighted_jaccard_mat(:), 'omitnan');

    end
end

for i=1:1:size(ClusterSharednessMatrix,1)
    Cluster1_root_id=ClusteringValueCounts{leafOrder(i),5};
    idx_1=ismember(ClusteringResult.root_id,Cluster1_root_id);
    Matrix_1=WantMatrix_Out(idx_1,:);

        % Initialize
        numerator = zeros(size(Matrix_1,1), size(Matrix_1,1));
        denominator = zeros(size(Matrix_1,1), size(Matrix_1,1));

        % Accumulate min/max over each column (postsynaptic neuron)
        for k = 1:size(Matrix_1,2)  % for each postsynaptic neuron
            a_col = Matrix_1(:,k);   % cluster A: 200x1
            b_col = Matrix_1(:,k)';  % cluster B: 1x10

            min_mat = min(a_col, b_col);  % (200x10) - min(a_ik, b_jk) for each (i,j)
            max_mat = max(a_col, b_col);  % (200x10)

            numerator = numerator + min_mat;
            denominator = denominator + max_mat;
        end
        % Weighted Jaccard matrix
        weighted_jaccard_mat = numerator ./ denominator;
        weighted_jaccard_mat(denominator == 0) = NaN;  % undefined case

        n_wj = size(weighted_jaccard_mat,1);
        weighted_jaccard_mat(1:n_wj+1:end) = NaN;  % set diagonal to NaN

        % cluster-to-cluster sharedness = mean similarity
        ClusterSharednessMatrix(i,i) = mean(weighted_jaccard_mat(:), 'omitnan');
end

%% figure 3: cluster-level projection similarity boxplot
% ClusterSharednessMatrix
M = ClusterSharednessMatrix * 100;
n = size(M,1);

% Diagonal entries (within-cluster)
diag_vals = diag(M);
diag_vals = diag_vals(~isnan(diag_vals)); % remove NaN for omitnan handling

% Off-diagonal entries (between-cluster, deduplicated)
off_diag_vals = M(triu(true(n),1));
off_diag_vals = off_diag_vals(~isnan(off_diag_vals)); % remove NaN for omitnan handling

% --- Build boxplot data and groups ---
data_to_plot = [diag_vals; off_diag_vals];
group_indices = [ones(size(diag_vals)); ...
                 2 * ones(size(off_diag_vals))];
group_labels = {'Within-cluster', 'Between-cluster'};

% --- Boxplot ---
figure(3); clf; set(gcf, 'Color', 'w');
hold on;

% horizontal boxplot
boxplot(data_to_plot, group_indices, ...
    'Orientation', 'horizontal', ... % horizontal orientation
    'Labels', group_labels, ...      % Y-axis labels
    'Notch', 'off', ...
    'Symbol','');                 % Notch off

% match colors to the original (optional)
% box color (blue)
h_boxes = findobj(gca, 'Tag', 'Box');
patch_color = [0.0000, 0.4470, 0.7410];
for j = 1:length(h_boxes)
    patch(get(h_boxes(j), 'XData'), get(h_boxes(j), 'YData'), patch_color, 'FaceAlpha', 0.5);
end
% line color and width (black, 1.5pt)
set(findobj(gca, 'Type', 'Line'), 'Color', 'k', 'LineWidth', 1.5);

% labels and settings
ylim([0.5 2.5]); % Y-axis range
xlim([0 65]);     % X-axis range
xlabel('Projection similarity (%)');
title('Cluster-level Projection Similarity');
set(gca, 'FontSize', 12, 'TickDir', 'out', 'Box', 'off');

% flip Y-axis so within-cluster is on top
set(gca, 'YDir','reverse');

set(gcf, 'Color', 'w');

% --- Rank-sum Test (Wilcoxon rank-sum test) ---
% d1: Within-cluster similarity
% d2: Between-cluster similarity
[p_val, h_stat, stats] = ranksum(diag_vals, off_diag_vals);

% print results
fprintf('\n--- Statistical Test Results ---\n');
fprintf('Method: Wilcoxon rank-sum test (Mann-Whitney U test)\n');
fprintf('P-value: %.4e\n', p_val);

if p_val < 0.05
    fprintf('Result: Statistically Significant (p < 0.05)\n');
else
    fprintf('Result: Not Statistically Significant\n');
end

% --- annotate the P-value on the plot (optional) ---
% place text at the top of the plot
text_str = sprintf('p = %.2e', p_val);
text(max(xlim)*0.7, 1.5, text_str, 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

%% SuperCluster projection similarity (weighted Jaccard)
SuperClusterOrder{1}=[21;24;28];
SuperClusterOrder{2}=[15;17;25;29];
SuperClusterOrder{3}=[9;18;19;26;27;30];
SuperClusterOrder{4}=[3;5;6;10;23];
SuperClusterOrder{5}=[4;7;12;13;14;22];
SuperClusterOrder{6}=[1;2;8;20];
SuperClusterOrder{7}=[11;16;31;32];

SuperClusterSharednessMatrix=zeros(size(SuperClusterOrder,2));

for i=1:1:size(SuperClusterSharednessMatrix,1)
    Cluster1_root_id=ClusteringValueCounts(SuperClusterOrder{i},5);
    Cluster1_root_id = vertcat(Cluster1_root_id{:});

    idx_1=ismember(ClusteringResult.root_id,Cluster1_root_id);
    Matrix_1=WantMatrix_Out(idx_1,:);

    for j=i+1:1:size(SuperClusterSharednessMatrix,2)
        Cluster2_root_id=ClusteringValueCounts(SuperClusterOrder{j},5);
        Cluster2_root_id = vertcat(Cluster2_root_id{:});

        idx_2=ismember(ClusteringResult.root_id,Cluster2_root_id);
        Matrix_2=WantMatrix_Out(idx_2,:);

        % Initialize
        numerator = zeros(size(Matrix_1,1), size(Matrix_2,1));
        denominator = zeros(size(Matrix_1,1), size(Matrix_2,1));

        % Accumulate min/max over each column (postsynaptic neuron)
        for k = 1:size(Matrix_1,2)  % for each postsynaptic neuron
            a_col = Matrix_1(:,k);   % cluster A: 200x1
            b_col = Matrix_2(:,k)';  % cluster B: 1x10

            min_mat = min(a_col, b_col);  % (200x10) - min(a_ik, b_jk) for each (i,j)
            max_mat = max(a_col, b_col);  % (200x10)

            numerator = numerator + min_mat;
            denominator = denominator + max_mat;
        end

        % Weighted Jaccard matrix
        weighted_jaccard_mat = numerator ./ denominator;
        weighted_jaccard_mat(denominator == 0) = NaN;  % undefined case

        % cluster-to-cluster sharedness = mean similarity
        SuperClusterSharednessMatrix(i,j) = mean(weighted_jaccard_mat(:), 'omitnan');
        SuperClusterSharednessMatrix(j,i) = mean(weighted_jaccard_mat(:), 'omitnan');

    end
end

for i=1:1:size(SuperClusterSharednessMatrix,1)
    Cluster1_root_id=ClusteringValueCounts(SuperClusterOrder{i},5);
    Cluster1_root_id = vertcat(Cluster1_root_id{:});

    idx_1=ismember(ClusteringResult.root_id,Cluster1_root_id);
    Matrix_1=WantMatrix_Out(idx_1,:);

    % Initialize
    numerator = zeros(size(Matrix_1,1), size(Matrix_1,1));
    denominator = zeros(size(Matrix_1,1), size(Matrix_1,1));

    % Accumulate min/max over each column (postsynaptic neuron)
    for k = 1:size(Matrix_1,2)  % for each postsynaptic neuron
        a_col = Matrix_1(:,k);   % cluster A: 200x1
        b_col = Matrix_1(:,k)';  % cluster B: 1x10

        min_mat = min(a_col, b_col);  % (200x10) - min(a_ik, b_jk) for each (i,j)
        max_mat = max(a_col, b_col);  % (200x10)

        numerator = numerator + min_mat;
        denominator = denominator + max_mat;
    end

    % Weighted Jaccard matrix
    weighted_jaccard_mat = numerator ./ denominator;
    weighted_jaccard_mat(denominator == 0) = NaN;  % undefined case

    n_wj = size(weighted_jaccard_mat,1);
    weighted_jaccard_mat(1:n_wj+1:end) = NaN;  % set diagonal to NaN

    % cluster-to-cluster sharedness = mean similarity
    SuperClusterSharednessMatrix(i,i) = mean(weighted_jaccard_mat(:), 'omitnan');
end

%% figure 4: supercluster-level projection similarity boxplot
% SuperClusterSharednessMatrix
M = SuperClusterSharednessMatrix * 100;
n = size(M,1);

% Diagonal entries (within-cluster)
diag_vals = diag(M);
diag_vals = diag_vals(~isnan(diag_vals)); % remove NaN for omitnan handling

% Off-diagonal entries (between-cluster, deduplicated)
off_diag_vals = M(triu(true(n),1));
off_diag_vals = off_diag_vals(~isnan(off_diag_vals)); % remove NaN for omitnan handling

% --- Build boxplot data and groups ---
data_to_plot = [diag_vals; off_diag_vals];
group_indices = [ones(size(diag_vals)); ...
                 2 * ones(size(off_diag_vals))];
group_labels = {'Within-cluster', 'Between-cluster'};

% --- Boxplot ---
figure(4); clf; set(gcf, 'Color', 'w');
hold on;

% horizontal boxplot
boxplot(data_to_plot, group_indices, ...
    'Orientation', 'horizontal', ... % horizontal orientation
    'Labels', group_labels, ...      % Y-axis labels
    'Notch', 'off', ...
    'Symbol','');                     % Notch off

% match colors to the original (optional)
% box color (blue)
h_boxes = findobj(gca, 'Tag', 'Box');
patch_color = [0.0000, 0.4470, 0.7410];
for j = 1:length(h_boxes)
    patch(get(h_boxes(j), 'XData'), get(h_boxes(j), 'YData'), patch_color, 'FaceAlpha', 0.5);
end
% line color and width (black, 1.5pt)
set(findobj(gca, 'Type', 'Line'), 'Color', 'k', 'LineWidth', 1.5);

% labels and settings
ylim([0.5 2.5]); % Y-axis range
xlim([0 65]);     % X-axis range
xlabel('Projection similarity (%)');
title('SuperCluster-level Projection Similarity');
set(gca, 'FontSize', 12, 'TickDir', 'out', 'Box', 'off');

% flip Y-axis so within-cluster is on top
set(gca, 'YDir','reverse');

set(gcf, 'Color', 'w');

% --- Rank-sum Test (Wilcoxon rank-sum test) ---
% d1: Within-cluster similarity
% d2: Between-cluster similarity
[p_val, h_stat, stats] = ranksum(diag_vals, off_diag_vals);

% print results
fprintf('\n--- Statistical Test Results ---\n');
fprintf('Method: Wilcoxon rank-sum test (Mann-Whitney U test)\n');
fprintf('P-value: %.4e\n', p_val);

if p_val < 0.05
    fprintf('Result: Statistically Significant (p < 0.05)\n');
else
    fprintf('Result: Not Statistically Significant\n');
end

% --- annotate the P-value on the plot (optional) ---
% place text at the top of the plot
text_str = sprintf('p = %.2e', p_val);
text(max(xlim)*0.7, 1.5, text_str, 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

%% Output neuropil targets per supercluster
% For each supercluster, tabulate the central-brain output synapses by neuropil
% (optic-lobe neuropils and VPN targets are excluded, matching the matrix construction).
opt = detectImportOptions(fullfile(baseDir,'Codex_Data','classification.csv'));
opt = setvartype(opt,'root_id','int64');
FAFBClassification = readtable(fullfile(baseDir,'Codex_Data','classification.csv'),opt);
VPN_root_id=FAFBClassification.root_id(strcmp(FAFBClassification.super_class,'visual_projection'));

SuperCluster_neuropil_targets = cell(size(SuperClusterOrder,2),1);
% Per-supercluster post-synaptic targets, read by Figures/fig_S1E_J_superclusters.m
SuperCluster_targets = struct('n_post_neurons',{}, 'post_type_counts',{}, ...
    'post_neuropil_syn',{}, 'post_type_syn',{});
for i=1:1:size(SuperClusterOrder,2)
    current_root_ids=ClusteringValueCounts(SuperClusterOrder{i},5);
    current_root_ids = vertcat(current_root_ids{:});

    post_idx=ismember(FAFBConnections.pre_root_id,current_root_ids);
    post_Connections=FAFBConnections(post_idx,:);
    % Same filtering as the matrix construction
    OpticR=ismember(post_Connections.neuropil,{'LA_R', 'ME_R','AME_R', 'LO_R','LOP_R','LA_L', 'ME_L','AME_L', 'LO_L','LOP_L'});
    post_Connections(OpticR,:)=[];
    post_Connections(ismember(post_Connections.post_root_id,VPN_root_id),:)=[];

    % Output synapses grouped by neuropil, sorted by total synapse count
    [Post_neuropils,~,ic]=unique(post_Connections.neuropil);
    for j=1:1:size(Post_neuropils,1)
        idx=ic==j;
        Post_neuropils{j,2}=sum(post_Connections.syn_count(idx));
    end
    Post_neuropils=sortrows(Post_neuropils,2,'descend');
    SuperCluster_neuropil_targets{i}=Post_neuropils;

    % Post-synaptic neurons: per-neuron synapse count, mapped to cell type
    [postNeurons,~,ic_n] = unique(post_Connections.post_root_id);
    synPerPost = accumarray(ic_n, post_Connections.syn_count);
    postTypes = strings(numel(postNeurons),1);
    [tf,loc] = ismember(postNeurons, FAFB_consolidated_cell_types.root_id);
    postTypes(tf)  = string(FAFB_consolidated_cell_types.primary_type(loc(tf)));
    postTypes(~tf) = "Untyped";

    [uTypes,~,ic_t] = unique(postTypes);
    SuperCluster_targets(i).n_post_neurons    = numel(postNeurons);
    SuperCluster_targets(i).post_type_counts  = sortrows([cellstr(uTypes), num2cell(accumarray(ic_t,1))],          2, 'descend'); % {type, neuron count}
    SuperCluster_targets(i).post_neuropil_syn = Post_neuropils;                                                                    % {neuropil, synapse count}
    SuperCluster_targets(i).post_type_syn     = sortrows([cellstr(uTypes), num2cell(accumarray(ic_t,synPerPost))], 2, 'descend');  % {type, synapse count}
end

save(fullfile(baseDir,'Processed_Data','FFP_supercluster_targets.mat'), 'SuperCluster_targets');

%% Mean pairwise distance between output synapses, within/between clusters and superclusters

opt = detectImportOptions(fullfile(baseDir,'Codex_Data','fafb_v783_princeton_synapse_table.csv'));
opt = setvartype(opt,'pre_root_id_720575940','int64');
opt = setvartype(opt,'post_root_id_720575940','int64');
FAFB_synapse_coordinates = readtable(fullfile(baseDir,'Codex_Data','fafb_v783_princeton_synapse_table.csv'),opt);
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

                % all pairwise distances
                D = pdist2(A, B, 'euclidean');   % K_a x K_b

                % 1) mean of all pairwise distances
                MeanPairDist(k,l) = mean(D(:),'omitnan');
            end
        end


        % cluster-to-cluster sharedness = mean similarity
        ClusterClosenessMatrix(i,j) = mean(MeanPairDist(:),'omitnan');
        ClusterClosenessMatrix(j,i) = mean(MeanPairDist(:),'omitnan');

    end
end

for i=1:1:size(ClusterClosenessMatrix,1)
    Cluster1_root_id=ClusteringValueCounts{leafOrder(i),5};
    n = numel(Cluster1_root_id);
    Cluster1_outSynapses = cell(n, 1);
    for k=1:1:size(Cluster1_root_id,1)
        current_root_id=Cluster1_root_id(k);
        out_syn_idx=FAFB_synapse_coordinates.pre_root_id==current_root_id&Central_syn_idx;
        Cluster1_outSynapses{k} = table2array(FAFB_synapse_coordinates(out_syn_idx, 1:3));
    end
        MeanPairDist=zeros(n,n);
        for k = 1:n
            A = Cluster1_outSynapses{k};   % K_a x 3

            for l = 1:n
                B = Cluster1_outSynapses{l};  % K_b x 3

                % all pairwise distances
                D = pdist2(A, B, 'euclidean');   % K_a x K_b

                % 1) mean of all pairwise distances
                MeanPairDist(k,l) = mean(D(:),'omitnan');
            end
        end
        n_md = size(MeanPairDist,1);
        MeanPairDist(1:n_md+1:end) = NaN;  % set diagonal to NaN

    % cluster-to-cluster sharedness = mean similarity
    ClusterClosenessMatrix(i,i) = mean(MeanPairDist(:),'omitnan');
end

%% figure 5: cluster-level outsynapse distance boxplot
% ClusterClosenessMatrix
M = ClusterClosenessMatrix;
n = size(M,1);

% Diagonal entries (within-cluster)
diag_vals = diag(M);
diag_vals = diag_vals(~isnan(diag_vals)); % remove NaN values

% Off-diagonal entries (between-cluster, deduplicated)
off_diag_vals = M(triu(true(n),1));
off_diag_vals = off_diag_vals(~isnan(off_diag_vals)); % remove NaN values

% --- Build boxplot data and groups ---
data_to_plot = [diag_vals; off_diag_vals];
group_indices = [ones(size(diag_vals)); ...
                 2 * ones(size(off_diag_vals))];
group_labels = {'Within-cluster', 'Between-cluster'};

% --- Boxplot ---
figure(5); clf; set(gcf, 'Color', 'w');
hold on;

% horizontal boxplot
boxplot(data_to_plot, group_indices, ...
    'Orientation', 'horizontal', ...
    'Labels', group_labels, ...
    'Notch', 'off', ...
    'Symbol','');                     % Notch off


% labels and settings
ylim([0.5 2.5]);
xlim([0 2.5e-4]);
xlabel('distance um');
title('Cluster-level outsynapse Distance');
set(gca, 'FontSize', 12, 'TickDir', 'out', 'Box', 'off');

% flip Y-axis so within-cluster is on top
set(gca, 'YDir','reverse');

set(gcf, 'Color', 'w');
% --- Rank-sum Test (Wilcoxon rank-sum test) ---
% d1: Within-cluster similarity
% d2: Between-cluster similarity
[p_val, h_stat, stats] = ranksum(diag_vals, off_diag_vals);

% print results
fprintf('\n--- Statistical Test Results ---\n');
fprintf('Method: Wilcoxon rank-sum test (Mann-Whitney U test)\n');
fprintf('P-value: %.4e\n', p_val);

if p_val < 0.05
    fprintf('Result: Statistically Significant (p < 0.05)\n');
else
    fprintf('Result: Not Statistically Significant\n');
end

% --- annotate the P-value on the plot (optional) ---
% place text at the top of the plot
text_str = sprintf('p = %.2e', p_val);
text(max(xlim)*0.7, 1.5, text_str, 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

%% SuperCluster mean distance (between clusters)
SuperClusterClosenessMatrix=zeros(size(SuperClusterOrder,2));

for i=1:1:size(SuperClusterClosenessMatrix,1)
    Cluster1_root_id=ClusteringValueCounts(SuperClusterOrder{i},5);
    Cluster1_root_id = vertcat(Cluster1_root_id{:});
    n = numel(Cluster1_root_id);
    Cluster1_outSynapses = cell(n, 1);

    for k=1:1:size(Cluster1_root_id,1)
        current_root_id=Cluster1_root_id(k);
        out_syn_idx=FAFB_synapse_coordinates.pre_root_id==current_root_id&Central_syn_idx;
         Cluster1_outSynapses{k} = table2array(FAFB_synapse_coordinates(out_syn_idx, 1:3));
    end

    for j = i+1:size(SuperClusterClosenessMatrix,2)

        Cluster2_root_id=ClusteringValueCounts(SuperClusterOrder{j},5);
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

                % all pairwise distances
                D = pdist2(A, B, 'euclidean');   % K_a x K_b

                % 1) mean of all pairwise distances
                MeanPairDist(k,l) = mean(D(:),'omitnan');
            end
        end


        % cluster-to-cluster sharedness = mean similarity
        SuperClusterClosenessMatrix(i,j) = mean(MeanPairDist(:),'omitnan');
        SuperClusterClosenessMatrix(j,i) = mean(MeanPairDist(:),'omitnan');

    end
end
%% SuperCluster mean distance (within cluster)
for i=1:1:size(SuperClusterClosenessMatrix,1)
    Cluster1_root_id=ClusteringValueCounts(SuperClusterOrder{i},5);
    Cluster1_root_id = vertcat(Cluster1_root_id{:});
    n = numel(Cluster1_root_id);
    Cluster1_outSynapses = cell(n, 1);

    for k=1:1:size(Cluster1_root_id,1)
        current_root_id=Cluster1_root_id(k);
        out_syn_idx=FAFB_synapse_coordinates.pre_root_id==current_root_id&Central_syn_idx;
         Cluster1_outSynapses{k} = table2array(FAFB_synapse_coordinates(out_syn_idx, 1:3));
    end


    MeanPairDist=zeros(n,n);
    for k = 1:n
        A = Cluster1_outSynapses{k};   % K_a x 3

        for l = 1:n
            B = Cluster1_outSynapses{l};  % K_b x 3

            % all pairwise distances
            D = pdist2(A, B, 'euclidean');   % K_a x K_b

            % 1) mean of all pairwise distances
            MeanPairDist(k,l) = mean(D(:),'omitnan');
        end
    end

    n_md = size(MeanPairDist,1);
    MeanPairDist(1:n_md+1:end) = NaN;  % set diagonal to NaN

    SuperClusterClosenessMatrix(i,i) = mean(MeanPairDist(:),'omitnan');

end

%% figure 6: supercluster-level outsynapse distance boxplot
% SuperClusterClosenessMatrix
M = SuperClusterClosenessMatrix;
n = size(M,1);

% Diagonal entries (within-cluster)
diag_vals = diag(M);
diag_vals = diag_vals(~isnan(diag_vals)); % remove NaN values

% Off-diagonal entries (between-cluster, deduplicated)
off_diag_vals = M(triu(true(n),1));
off_diag_vals = off_diag_vals(~isnan(off_diag_vals)); % remove NaN values

% --- Build boxplot data and groups ---
data_to_plot = [diag_vals; off_diag_vals];
group_indices = [ones(size(diag_vals)); ...
                 2 * ones(size(off_diag_vals))];
group_labels = {'Within-cluster', 'Between-cluster'};

% --- Boxplot ---
figure(6); clf; set(gcf, 'Color', 'w');
hold on;

% horizontal boxplot
boxplot(data_to_plot, group_indices, ...
    'Orientation', 'horizontal', ...
    'Labels', group_labels, ...
    'Notch', 'off', ...
    'Symbol','');

% labels and settings
ylim([0.5 2.5]);
xlim([0 2.5e-4]);
xlabel('distance um');
title('Super-Cluster-level outsynapse Distance');
set(gca, 'FontSize', 12, 'TickDir', 'out', 'Box', 'off');

% flip Y-axis so within-cluster is on top
set(gca, 'YDir','reverse');

set(gcf, 'Color', 'w');
% --- Rank-sum Test (Wilcoxon rank-sum test) ---
% d1: Within-cluster similarity
% d2: Between-cluster similarity
[p_val, h_stat, stats] = ranksum(diag_vals, off_diag_vals);

% print results
fprintf('\n--- Statistical Test Results ---\n');
fprintf('Method: Wilcoxon rank-sum test (Mann-Whitney U test)\n');
fprintf('P-value: %.4e\n', p_val);

if p_val < 0.05
    fprintf('Result: Statistically Significant (p < 0.05)\n');
else
    fprintf('Result: Not Statistically Significant\n');
end

% --- annotate the P-value on the plot (optional) ---
% place text at the top of the plot
text_str = sprintf('p = %.2e', p_val);
text(max(xlim)*0.7, 1.5, text_str, 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

%% Per-type cluster consistency
% Are neurons of the same type assigned to the same cluster?
types = ClusteringResult.type;         % e.g. {'LC9','LC9','LC10',...}
clusters = ClusteringResult.Cluster;   % e.g. [1,1,2,...]

unique_types = unique(types);
n_types = length(unique_types);

type_mean_ratios = nan(n_types, 1);  % mean consistency per type

for t = 1:n_types
    type_name = unique_types{t};
    idx = strcmp(types, type_name);      % indices of neurons of this type
    type_cluster = clusters(idx);        % cluster assignments of these neurons

    % cannot compare if only one neuron
    if sum(idx) <= 1
        continue;
    end

    % count pairs within the type that share a cluster
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

    % store mean consistency for this type
    type_mean_ratios(t) = count_same_cluster / total_pairs;
end

% number of neurons per type
type_counts = cellfun(@(x) sum(strcmp(types, x)), unique_types);

%% figure 7: consistency vs type size (scatter) and CDF, shared X-axis
figure(7);
set(gcf, 'Color', 'w');

% 'Intra-type clustering consistency' is the shared X-axis.

% --------------------
% Left Y-axis (Y1): neuron count (scatter)
% --------------------
yyaxis left;
scatter(type_mean_ratios, type_counts, 50, 'filled', ...
        'MarkerFaceColor', [0 0.4470 0.7410], ...
        'DisplayName', 'Type Size');
ylabel('Number of neurons per type');

% left Y-axis color and grid
ax = gca;
ax.YColor = [0 0.4470 0.7410]; % match left Y-axis color to scatter
grid on;

% --------------------
% Right Y-axis (Y2): cumulative probability (CDF)
% --------------------
yyaxis right;
% compute the empirical CDF
% x = consistency values (X), f = cumulative probability (Y)
[f, x] = ecdf(type_mean_ratios);

plot(x, f, 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2, ...
     'DisplayName', 'CDF');
ylabel('Cumulative Probability');

% right Y-axis color and range
ax.YColor = [0.8500 0.3250 0.0980]; % match right Y-axis color to CDF
ylim([0, 1]); % CDF Y-axis is always 0~1

% --------------------
% shared axis settings
% --------------------
xlabel('Intra-type clustering consistency');
title('Consistency vs. Type Size and Cumulative Distribution');
xlim([0-0.005, 1.005]); % shared X-axis 0~1
legend('Location', 'northwest'); % add legend
set(gca,'TickDir','out','Box','off')

function leafOrder = get_leaf_order_from_linkage(Z)
    % Z: linkage matrix (size: [n-1, 3])
    % return: leafOrder - reordered list of leaf indices

    n = size(Z, 1) + 1;
    total_nodes = 2 * n - 1;

    % record the leaf nodes contained in each cluster index
    cluster_leaves = cell(total_nodes, 1);

    % initial leaf nodes: each contains only itself
    for i = 1:n
        cluster_leaves{i} = i;
    end

    % build leaf sets by merging along the linkage rows
    for i = 1:n-1
        c1 = Z(i,1);
        c2 = Z(i,2);
        new_cluster = n + i;
        cluster_leaves{new_cluster} = [cluster_leaves{c1}; cluster_leaves{c2}];
    end

    % return the leaf order from the final tree root node
    leafOrder = cluster_leaves{2*n - 1};
end
