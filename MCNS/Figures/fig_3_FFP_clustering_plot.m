%% Load data
% MCNS analogue of the FAFB Figures/fig_3C_D_E_F_G_FFP_clustering_plot.m.
% Plots the FFP clustering result: the raw-similarity dendrogram and the
% inter-cluster connectivity matrix (both reordered by dendrogram leaf order).
clear all; close all; clc;

baseDir = fileparts(fileparts(mfilename('fullpath')));  % resolved from this script's location (portable; no absolute paths)

% Leiden cluster assignment of the FFP neurons (rows), produced by
% Data_Processing/s08_FFP_leiden_clustering.ipynb (singleton-row clusters removed)
opt = detectImportOptions(fullfile(baseDir,'Processed_Data','leiden_right_FFP_output_CB_no_VPN_thr0_resol3.6_100000_row_singletons_removed.csv'));
opt = setvartype(opt,'root_id','int64');
ClusteringResult=readtable(fullfile(baseDir,'Processed_Data','leiden_right_FFP_output_CB_no_VPN_thr0_resol3.6_100000_row_singletons_removed.csv'),opt);

% MCNS connectivity
opt = detectImportOptions(fullfile(baseDir,'MCNS_Data','male-cns-v0.9-all-connections.csv'));
opt = setvartype(opt,'pre_root_id','int64');
opt = setvartype(opt,'post_root_id','int64');
MCNSConnections = readtable(fullfile(baseDir,'MCNS_Data','male-cns-v0.9-all-connections.csv'),opt);

opt = detectImportOptions(fullfile(baseDir,'MCNS_Data','male-cns-v0.9-primary-types.csv'));
opt = setvartype(opt,'root_id','int64');
MCNSConsolidatedTypes = readtable(fullfile(baseDir,'MCNS_Data','male-cns-v0.9-primary-types.csv'),opt);

% Inter-cluster (bipartite) connectivity matrix, produced by s08
opt = detectImportOptions(fullfile(baseDir,'Processed_Data','leiden_right_FFP_intercluster_connectivity_CB_no_VPN_thr0_resol3.6_100000_row_singletons_removed.csv'));
BipartiteConnectivity = readmatrix(fullfile(baseDir,'Processed_Data','leiden_right_FFP_intercluster_connectivity_CB_no_VPN_thr0_resol3.6_100000_row_singletons_removed.csv'),opt);
BipartiteConnectivity = reshape(BipartiteConnectivity,[max(ClusteringResult.Cluster)+1,max(ClusteringResult.Cluster)+1]);

% Raw-similarity dendrogram (linkage matrix), produced by s08
opt = detectImportOptions(fullfile(baseDir,'Processed_Data','dendrogram_raw_CB_no_VPN_thr0_resol3.6_100000_row_singletons_removed.csv'));
RightFFP_Dendrogram = readmatrix(fullfile(baseDir,'Processed_Data','dendrogram_raw_CB_no_VPN_thr0_resol3.6_100000_row_singletons_removed.csv'),opt);
RightFFP_Dendrogram(:,4)=[];
RightFFP_Dendrogram(:,1)=RightFFP_Dendrogram(:,1)+1;
RightFFP_Dendrogram(:,2)=RightFFP_Dendrogram(:,2)+1;

% Post-neuron root_ids (provides Post_Want), produced by s07
load(fullfile(baseDir,'Processed_Data','post_neurons_FFP_opticlobe_central_no_VPN_thr0.mat'))
% Kept post-neuron columns after singleton-row removal, produced by s08
KeptCols = readtable(fullfile(baseDir,'Processed_Data','leiden_post_right_FFP_kept_cols_CB_no_VPN_thr0_100000.csv'));
% Python indices are 0-based, so add 1 for MATLAB
col_idx = KeptCols.original_index + 1;
Post_Want_Filtered = Post_Want(col_idx);

% Leiden cluster assignment of the post neurons (columns), produced by s08
ClusteringResult_PostNeurons=readtable(fullfile(baseDir,'Processed_Data','leiden_post_right_FFP_CB_no_VPN_thr0_100000_row_singletons_removed.csv'));
ClusteringResult_PostNeurons.Properties.VariableNames(1) = "root_id";
ClusteringResult_PostNeurons.Properties.VariableNames(2) = "Cluster";
ClusteringResult_PostNeurons.root_id = Post_Want_Filtered;

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
    idx_Consol=find(MCNSConsolidatedTypes.root_id==ClusteringResult_PostNeurons.root_id(i));
    if ~isempty(idx_Consol)
        ClusteringResult_PostNeurons.type{i}=MCNSConsolidatedTypes.primary_type{idx_Consol};
    else
        ClusteringResult_PostNeurons.type{i}='';
    end
end

for i=1:1:size(ClusteringValueCounts_PostNeurons,1)
    AllType=ClusteringResult_PostNeurons.type(ClusteringValueCounts_PostNeurons{i,3});
    UniqueType=unique(AllType);
    current_FFP_Cluster=ClusteringValueCounts{i,6};
    for j=1:1:size(UniqueType,1)
        UniqueType{j,2}=sum(strcmp(AllType,UniqueType{j,1}));
        UniqueType{j,3}=ClusteringValueCounts_PostNeurons{i,5}(strcmp(AllType,UniqueType{j,1}));
        UniqueType{j,4}=sum(MCNSConnections.syn_count(ismember(MCNSConnections.post_root_id,UniqueType{j,3})&ismember(MCNSConnections.pre_root_id,ClusteringValueCounts{i, 5})));
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
    'Reorder', leafOrder, 'ColorThreshold', 0.51, 'Labels', temp_RightFFP_labels);
set(H,'LineWidth',2);
set(gca,'TickDir','out');
ylim([0 5]);

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

%% Output neuropil targets per supercluster
% Group the 32 FFP clusters into 11 hand-defined superclusters and, for each,
% tabulate the central-brain output synapses by neuropil (optic-lobe neuropils
% and VPN targets excluded, matching the matrix construction).
SuperClusterOrder{1}=[2;20;21;27;30];
SuperClusterOrder{2}=[25;29];
SuperClusterOrder{3}=[16;19;31];
SuperClusterOrder{4}=[28];
SuperClusterOrder{5}=[23];
SuperClusterOrder{6}=[5;7;17];
SuperClusterOrder{7}=[3;11;32];
SuperClusterOrder{8}=[8;9;12;13;14;18;22];
SuperClusterOrder{9}=[1;4;6;24];
SuperClusterOrder{10}=[15;26];
SuperClusterOrder{11}=[10];

opt = detectImportOptions(fullfile(baseDir,'MCNS_Data','male-cns-v0.9-classification.csv'));
opt = setvartype(opt,'root_id','int64');
MCNSClassification = readtable(fullfile(baseDir,'MCNS_Data','male-cns-v0.9-classification.csv'),opt);
VPN_root_id=MCNSClassification.root_id(strcmp(MCNSClassification.super_class,'visual_projection'));

SuperCluster_neuropil_targets = cell(size(SuperClusterOrder,2),1);
% Per-supercluster post-synaptic targets
SuperCluster_targets = struct('n_post_neurons',{}, 'post_type_counts',{}, ...
    'post_neuropil_syn',{}, 'post_type_syn',{});
for i=1:1:size(SuperClusterOrder,2)
    current_root_ids=ClusteringValueCounts(SuperClusterOrder{i},5);
    current_root_ids = vertcat(current_root_ids{:});

    post_idx=ismember(MCNSConnections.pre_root_id,current_root_ids);
    post_Connections=MCNSConnections(post_idx,:);
    % Same filtering as the matrix construction
    OpticR=ismember(post_Connections.neuropil,{'LA(R)', 'ME(R)','AME(R)', 'LO(R)','LOP(R)','LA(L)', 'ME(L)','AME(L)', 'LO(L)','LOP(L)'});
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
    [tf,loc] = ismember(postNeurons, MCNSConsolidatedTypes.root_id);
    postTypes(tf)  = string(MCNSConsolidatedTypes.primary_type(loc(tf)));
    postTypes(~tf) = "Untyped";

    [uTypes,~,ic_t] = unique(postTypes);
    SuperCluster_targets(i).n_post_neurons    = numel(postNeurons);
    SuperCluster_targets(i).post_type_counts  = sortrows([cellstr(uTypes), num2cell(accumarray(ic_t,1))],          2, 'descend'); % {type, neuron count}
    SuperCluster_targets(i).post_neuropil_syn = Post_neuropils;                                                                    % {neuropil, synapse count}
    SuperCluster_targets(i).post_type_syn     = sortrows([cellstr(uTypes), num2cell(accumarray(ic_t,synPerPost))], 2, 'descend');  % {type, synapse count}
end

save(fullfile(baseDir,'Processed_Data','FFP_supercluster_targets.mat'), 'SuperCluster_targets');

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
