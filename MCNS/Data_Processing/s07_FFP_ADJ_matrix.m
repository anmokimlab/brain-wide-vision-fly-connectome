%% s07_FFP_ADJ_matrix (MCNS)
% MCNS analogue of the FAFB Data_Processing/s08_FFP_ADJ_matrix.m.
%
% Build the right-FFP -> central-brain output adjacency matrix: keep only the
% central-brain outputs by dropping connections inside the optic lobe (right +
% left) and outputs onto VPN neurons. Precursor for the Leiden clustering in
% Data_Processing/s08_FFP_leiden_clustering.ipynb.
clear all; close all; clc;

baseDir = fileparts(fileparts(mfilename('fullpath')));  % resolved from this script's location (portable; no absolute paths)

% RightFFP_NPIs from MCNS/Figures/fig_1_postPI_prePI_right_neurons.m
load(fullfile(baseDir, 'Processed_Data', 'right_neurons_thr0.mat'), 'RightFFP_NPIs')

opt = detectImportOptions(fullfile(baseDir, 'MCNS_Data', 'male-cns-v0.9-all-connections.csv'));
opt = setvartype(opt,'pre_root_id','int64');
opt = setvartype(opt,'post_root_id','int64');
MCNSConnections = readtable(fullfile(baseDir, 'MCNS_Data', 'male-cns-v0.9-all-connections.csv'),opt);

opt = detectImportOptions(fullfile(baseDir, 'MCNS_Data', 'male-cns-v0.9-classification.csv'));
opt = setvartype(opt,'root_id','int64');
MCNSClassification = readtable(fullfile(baseDir, 'MCNS_Data', 'male-cns-v0.9-classification.csv'),opt);

VPN_root_id = MCNSClassification.root_id(strcmp(MCNSClassification.super_class,'visual_projection'));

%% Build the output adjacency matrix of right FFP neurons
WantToSee = RightFFP_NPIs.root_id;
FromWant = ismember(MCNSConnections.pre_root_id,WantToSee);
WantConnection_Out = MCNSConnections(FromWant,:);
% Drop connections inside the optic lobe (right + left) and onto VPN neurons,
% keeping only the central-brain outputs.
OpticR = ismember(WantConnection_Out.neuropil,{'LA(R)', 'ME(R)','AME(R)', 'LO(R)','LOP(R)','LA(L)', 'ME(L)','AME(L)', 'LO(L)','LOP(L)'});
WantConnection_Out(OpticR,:) = [];
WantConnection_Out(ismember(WantConnection_Out.post_root_id,VPN_root_id),:) = [];

%% Making matrix
Post_Want = unique(WantConnection_Out.post_root_id);
WantMatrix_Out = zeros(size(WantToSee,1),size(Post_Want,1));

for i=1:1:size(WantConnection_Out,1)
    preIdx = find(WantToSee==WantConnection_Out.pre_root_id(i));
    postIdx = find(Post_Want==WantConnection_Out.post_root_id(i));
    WantMatrix_Out(preIdx,postIdx) = WantMatrix_Out(preIdx,postIdx)+WantConnection_Out.syn_count(i); % How to account for NT (neurotransmitter)?
end

%% Save outputs
writetable(RightFFP_NPIs, fullfile(baseDir,'Processed_Data','right_FFP_table.csv'))
csvwrite(fullfile(baseDir,'Processed_Data','right_FFP_output_matrix_CB_no_VPN_thr0.csv'),WantMatrix_Out)
writematrix(Post_Want, fullfile(baseDir,'Processed_Data','post_right_FFP_CB_no_VPN_thr0.csv'),'Delimiter',',')
save(fullfile(baseDir,'Processed_Data','right_FFP_output_matrix_CB_no_VPN_thr0.mat'),"WantMatrix_Out");
save(fullfile(baseDir,'Processed_Data','post_neurons_FFP_opticlobe_central_no_VPN_thr0.mat'),"Post_Want");
