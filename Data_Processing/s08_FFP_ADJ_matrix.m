clear all; close all; clc;

baseDir = fileparts(fileparts(mfilename('fullpath')));  % resolved from this script's location (portable; no absolute paths)

% FFP / FBP neuron classification (provides RightFFP_NPIs), saved by
% Figures/fig_1D_E_postPI_prePI_right_neurons.m
load(fullfile(baseDir, 'Processed_Data', 'right_neurons_thr0.mat'));

% Load Codex data
opt = detectImportOptions(fullfile(baseDir, 'Codex_Data', 'connections_princeton.csv'));
opt = setvartype(opt,'pre_root_id','int64');
opt = setvartype(opt,'post_root_id','int64');
FAFBConnections = readtable(fullfile(baseDir, 'Codex_Data', 'connections_princeton.csv'),opt);

opt = detectImportOptions(fullfile(baseDir, 'Codex_Data', 'classification.csv'));
opt = setvartype(opt,'root_id','int64');
FAFBClassification = readtable(fullfile(baseDir, 'Codex_Data', 'classification.csv'),opt);

Central_root_id = FAFBClassification.root_id(strcmp(FAFBClassification.super_class,'central'));
Optic_root_id   = FAFBClassification.root_id(strcmp(FAFBClassification.super_class,'optic'));
VPN_root_id     = FAFBClassification.root_id(strcmp(FAFBClassification.super_class,'visual_projection'));
VCN_root_id     = FAFBClassification.root_id(strcmp(FAFBClassification.super_class,'visual_centrifugal'));


%% Build the output adjacency matrix of right FFP neurons
WantToSee = RightFFP_NPIs.root_id;
FromWant = ismember(FAFBConnections.pre_root_id,WantToSee);
WantConnection_Out = FAFBConnections(FromWant,:);
% Drop connections inside the optic lobe (right + left) and onto VPN neurons,
% keeping only the central-brain outputs.
OpticR = ismember(WantConnection_Out.neuropil,{'LA_R', 'ME_R','AME_R', 'LO_R','LOP_R','LA_L', 'ME_L','AME_L', 'LO_L','LOP_L'});
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
