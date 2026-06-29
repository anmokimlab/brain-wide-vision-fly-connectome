% Figure 1B: render the LC9 neuron mesh together with its pre-/post-synaptic
% sites, and plot the distribution of those synapses along the x-axis.
clear all; clc; close all

baseDir = fileparts(fileparts(mfilename('fullpath')));  % resolved from this script's location (portable; no absolute paths)

% Synapse coordinate table (Codex). pre/post root_id columns are stored as
% offsets from 720575940000000000, so restore the full int64 root_id below.
opt = detectImportOptions(fullfile(baseDir, 'Codex_Data', 'fafb_v783_princeton_synapse_table.csv'));
opt = setvartype(opt,'pre_root_id_720575940','int64');
opt = setvartype(opt,'post_root_id_720575940','int64');
FAFB_synapse_coordinates = readtable(fullfile(baseDir, 'Codex_Data', 'fafb_v783_princeton_synapse_table.csv'),opt);
FAFB_synapse_coordinates.Properties.VariableNames(11) = "pre_root_id";
FAFB_synapse_coordinates.Properties.VariableNames(12) = "post_root_id";
FAFB_synapse_coordinates.pre_root_id=FAFB_synapse_coordinates.pre_root_id+int64(720575940000000000);
FAFB_synapse_coordinates.post_root_id=FAFB_synapse_coordinates.post_root_id+int64(720575940000000000);



%%

NeuronFaces = readmatrix(fullfile(baseDir, 'Processed_Data', 'fig1b_neuron_mesh', '720575940635224943_faces.csv'));
NeuronVertices=readmatrix(fullfile(baseDir, 'Processed_Data', 'fig1b_neuron_mesh', '720575940635224943_vertices.csv'));
NeuronFaces=NeuronFaces+1;
Neuron_root_id=int64(720575940635224943);


%%%

PreSynapse_x=FAFB_synapse_coordinates.pre_x(FAFB_synapse_coordinates.pre_root_id==Neuron_root_id);
PreSynapse_y=FAFB_synapse_coordinates.pre_y(FAFB_synapse_coordinates.pre_root_id==Neuron_root_id);
PreSynapse_z=FAFB_synapse_coordinates.pre_z(FAFB_synapse_coordinates.pre_root_id==Neuron_root_id);

Output_location=[PreSynapse_x  PreSynapse_y PreSynapse_z];

%%%
PostSynapse_x=FAFB_synapse_coordinates.post_x(FAFB_synapse_coordinates.post_root_id==Neuron_root_id);
PostSynapse_y=FAFB_synapse_coordinates.post_y(FAFB_synapse_coordinates.post_root_id==Neuron_root_id);
PostSynapse_z=FAFB_synapse_coordinates.post_z(FAFB_synapse_coordinates.post_root_id==Neuron_root_id);

Input_location=[PostSynapse_x  PostSynapse_y PostSynapse_z];

%%%
figure(1); set(gcf,'color','w')
trisurf(NeuronFaces, NeuronVertices(:,1),NeuronVertices(:,2),NeuronVertices(:,3),'EdgeColor','none','FaceAlpha',0.8,'FaceColor',[0.2 0.2 0.2]); hold on;
scatter3(Input_location(:,1),Input_location(:,2),Input_location(:,3),15,'filled','MarkerFaceAlpha',0.8)
scatter3(Output_location(:,1),Output_location(:,2),Output_location(:,3),15,'filled','MarkerFaceAlpha',0.8)
hold off;
axis equal
xlabel('x')
ylabel('y')
zlabel('z')
view([0 0 1]); 

%%%
minSynapse=min(min(Input_location(:,1)),min(Output_location(:,1)));
maxSynapse=max(max(Input_location(:,1)),max(Output_location(:,1)));

%%
figure(2); set(gcf,'color','w')
histogram(Input_location(:,1),linspace(minSynapse,maxSynapse,100),'EdgeColor','none','FaceAlpha',0.7); hold on;
histogram(Output_location(:,1),linspace(minSynapse,maxSynapse,100),'EdgeColor','none','FaceAlpha',0.7)
set(gca,'XColor','none','YColor','none','ZColor','none')


