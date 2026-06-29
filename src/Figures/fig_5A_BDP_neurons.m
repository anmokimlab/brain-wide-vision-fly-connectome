%% 1. Load data and initialize
clear all; clc; close all

baseDir = fileparts(fileparts(mfilename('fullpath')));  % resolved from this script's location (portable; no absolute paths)

opt = detectImportOptions(fullfile(baseDir, 'Codex_Data', 'fafb_v783_princeton_synapse_table.csv'));
opt = setvartype(opt,'pre_root_id_720575940','int64');
opt = setvartype(opt,'post_root_id_720575940','int64');
FAFB_synapse_coordinates = readtable(fullfile(baseDir, 'Codex_Data', 'fafb_v783_princeton_synapse_table.csv'),opt);

FAFB_synapse_coordinates.Properties.VariableNames(11) = "pre_root_id";
FAFB_synapse_coordinates.Properties.VariableNames(12) = "post_root_id";
FAFB_synapse_coordinates.pre_root_id  = FAFB_synapse_coordinates.pre_root_id  + int64(720575940000000000);
FAFB_synapse_coordinates.post_root_id = FAFB_synapse_coordinates.post_root_id + int64(720575940000000000);

%% 2. Neurons to plot
% Neuron meshes are downloaded by Data_Processing/s01_fetch_neuropil_neuron.ipynb
% into Processed_Data/fig5a_neuron_mesh.
% Each entry: name, root_id, and the per-neuron view angle [az el] in PCA space.
meshDir = fullfile(baseDir, 'Processed_Data', 'fig5a_neuron_mesh');

neurons = struct( ...
    'name',    {'LC9', 'LT43', 'LT52'}, ...
    'root_id', {int64(720575940635224943), int64(720575940629376528), int64(720575940631938647)}, ...
    'view',    {[3.769370935261506, -90], ...
                [1.748570307862659e+02, -9.692724476447463], ...
                [-0.538206501722826, -63.563181048454595]});

%% 3. Plot each neuron morphology (Figure 5A)
for n = 1:numel(neurons)
    Neuron_root_id = neurons(n).root_id;

    NeuronFaces    = readmatrix(fullfile(meshDir, sprintf('%d_faces.csv',    Neuron_root_id)));
    NeuronVertices = readmatrix(fullfile(meshDir, sprintf('%d_vertices.csv', Neuron_root_id)));
    NeuronFaces    = NeuronFaces + 1;   % MATLAB indexing correction

    % --- Synapse coordinates ---
    PreSynapse_x = FAFB_synapse_coordinates.pre_x(FAFB_synapse_coordinates.pre_root_id == Neuron_root_id);
    PreSynapse_y = FAFB_synapse_coordinates.pre_y(FAFB_synapse_coordinates.pre_root_id == Neuron_root_id);
    PreSynapse_z = FAFB_synapse_coordinates.pre_z(FAFB_synapse_coordinates.pre_root_id == Neuron_root_id);
    Output_location = [PreSynapse_x PreSynapse_y PreSynapse_z];

    PostSynapse_x = FAFB_synapse_coordinates.post_x(FAFB_synapse_coordinates.post_root_id == Neuron_root_id);
    PostSynapse_y = FAFB_synapse_coordinates.post_y(FAFB_synapse_coordinates.post_root_id == Neuron_root_id);
    PostSynapse_z = FAFB_synapse_coordinates.post_z(FAFB_synapse_coordinates.post_root_id == Neuron_root_id);
    Input_location = [PostSynapse_x PostSynapse_y PostSynapse_z];

    % --- Align mesh and synapses to the synapse-cloud PCA axes ---
    Ref_mean = mean([Input_location; Output_location], 1);
    coeffs   = pca([Input_location; Output_location]);

    NeuronVertices_PCA  = (NeuronVertices  - Ref_mean) * coeffs;
    Input_location_PCA  = (Input_location  - Ref_mean) * coeffs;
    Output_location_PCA = (Output_location - Ref_mean) * coeffs;

    % --- Draw ---
    figure(n); set(gcf,'color','w','Position',[444,50,1209,945])
    trisurf(NeuronFaces, NeuronVertices_PCA(:,1), NeuronVertices_PCA(:,2), NeuronVertices_PCA(:,3), ...
        'EdgeColor','none', 'FaceAlpha',0.8, 'FaceColor',[0.2 0.2 0.2]); hold on;
    scatter3(Input_location_PCA(:,1),  Input_location_PCA(:,2),  Input_location_PCA(:,3),  15, 'filled', 'MarkerFaceAlpha',0.8)
    scatter3(Output_location_PCA(:,1), Output_location_PCA(:,2), Output_location_PCA(:,3), 15, 'filled', 'MarkerFaceAlpha',0.8)
    hold off;
    axis equal
    xlabel('x'); ylabel('y'); zlabel('z');
    set(gca,'XColor','none','YColor','none','ZColor','none')

    % Constant grid spacing
    xL = xlim; yL = ylim; zL = zlim;
    tick_step = 10000;
    xticks(xL(1):tick_step:xL(2));
    yticks(yL(1):tick_step:yL(2));
    zticks(zL(1):tick_step:zL(2));

    view(neurons(n).view)   % per-neuron PCA-space view angle

end
