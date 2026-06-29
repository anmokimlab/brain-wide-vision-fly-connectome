%% 1. Load data and initialize
clear all; clc; close all

baseDir = fileparts(fileparts(mfilename('fullpath')));  % resolved from this script's location (portable; no absolute paths)

opt = detectImportOptions(fullfile(baseDir, 'Codex_Data', 'fafb_v783_princeton_synapse_table.csv'));
opt = setvartype(opt,'pre_root_id_720575940','int64');
opt = setvartype(opt,'post_root_id_720575940','int64');
FAFB_synapse_coordinates = readtable(fullfile(baseDir, 'Codex_Data', 'fafb_v783_princeton_synapse_table.csv'),opt);

FAFB_synapse_coordinates.Properties.VariableNames(11) = "pre_root_id";
FAFB_synapse_coordinates.Properties.VariableNames(12) = "post_root_id";
FAFB_synapse_coordinates.pre_root_id=FAFB_synapse_coordinates.pre_root_id+int64(720575940000000000);
FAFB_synapse_coordinates.post_root_id=FAFB_synapse_coordinates.post_root_id+int64(720575940000000000);

%% 2. Load and optimize the mesh data
close all
clearvars -except FAFB_synapse_coordinates baseDir
clc

% Neuron meshes are downloaded by Data_Processing/s01_fetch_neuropil_neuron.ipynb
% into Processed_Data/fig1h_neuron_mesh. Available neurons for Figure 1H:
%   720575940623047629  LPLC2
%   720575940642423797  cLP01
%   720575940635052569  LT52
%   720575940628406538  LT58
%   720575940633124143  PLP215
Neuron_root_id = int64(720575940628406538);   % select the neuron to plot (LT58)
meshDir = fullfile(baseDir, 'Processed_Data', 'fig1h_neuron_mesh');

NeuronFaces = readmatrix(fullfile(meshDir, sprintf('%d_faces.csv', Neuron_root_id)));
NeuronVertices = readmatrix(fullfile(meshDir, sprintf('%d_vertices.csv', Neuron_root_id)));
NeuronFaces = NeuronFaces + 1;   % MATLAB indexing correction


%% 2. Mesh simplification (decimation)
% Reduce the number of faces to ~2% (0.02). Adjust as needed (e.g. 0.05).
reduction_ratio = 0.02;
disp(['Original face count: ', num2str(size(NeuronFaces, 1))]);

[ReducedFaces, ReducedVertices] = reducepatch(NeuronFaces, NeuronVertices, reduction_ratio);

disp(['Face count after reduction: ', num2str(size(ReducedFaces, 1))]);
%
% % 3. Check the result (optional)
% patch('Faces', ReducedFaces, 'Vertices', ReducedVertices, ...
%       'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none');
% axis equal;

%% 3. Extract synapse coordinates
PreSynapse_x = FAFB_synapse_coordinates.pre_x(FAFB_synapse_coordinates.pre_root_id == Neuron_root_id);
PreSynapse_y = FAFB_synapse_coordinates.pre_y(FAFB_synapse_coordinates.pre_root_id == Neuron_root_id);
PreSynapse_z = FAFB_synapse_coordinates.pre_z(FAFB_synapse_coordinates.pre_root_id == Neuron_root_id);
Output_location = [PreSynapse_x PreSynapse_y PreSynapse_z];

PostSynapse_x = FAFB_synapse_coordinates.post_x(FAFB_synapse_coordinates.post_root_id == Neuron_root_id);
PostSynapse_y = FAFB_synapse_coordinates.post_y(FAFB_synapse_coordinates.post_root_id == Neuron_root_id);
PostSynapse_z = FAFB_synapse_coordinates.post_z(FAFB_synapse_coordinates.post_root_id == Neuron_root_id);
Input_location = [PostSynapse_x PostSynapse_y PostSynapse_z];

%% 4. Figure 1: neuron and synapse visualization
figure(1); set(gcf,'color','w')
% Draw using the reduced faces and vertices.
trisurf(ReducedFaces, ReducedVertices(:,1), ReducedVertices(:,2), ReducedVertices(:,3), ...
    'EdgeColor','none', 'FaceAlpha',0.8, 'FaceColor',[0.2 0.2 0.2]);
hold on;
%
scatter3(Input_location(:,1), Input_location(:,2), Input_location(:,3), 15, 'filled', 'MarkerFaceAlpha', 0.8)
scatter3(Output_location(:,1), Output_location(:,2), Output_location(:,3), 15, 'filled', 'MarkerFaceAlpha', 0.8)

hold off;
axis equal; xlabel('x'); ylabel('y'); zlabel('z');
view([1.192345471906226e+02 -57.951351351351342])

% Save as EPS (vector)

%% 5. Transformation matrix and histogram
target_az = 1.192345471906226e+02;
target_el = -57.951351351351342;
T = viewmtx(target_az, target_el);

Input_homo = [Input_location, ones(size(Input_location, 1), 1)];
Output_homo = [Output_location, ones(size(Output_location, 1), 1)];

Input_screen = Input_homo * T;
Output_screen = Output_homo * T;
Input_screen_y = Input_screen(:, 2);
Output_screen_y = Output_screen(:, 2);

figure(2); clf; set(gcf,'color','w')
histogram(Input_screen_y, 100, 'EdgeColor','none','FaceAlpha',0.7, 'DisplayName', 'Input'); hold on;
histogram(Output_screen_y, 100, 'EdgeColor','none','FaceAlpha',0.7, 'DisplayName', 'Output');
legend('show');
title(sprintf('Screen-Y Distribution (Az: %.1f, El: %.1f)', target_az, target_el));
xlabel('Screen Vertical Position'); ylabel('Count');

%% 6. Direction-axis visualization
figure(3); clf; set(gcf,'color','w')
origin = [0 0 0 1];
axes_3d = [1 0 0 1; 0 1 0 1; 0 0 1 1];
origin_screen = origin * T;
axes_screen = axes_3d * T;
dirs = axes_screen(:, 1:2) - origin_screen(1:2);
colors = [1 0 0; 0 1 0; 0 0 1];
axis_names = {'X-axis', 'Y-axis', 'Z-axis'};

for i = 1:3
    quiver(0, 0, dirs(i,1), dirs(i,2), 0, 'Color', colors(i,:), 'LineWidth', 2, 'MaxHeadSize', 0.5);
    hold on;
    text(dirs(i,1)*1.2, dirs(i,2)*1.2, axis_names{i}, 'Color', colors(i,:), 'FontWeight', 'bold');
end
grid on; axis equal; xlim([-1.5 1.5]); ylim([-1.5 1.5]);

%%


%% 1. Project the 3D coordinates to 2D for the current view angle
az = 119.2345;
el = -57.9514;

% Get the rotation matrix for the current view (angle-based)
view_mat = viewmtx(az, el);

% Project the neuron vertices (4D homogeneous-coordinate transform)
V4 = [NeuronVertices, ones(size(NeuronVertices,1), 1)] * view_mat';
V2d = V4(:, 1:2);   % projected 2D x, y coordinates

% Project the synapse positions
In4 = [Input_location, ones(size(Input_location,1), 1)] * view_mat';
Out4 = [Output_location, ones(size(Output_location,1), 1)] * view_mat';
In2d = In4(:, 1:2);
Out2d = Out4(:, 1:2);

%% 2. Extract the silhouette (substitute for the Illustrator pathfinder)
% alphaShape builds a boundary that wraps the point cloud.
% Adjust the 'Alpha' value (radius) to control the detail: too large -> blunt,
% too small -> holes appear.
shp = alphaShape(V2d(:,1), V2d(:,2), 10);   % 10 may need tuning to the data scale
[~, V_boundary] = boundaryFacets(shp);
