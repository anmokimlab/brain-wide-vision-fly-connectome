%% Load data / synapse coordinates / contralateral-L root ids
clear all; clc; close all

baseDir = fileparts(fileparts(mfilename('fullpath')));  % resolved from this script's location (portable; no absolute paths)

% Synapse coordinate table (x/y/z of every pre- and post-synaptic site)
opt = detectImportOptions(fullfile(baseDir, 'Codex_Data', 'fafb_v783_princeton_synapse_table.csv'));
opt = setvartype(opt,'pre_root_id_720575940','int64');
opt = setvartype(opt,'post_root_id_720575940','int64');
FAFB_synapse_coordinates = readtable(fullfile(baseDir, 'Codex_Data', 'fafb_v783_princeton_synapse_table.csv'),opt);
FAFB_synapse_coordinates.Properties.VariableNames(11) = "pre_root_id";
FAFB_synapse_coordinates.Properties.VariableNames(12) = "post_root_id";
FAFB_synapse_coordinates.pre_root_id=FAFB_synapse_coordinates.pre_root_id+int64(720575940000000000);
FAFB_synapse_coordinates.post_root_id=FAFB_synapse_coordinates.post_root_id+int64(720575940000000000);
FAFB_synapse_coordinates.pre_x=FAFB_synapse_coordinates.pre_x*1e-9;
FAFB_synapse_coordinates.pre_y=FAFB_synapse_coordinates.pre_y*1e-9;
FAFB_synapse_coordinates.pre_z=FAFB_synapse_coordinates.pre_z*1e-9;

opt = detectImportOptions(fullfile(baseDir, 'Codex_Data', 'consolidated_cell_types.csv'));
opt = setvartype(opt,'root_id','int64');
FAFBConsolidated_type = readtable(fullfile(baseDir, 'Codex_Data', 'consolidated_cell_types.csv'),opt);

% Optic-lobe neuropil meshes (faces are 0-indexed, vertices in nm -> m)
% Downloaded by Data_Processing/s01_fetch_neuropil_neuron.ipynb.
neuropilDir = fullfile(baseDir, 'Processed_Data', 'optic_lobe_neuropil_mesh');

LOP_R_Faces = readmatrix(fullfile(neuropilDir, 'LoP_R_faces.csv'));
LOP_R_Vertices=readmatrix(fullfile(neuropilDir, 'LoP_R_vertices.csv'));
LOP_R_Faces=LOP_R_Faces+1;
LOP_R_Vertices=LOP_R_Vertices*1e-9;

LOP_L_Faces = readmatrix(fullfile(neuropilDir, 'LoP_L_faces.csv'));
LOP_L_Vertices=readmatrix(fullfile(neuropilDir, 'LoP_L_vertices.csv'));
LOP_L_Faces=LOP_L_Faces+1;
LOP_L_Vertices=LOP_L_Vertices*1e-9;

LO_R_Faces = readmatrix(fullfile(neuropilDir, 'Lo_R_faces.csv'));
LO_R_Vertices=readmatrix(fullfile(neuropilDir, 'Lo_R_vertices.csv'));
LO_R_Faces=LO_R_Faces+1;
LO_R_Vertices=LO_R_Vertices*1e-9;

LO_L_Faces = readmatrix(fullfile(neuropilDir, 'Lo_L_faces.csv'));
LO_L_Vertices=readmatrix(fullfile(neuropilDir, 'Lo_L_vertices.csv'));
LO_L_Faces=LO_L_Faces+1;
LO_L_Vertices=LO_L_Vertices*1e-9;

ME_R_Faces = readmatrix(fullfile(neuropilDir, 'Me_R_faces.csv'));
ME_R_Vertices=readmatrix(fullfile(neuropilDir, 'Me_R_vertices.csv'));
ME_R_Faces=ME_R_Faces+1;
ME_R_Vertices=ME_R_Vertices*1e-9;

ME_L_Faces = readmatrix(fullfile(neuropilDir, 'Me_L_faces.csv'));
ME_L_Vertices=readmatrix(fullfile(neuropilDir, 'Me_L_vertices.csv'));
ME_L_Faces=ME_L_Faces+1;
ME_L_Vertices=ME_L_Vertices*1e-9;

%% Custom colormaps
% Number of colormap steps
n = 256;

% Black -> orange gradient (used for the innervation-depth heatmaps below)
startColor = [0, 0, 0];
endColor = [0.8500, 0.3250, 0.0980];
customColormapOut = zeros(n, 3);
for i = 1:3
    customColormapOut(:, i) = linspace(startColor(i), endColor(i), n);
end

%% Group the FBP neurons by cell type
% RightFBP_NPIs comes from the FFP/FBP/BDP classification saved by
% Figures/fig_1D_E_postPI_prePI_right_neurons.m
load(fullfile(baseDir, 'Processed_Data', 'right_neurons_thr0.mat'),'RightFBP_NPIs')

[RightFBP_type,~,ic]=unique(RightFBP_NPIs.type);

RightFBP_type=table(RightFBP_type,'VariableNames',{'type'});

for i=1:1:size(RightFBP_type,1)
    idx=ic==i;
    RightFBP_type.root_id{i}=RightFBP_NPIs.root_id(idx);
end

% Group FBP types by their target optic lobe (Me / Lo / Lop / Multi).
% (index lists into RightFBP_type, from the FBP output-neuropil classification, s10)
Me_FBP_idx=[4 22 25 26 29 30 31 32 33 68 69 70 71 72 73 74 75 76 78 81 82 83 84 85 86 87 93];
Lo_FBP_idx=[2 10 11 12 13 14 15 16 17 18 19 21 34 35 36 37 38 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 95];
Lop_FBP_idx=[3 7 8 9 28 62 63 64 65 66];
Multi_FBP_idx=[1 5 6 23 24 27 39 60 61 67 77 79 80 88 89 90 91 92 94 96];

FBP_Me=RightFBP_type(Me_FBP_idx,:);
FBP_Lo=RightFBP_type(Lo_FBP_idx,:);
FBP_Lop=RightFBP_type(Lop_FBP_idx,:);
FBP_Multi=RightFBP_type(Multi_FBP_idx,:);

%% Analysis configurations
% For each optic lobe, fit a reference layer (Ref) and profile the depth of the
% FBP group's OUTPUT (axon) synapses relative to that layer.
%   Me  : reference Dm6 dendrite,  target FBP_Me  axon
%   Lo  : reference LT1a dendrite, target FBP_Lo  axon
%   LoP : reference T5a  axon,     target FBP_Lop axon
configs = struct( ...
    'RefNeuron', {{'Dm6'},      {'LT1a'},     {'T5a'}}, ...
    'RefWhat',   {'dendrite',   'dendrite',   'axon'}, ...
    'RefWhere',  {'Me',         'Lo',         'LoP'}, ...
    'FBPgroup',  {FBP_Me,       FBP_Lo,       FBP_Lop}, ...
    'tagR',      {'FBP_Synapse_ME_R', 'FBP_Synapse_LO_R', 'FBP_Synapse_LOP_R'}, ...
    'tagL',      {'FBP_Synapse_ME_L', 'FBP_Synapse_LO_L', 'FBP_Synapse_LOP_L'});

% Depth bins shared across all configurations
edges = -7e-5:1e-6:5e-5;

% Output container; each field becomes a top-level variable in the saved .mat
FBP_output = struct();

close all
for c = 1:numel(configs)
    figBase = (c-1)*10;   % figure offset so all three runs stay on screen

    Ref_Neuron = configs(c).RefNeuron;
    Ref_What   = configs(c).RefWhat;
    Ref_Where  = configs(c).RefWhere;

    % Target = this group's FBP neurons, profiled by their output (axon) synapses
    grp = configs(c).FBPgroup;
    Target_Neuron = grp.type;
    Target_What   = 'axon';
    Target_Where  = Ref_Where;

    %% Reference and target root ids
    Ref_Table=table(Ref_Neuron,'VariableNames',{'cell_type'});
    for i=1:1:size(Ref_Table,1)
        idx=strcmpi(FAFBConsolidated_type.primary_type,Ref_Table.cell_type{i});
        Wantrootids=FAFBConsolidated_type.root_id(idx);
        Ref_Table.root_ids{i}=Wantrootids;
        Ref_Table.N{i}=sum(idx);
    end

    Target_Table=table(Target_Neuron,'VariableNames',{'cell_type'});
    Target_Table.root_ids = grp.root_id;   % root_ids already grouped in RightFBP_type

    %% Reference synapse locations
    for i=1:1:size(Ref_Table,1)
        ref_syn_loc=[];
        if strcmpi(Ref_What,'axon')
            Ref_Syn_pre_root_idx=find(ismember(FAFB_synapse_coordinates.pre_root_id,Ref_Table.root_ids{i}));
            ref_syn_loc= FAFB_synapse_coordinates(Ref_Syn_pre_root_idx,1:3);
        elseif strcmpi(Ref_What,'dendrite')
            Ref_Syn_post_root_idx=find(ismember(FAFB_synapse_coordinates.post_root_id,Ref_Table.root_ids{i}));
            ref_syn_loc=FAFB_synapse_coordinates(Ref_Syn_post_root_idx,1:3);
        else
            error('axon or dendrite')
        end

        ref_syn_loc=ref_syn_loc{:,:};
        if strcmpi(Ref_Where,'Lop')
            Ref_syn_loc_in_R=ref_syn_loc(intriangulation(LOP_R_Vertices,LOP_R_Faces,ref_syn_loc),:);
            Ref_syn_loc_in_L=ref_syn_loc(intriangulation(LOP_L_Vertices,LOP_L_Faces,ref_syn_loc),:);
        elseif strcmpi(Ref_Where,'Lo')
            Ref_syn_loc_in_R=ref_syn_loc(intriangulation(LO_R_Vertices,LO_R_Faces,ref_syn_loc),:);
            Ref_syn_loc_in_L=ref_syn_loc(intriangulation(LO_L_Vertices,LO_L_Faces,ref_syn_loc),:);
        elseif strcmpi(Ref_Where,'Me')
            Ref_syn_loc_in_R=ref_syn_loc(intriangulation(ME_R_Vertices,ME_R_Faces,ref_syn_loc),:);
            Ref_syn_loc_in_L=ref_syn_loc(intriangulation(ME_L_Vertices,ME_L_Faces,ref_syn_loc),:);
        else
            error('Lop or Lo or Me')
        end
        %%%% Remove outliers
        Ref_syn_loc_in_L_pointCloud = pointCloud(Ref_syn_loc_in_L);
        Ref_syn_loc_in_L_pointCloud_Denoise = pcdenoise(Ref_syn_loc_in_L_pointCloud,"NumNeighbors",20,"Threshold",2.5);
        Ref_syn_loc_in_L_Denoise=Ref_syn_loc_in_L_pointCloud_Denoise.Location;

        Ref_syn_loc_in_R_pointCloud = pointCloud(Ref_syn_loc_in_R);
        Ref_syn_loc_in_R_pointCloud_Denoise = pcdenoise(Ref_syn_loc_in_R_pointCloud,"NumNeighbors",20,"Threshold",2.5);
        Ref_syn_loc_in_R_Denoise=Ref_syn_loc_in_R_pointCloud_Denoise.Location;

        Ref_Table.syn_loc{i}=ref_syn_loc;
        Ref_Table.syn_loc_L{i}=Ref_syn_loc_in_L;
        Ref_Table.syn_loc_R{i}=Ref_syn_loc_in_R;
        Ref_Table.syn_loc_L_Denoise{i}=Ref_syn_loc_in_L_Denoise;
        Ref_Table.syn_loc_R_Denoise{i}=Ref_syn_loc_in_R_Denoise;
    end

    %% PCA of the reference layer
    Ref_x_mean_R=mean(Ref_Table.syn_loc_R_Denoise{1}(:,1));
    Ref_y_mean_R=mean(Ref_Table.syn_loc_R_Denoise{1}(:,2));
    Ref_z_mean_R=mean(Ref_Table.syn_loc_R_Denoise{1}(:,3));
    [coeffs_R, transformedData_R] = pca(Ref_Table.syn_loc_R_Denoise{1});

    Ref_x_mean_L=mean(Ref_Table.syn_loc_L_Denoise{1}(:,1));
    Ref_y_mean_L=mean(Ref_Table.syn_loc_L_Denoise{1}(:,2));
    Ref_z_mean_L=mean(Ref_Table.syn_loc_L_Denoise{1}(:,3));
    [coeffs_L, transformedData_L] = pca(Ref_Table.syn_loc_L_Denoise{1});
    reverse='false';
    if strcmpi(reverse,'true')
        transformedData_L(:,3) = -transformedData_L(:,3);   % flip only the side (R or L) that needs it
        coeffs_L(:,3) = -coeffs_L(:,3);
    end

    %% Fit a parabolic surface to the reference layer
    [xData_R, yData_R, zData_R] = prepareSurfaceData(transformedData_R(:,1), transformedData_R(:,2), transformedData_R(:,3));
    [xData_L, yData_L, zData_L] = prepareSurfaceData(transformedData_L(:,1), transformedData_L(:,2), transformedData_L(:,3));

    if strcmpi(Ref_Where,'Lop')
        ft = fittype( 'poly33' );
        [fitresult_R, gof_R] = fit( [xData_R, yData_R], zData_R, ft, 'Normalize', 'off'  );
        paraboloid_surface_R = @(x, y) fitresult_R.p00+fitresult_R.p10*x+fitresult_R.p01*y...
            +fitresult_R.p20*(x.^2)+fitresult_R.p11*x.*y+fitresult_R.p02*(y.^2)...
            +fitresult_R.p30*(x.^3)+fitresult_R.p21*(x.^2.*y)+fitresult_R.p12*(x.*y.^2)+fitresult_R.p03*(y.^3);
        fprintf('[%s] R-squared_R: %f\n', Ref_Where, gof_R.adjrsquare);

        [fitresult_L, gof_L] = fit( [xData_L, yData_L], zData_L, ft, 'Normalize', 'off'  );
        paraboloid_surface_L = @(x, y) fitresult_L.p00+fitresult_L.p10*x+fitresult_L.p01*y...
            +fitresult_L.p20*(x.^2)+fitresult_L.p11*x.*y+fitresult_L.p02*(y.^2)...
            +fitresult_L.p30*(x.^3)+fitresult_L.p21*(x.^2.*y)+fitresult_L.p12*(x.*y.^2)+fitresult_L.p03*(y.^3);
        fprintf('[%s] R-squared_L: %f\n', Ref_Where, gof_L.adjrsquare);

    elseif strcmpi(Ref_Where,'Lo') || strcmpi(Ref_Where,'Me')
        ft = fittype( 'poly22' );
        [fitresult_R, gof_R] = fit( [xData_R, yData_R], zData_R, ft, 'Normalize', 'off'  );
        paraboloid_surface_R = @(x, y) fitresult_R.p00+fitresult_R.p10*x+fitresult_R.p01*y...
            +fitresult_R.p20*x.^2+fitresult_R.p11*x.*y+fitresult_R.p02*y.^2;
        fprintf('[%s] R-squared_R: %f\n', Ref_Where, gof_R.adjrsquare);

        [fitresult_L, gof_L] = fit( [xData_L, yData_L], zData_L, ft, 'Normalize', 'off'  );
        paraboloid_surface_L = @(x, y) fitresult_L.p00+fitresult_L.p10*x+fitresult_L.p01*y...
            +fitresult_L.p20*x.^2+fitresult_L.p11*x.*y+fitresult_L.p02*y.^2;
        fprintf('[%s] R-squared_L: %f\n', Ref_Where, gof_L.adjrsquare);
    else
        error('Lop or Lo or Me')
    end

    %% Figure (figBase+1 & +2): right-hemisphere reference layer in PCA space
    figure(figBase+1);set(gcf,'Color','w');
    hold on;
    [X,Y] = meshgrid(min(transformedData_R(:,1)):5e-7:max(transformedData_R(:,1)),min(transformedData_R(:,2)):5e-7:max(transformedData_R(:,2)));
    Z = paraboloid_surface_R(X,Y);
    surf(X,Y,Z,'EdgeColor','none')
    scatter3(transformedData_R(:,1),transformedData_R(:,2),transformedData_R(:,3),20,[0.8500 0.3250 0.0980],'filled');
    grid on;

    if strcmpi(Ref_Where,'Me')
        ME_R_Vertices_PCA=ME_R_Vertices;
        ME_R_Vertices_PCA(:,1)=ME_R_Vertices_PCA(:,1)-Ref_x_mean_R;
        ME_R_Vertices_PCA(:,2)=ME_R_Vertices_PCA(:,2)-Ref_y_mean_R;
        ME_R_Vertices_PCA(:,3)=ME_R_Vertices_PCA(:,3)-Ref_z_mean_R;
        ME_R_Vertices_PCA=ME_R_Vertices_PCA*coeffs_R;
        trisurf(ME_R_Faces,ME_R_Vertices_PCA(:,1),ME_R_Vertices_PCA(:,2),ME_R_Vertices_PCA(:,3),'FaceAlpha',0.2,'EdgeColor','None');
    elseif strcmpi(Ref_Where,'Lo')
        LO_R_Vertices_PCA=LO_R_Vertices;
        LO_R_Vertices_PCA(:,1)=LO_R_Vertices_PCA(:,1)-Ref_x_mean_R;
        LO_R_Vertices_PCA(:,2)=LO_R_Vertices_PCA(:,2)-Ref_y_mean_R;
        LO_R_Vertices_PCA(:,3)=LO_R_Vertices_PCA(:,3)-Ref_z_mean_R;
        LO_R_Vertices_PCA=LO_R_Vertices_PCA*coeffs_R;
        trisurf(LO_R_Faces,LO_R_Vertices_PCA(:,1),LO_R_Vertices_PCA(:,2),LO_R_Vertices_PCA(:,3),'FaceAlpha',0.2,'EdgeColor','None');
    elseif strcmpi(Ref_Where,'Lop')
        LOP_R_Vertices_PCA=LOP_R_Vertices;
        LOP_R_Vertices_PCA(:,1)=LOP_R_Vertices_PCA(:,1)-Ref_x_mean_R;
        LOP_R_Vertices_PCA(:,2)=LOP_R_Vertices_PCA(:,2)-Ref_y_mean_R;
        LOP_R_Vertices_PCA(:,3)=LOP_R_Vertices_PCA(:,3)-Ref_z_mean_R;
        LOP_R_Vertices_PCA=LOP_R_Vertices_PCA*coeffs_R;
        trisurf(LOP_R_Faces,LOP_R_Vertices_PCA(:,1),LOP_R_Vertices_PCA(:,2),LOP_R_Vertices_PCA(:,3),'FaceAlpha',0.2,'EdgeColor','None');
    end

    axis("equal")
    xlabel('PC1','Interpreter','none'); ylabel('PC2','Interpreter','none'); zlabel('PC3','Interpreter','none');
    title(sprintf('%s ref - R', Ref_Where));
    view([-31.459612469645986,5.86855640905495])

    ax1 = gca;
    x_limits = xlim(ax1); y_limits = ylim(ax1); z_limits = zlim(ax1);

    figure(figBase+2);set(gcf,'Color','w');
    hold on;
    xyz_axes = eye(3);               % identity matrix (each column is one axis)
    arrow_length = 2e-5;
    new_axes = arrow_length * (xyz_axes * coeffs_R);   % pre-scaled
    quiver3(0, 0, 0, new_axes(1,1), new_axes(2,1), new_axes(3,1), 'r', 'LineWidth', 2); % X axis
    quiver3(0, 0, 0, new_axes(1,2), new_axes(2,2), new_axes(3,2), 'g', 'LineWidth', 2); % Y axis
    quiver3(0, 0, 0, new_axes(1,3), new_axes(2,3), new_axes(3,3), 'b', 'LineWidth', 2); % Z axis
    ax2 = gca;
    xlim(ax2, x_limits); ylim(ax2, y_limits); zlim(ax2, z_limits);
    xlabel('PC1','Interpreter','none'); ylabel('PC2','Interpreter','none'); zlabel('PC3','Interpreter','none');
    title(sprintf('%s ref axes - R', Ref_Where));
    view([-31.459612469645986,5.86855640905495])

    %% Figure (figBase+3 & +4): left-hemisphere reference layer in PCA space
    figure(figBase+3); set(gcf, 'Color', 'w');
    hold on;
    [X,Y] = meshgrid(min(transformedData_L(:,1)):5e-7:max(transformedData_L(:,1)), min(transformedData_L(:,2)):5e-7:max(transformedData_L(:,2)));
    Z = paraboloid_surface_L(X,Y);
    surf(X,Y,Z,'EdgeColor','none')
    scatter3(transformedData_L(:,1), transformedData_L(:,2), transformedData_L(:,3), 20, [0.8500 0.3250 0.0980], 'filled');
    grid on;

    if strcmpi(Ref_Where,'Me')
        ME_L_Vertices_PCA = ME_L_Vertices;
        ME_L_Vertices_PCA(:,1) = ME_L_Vertices_PCA(:,1) - Ref_x_mean_L;
        ME_L_Vertices_PCA(:,2) = ME_L_Vertices_PCA(:,2) - Ref_y_mean_L;
        ME_L_Vertices_PCA(:,3) = ME_L_Vertices_PCA(:,3) - Ref_z_mean_L;
        ME_L_Vertices_PCA = ME_L_Vertices_PCA * coeffs_L;
        trisurf(ME_L_Faces, ME_L_Vertices_PCA(:,1), ME_L_Vertices_PCA(:,2), ME_L_Vertices_PCA(:,3), 'FaceAlpha', 0.2, 'EdgeColor', 'None');
    elseif strcmpi(Ref_Where,'Lo')
        LO_L_Vertices_PCA = LO_L_Vertices;
        LO_L_Vertices_PCA(:,1) = LO_L_Vertices_PCA(:,1) - Ref_x_mean_L;
        LO_L_Vertices_PCA(:,2) = LO_L_Vertices_PCA(:,2) - Ref_y_mean_L;
        LO_L_Vertices_PCA(:,3) = LO_L_Vertices_PCA(:,3) - Ref_z_mean_L;
        LO_L_Vertices_PCA = LO_L_Vertices_PCA * coeffs_L;
        trisurf(LO_L_Faces, LO_L_Vertices_PCA(:,1), LO_L_Vertices_PCA(:,2), LO_L_Vertices_PCA(:,3), 'FaceAlpha', 0.2, 'EdgeColor', 'None');
    elseif strcmpi(Ref_Where,'Lop')
        LOP_L_Vertices_PCA = LOP_L_Vertices;
        LOP_L_Vertices_PCA(:,1) = LOP_L_Vertices_PCA(:,1) - Ref_x_mean_L;
        LOP_L_Vertices_PCA(:,2) = LOP_L_Vertices_PCA(:,2) - Ref_y_mean_L;
        LOP_L_Vertices_PCA(:,3) = LOP_L_Vertices_PCA(:,3) - Ref_z_mean_L;
        LOP_L_Vertices_PCA = LOP_L_Vertices_PCA * coeffs_L;
        trisurf(LOP_L_Faces, LOP_L_Vertices_PCA(:,1), LOP_L_Vertices_PCA(:,2), LOP_L_Vertices_PCA(:,3), 'FaceAlpha', 0.2, 'EdgeColor', 'None');
    end

    axis("equal")
    xlabel('PC1','Interpreter','none'); ylabel('PC2','Interpreter','none'); zlabel('PC3','Interpreter','none');
    title(sprintf('%s ref - L', Ref_Where));
    view([-31.459612469645986,5.86855640905495])

    ax1 = gca;
    x_limits = xlim(ax1); y_limits = ylim(ax1); z_limits = zlim(ax1);

    figure(figBase+4); set(gcf,'Color','w');
    hold on;
    xyz_axes = eye(3);
    arrow_length = 2e-5;
    new_axes = arrow_length * (xyz_axes * coeffs_L);
    quiver3(0, 0, 0, new_axes(1,1), new_axes(2,1), new_axes(3,1), 'r', 'LineWidth', 2);
    quiver3(0, 0, 0, new_axes(1,2), new_axes(2,2), new_axes(3,2), 'g', 'LineWidth', 2);
    quiver3(0, 0, 0, new_axes(1,3), new_axes(2,3), new_axes(3,3), 'b', 'LineWidth', 2);
    ax2 = gca;
    xlim(ax2, x_limits); ylim(ax2, y_limits); zlim(ax2, z_limits);
    xlabel('PC1','Interpreter','none'); ylabel('PC2','Interpreter','none'); zlabel('PC3','Interpreter','none');
    title(sprintf('%s ref axes - L', Ref_Where));
    view([-31.459612469645986,5.86855640905495])

    %% Target synapse locations
    Target_Table.syn_zDist_from_Parabola_R = cell(height(Target_Table),1);
    Target_Table.syn_zDist_from_Parabola_L = cell(height(Target_Table),1);

    for i=1:1:size(Target_Table,1)
        target_syn_loc=[];
        if strcmpi(Target_What,'axon')
            Target_Syn_pre_root_idx=find(ismember(FAFB_synapse_coordinates.pre_root_id,Target_Table.root_ids{i}));
            target_syn_loc=FAFB_synapse_coordinates(Target_Syn_pre_root_idx,1:3);
        elseif strcmpi(Target_What,'dendrite')
            Target_Syn_post_root_idx=find(ismember(FAFB_synapse_coordinates.post_root_id,Target_Table.root_ids{i}));
            target_syn_loc=FAFB_synapse_coordinates(Target_Syn_post_root_idx,1:3);
        else
            error('axon or dendrite')
        end
        if isempty(target_syn_loc)
            continue;
        end
        target_syn_loc=target_syn_loc{:,:};
        if strcmpi(Target_Where,'Lop')
            Target_syn_loc_in_R=target_syn_loc(intriangulation(LOP_R_Vertices,LOP_R_Faces,target_syn_loc),:);
            Target_syn_loc_in_L=target_syn_loc(intriangulation(LOP_L_Vertices,LOP_L_Faces,target_syn_loc),:);
        elseif strcmpi(Target_Where,'Lo')
            Target_syn_loc_in_R=target_syn_loc(intriangulation(LO_R_Vertices,LO_R_Faces,target_syn_loc),:);
            Target_syn_loc_in_L=target_syn_loc(intriangulation(LO_L_Vertices,LO_L_Faces,target_syn_loc),:);
        elseif strcmpi(Target_Where,'Me')
            Target_syn_loc_in_R=target_syn_loc(intriangulation(ME_R_Vertices,ME_R_Faces,target_syn_loc),:);
            Target_syn_loc_in_L=target_syn_loc(intriangulation(ME_L_Vertices,ME_L_Faces,target_syn_loc),:);
        else
            error('Lop or Lo or Me')
        end

        % Shift to the reference centroid and rotate into the reference PCA frame
        Target_syn_loc_in_L_Shift=Target_syn_loc_in_L;
        Target_syn_loc_in_L_Shift(:,1)=Target_syn_loc_in_L_Shift(:,1)-Ref_x_mean_L;
        Target_syn_loc_in_L_Shift(:,2)=Target_syn_loc_in_L_Shift(:,2)-Ref_y_mean_L;
        Target_syn_loc_in_L_Shift(:,3)=Target_syn_loc_in_L_Shift(:,3)-Ref_z_mean_L;
        Target_syn_loc_in_L_Shift_PCA=Target_syn_loc_in_L_Shift*coeffs_L;

        Target_syn_loc_in_R_Shift=Target_syn_loc_in_R;
        Target_syn_loc_in_R_Shift(:,1)=Target_syn_loc_in_R_Shift(:,1)-Ref_x_mean_R;
        Target_syn_loc_in_R_Shift(:,2)=Target_syn_loc_in_R_Shift(:,2)-Ref_y_mean_R;
        Target_syn_loc_in_R_Shift(:,3)=Target_syn_loc_in_R_Shift(:,3)-Ref_z_mean_R;
        Target_syn_loc_in_R_Shift_PCA=Target_syn_loc_in_R_Shift*coeffs_R;

        % Signed distance of each target synapse from the reference parabola (= innervation depth)
        z_hat_R=paraboloid_surface_R(Target_syn_loc_in_R_Shift_PCA(:,1),Target_syn_loc_in_R_Shift_PCA(:,2));
        Target_Table.syn_zDist_from_Parabola_R{i}=Target_syn_loc_in_R_Shift_PCA(:,3)-z_hat_R;

        z_hat_L=paraboloid_surface_L(Target_syn_loc_in_L_Shift_PCA(:,1),Target_syn_loc_in_L_Shift_PCA(:,2));
        Target_Table.syn_zDist_from_Parabola_L{i}=Target_syn_loc_in_L_Shift_PCA(:,3)-z_hat_L;
    end

    %% Bin the depth distribution (per target type), normalized to its own peak
    synapseBox_z_R=zeros(size(edges,2)-1,size(Target_Table,1));
    synapseBox_z_L=zeros(size(edges,2)-1,size(Target_Table,1));

    for i=1:1:size(Target_Table,1)
        temp=Target_Table.syn_zDist_from_Parabola_R{i};
        counts = histcounts(temp, 'BinEdges', edges);
        mx = max(counts); if mx==0, mx=1; end
        synapseBox_z_R(:,i)=counts / mx;

        temp=Target_Table.syn_zDist_from_Parabola_L{i};
        counts = histcounts(temp, 'BinEdges', edges);
        mx = max(counts); if mx==0, mx=1; end
        synapseBox_z_L(:,i)=counts / mx;
    end

    %% Figure (figBase+5 & +6): innervation-depth heatmaps (R / L)
    figure(figBase+5); set(gcf,'Color','w')
    imagesc(synapseBox_z_R)
    set(gca,'YTick',0.5:1:size(synapseBox_z_R,1)+0.5,'YTickLabel',edges,'XTick',1:1:size(Target_Table,1),'XTickLabel',Target_Table.cell_type,'TickDir','out','Box','off','YDir','normal')
    xlabel('Neurons'); ylabel('Innvervaiton depth (μm)');
    colormap(customColormapOut)
    title(sprintf('%s - R', Ref_Where))
    xtickangle(45)

    figure(figBase+6); set(gcf,'Color','w')
    imagesc(synapseBox_z_L)
    set(gca,'YTick',0.5:1:size(synapseBox_z_L,1)+0.5,'YTickLabel',edges,'XTick',1:1:size(Target_Table,1),'XTickLabel',Target_Table.cell_type,'TickDir','out','Box','off','YDir','normal')
    xlabel('Neurons'); ylabel('Innvervaiton depth (μm)');
    colormap(customColormapOut)
    title(sprintf('%s - L', Ref_Where))
    xtickangle(45)

    %% Store the depth-profile matrices under the neuropil-specific variable names
    FBP_output.(configs(c).tagR) = synapseBox_z_R;
    FBP_output.(configs(c).tagL) = synapseBox_z_L;
end

%% Save all six FBP output-synapse depth profiles into one .mat
% Variables: FBP_Synapse_ME_R/L, FBP_Synapse_LO_R/L, FBP_Synapse_LOP_R/L
save(fullfile(baseDir, 'Processed_Data', 'FBP_output_synapse_layer.mat'), '-struct', 'FBP_output');
