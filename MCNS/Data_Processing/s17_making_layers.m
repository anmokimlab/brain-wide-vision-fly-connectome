%% s17 — Reference-layer PCA bases for the optic-lobe neuropils  [MCNS]
% MCNS analogue of the FAFB Data_Processing/s17_making_layers.m.
% For each reference cell type (Me = Dm6 dendrite, Lo = LT1a dendrite, LoP = T5a axon),
% it collects the reference synapses inside the target neuropil (right and left), denoises
% them, runs PCA, and fits a parabolic surface (poly22 for Me/Lo, poly33 for LoP). The
% per-hemisphere mean, PCA coefficients, and surface fit are stored in PCA_Basis and saved
% to Processed_Data/neuropil_PCA_basis.mat, read by Figures/fig_6A_BLP_RFs_PFs.m.
%
% Unlike the FAFB version (which assigns synapses to a neuropil via mesh intriangulation),
% the MCNS synapse table carries a `neuropil` column, so synapses are assigned to a region
% by that column; the ROI meshes are used only for the visualization.

clear all; clc; close all

baseDir = fileparts(fileparts(mfilename('fullpath')));  % resolved from this script's location (portable; no absolute paths)

%% Load data
% Synapse coordinates (8 nm voxels -> metres; columns 1:3 = pre_x/pre_y/pre_z)
opt = detectImportOptions(fullfile(baseDir, 'MCNS_Data', 'male-cns-v0.9-synapse-coordinates.csv'));
opt = setvartype(opt,'pre_root_id','int64');
opt = setvartype(opt,'post_root_id','int64');
MCNS_synapse_coordinates = readtable(fullfile(baseDir, 'MCNS_Data', 'male-cns-v0.9-synapse-coordinates.csv'),opt);
MCNS_synapse_coordinates.pre_x = MCNS_synapse_coordinates.pre_x*8e-9;
MCNS_synapse_coordinates.pre_y = MCNS_synapse_coordinates.pre_y*8e-9;
MCNS_synapse_coordinates.pre_z = MCNS_synapse_coordinates.pre_z*8e-9;

opt = detectImportOptions(fullfile(baseDir, 'MCNS_Data', 'male-cns-v0.9-primary-types.csv'));
opt = setvartype(opt,'root_id','int64');
MCNSConsolidatedTypes = readtable(fullfile(baseDir, 'MCNS_Data', 'male-cns-v0.9-primary-types.csv'),opt);

% Optic-lobe ROI meshes (faces 0-indexed, vertices in 8 nm voxels -> m),
% exported by Data_Processing/s10_export_roi_mesh_csv.py.
meshDir = fullfile(baseDir, 'Processed_Data', 'roi_mesh_csv');
[LOP_R_Faces, LOP_R_Vertices] = load_roi_mesh(meshDir, 'LOP(R)');
[LOP_L_Faces, LOP_L_Vertices] = load_roi_mesh(meshDir, 'LOP(L)');
[LO_R_Faces,  LO_R_Vertices]  = load_roi_mesh(meshDir, 'LO(R)');
[LO_L_Faces,  LO_L_Vertices]  = load_roi_mesh(meshDir, 'LO(L)');
[ME_R_Faces,  ME_R_Vertices]  = load_roi_mesh(meshDir, 'ME(R)');
[ME_L_Faces,  ME_L_Vertices]  = load_roi_mesh(meshDir, 'ME(L)');

%% Batch over the three reference layers (Me / Lo / LoP)
BatchConfigs = { ...
    'Dm6',  'dendrite', 'Me' ;  ...   % Medulla reference
    'LT1a', 'dendrite', 'Lo' ;  ...   % Lobula reference
    'T5a',  'axon',     'Lop'};       % Lobula plate reference
reversePC3 = {'true', 'false', 'false'};   % flip PC3 (both hemispheres) for this config

PCA_Basis = struct();

for bi = 1:size(BatchConfigs,1)
    Ref_Neuron = {BatchConfigs{bi,1}};   % e.g. {'Dm6'}
    Ref_What   = BatchConfigs{bi,2};     % 'dendrite' or 'axon'
    Ref_Where  = BatchConfigs{bi,3};     % 'Me' / 'Lo' / 'Lop' (strcmpi)

    fprintf('\n=== [%d/%d] Ref: %s | What: %s | Where: %s ===\n', ...
        bi, size(BatchConfigs,1), Ref_Neuron{1}, Ref_What, Ref_Where);

    % --- Reference neuron root_ids ---
    Ref_Table = table(Ref_Neuron, 'VariableNames', {'cell_type'});
    for i = 1:size(Ref_Table,1)
        idx = strcmpi(MCNSConsolidatedTypes.primary_type, strtrim(Ref_Table.cell_type{i}));
        Ref_Table.root_ids{i} = MCNSConsolidatedTypes.root_id(idx);
        Ref_Table.N{i} = sum(idx);
    end

    % --- Reference synapse locations inside the target neuropil (R / L), denoised ---
    % Synapses are assigned to a neuropil by the synapse table's `neuropil` column.
    for i = 1:size(Ref_Table,1)
        if strcmpi(strtrim(Ref_What),'axon')
            ref_syn_loc = MCNS_synapse_coordinates(ismember(MCNS_synapse_coordinates.pre_root_id, Ref_Table.root_ids{i}), :);
        elseif strcmpi(strtrim(Ref_What),'dendrite')
            ref_syn_loc = MCNS_synapse_coordinates(ismember(MCNS_synapse_coordinates.post_root_id, Ref_Table.root_ids{i}), :);
        else
            error('axon or dendrite');
        end

        if strcmpi(Ref_Where,'Lop')
            Ref_syn_loc_in_R = ref_syn_loc(strcmpi(ref_syn_loc.neuropil,'LOP(R)'),:);
            Ref_syn_loc_in_L = ref_syn_loc(strcmpi(ref_syn_loc.neuropil,'LOP(L)'),:);
        elseif strcmpi(Ref_Where,'Lo')
            Ref_syn_loc_in_R = ref_syn_loc(strcmpi(ref_syn_loc.neuropil,'LO(R)'),:);
            Ref_syn_loc_in_L = ref_syn_loc(strcmpi(ref_syn_loc.neuropil,'LO(L)'),:);
        elseif strcmpi(Ref_Where,'Me')
            Ref_syn_loc_in_R = ref_syn_loc(strcmpi(ref_syn_loc.neuropil,'ME(R)'),:);
            Ref_syn_loc_in_L = ref_syn_loc(strcmpi(ref_syn_loc.neuropil,'ME(L)'),:);
        else
            error('Lop or Lo or Me');
        end

        Ref_syn_loc_in_R = Ref_syn_loc_in_R{:,1:3};
        Ref_syn_loc_in_L = Ref_syn_loc_in_L{:,1:3};

        % Outlier removal
        Ref_syn_loc_in_L_Denoise = pcdenoise(pointCloud(Ref_syn_loc_in_L), 'NumNeighbors',20, 'Threshold',2.5).Location;
        Ref_syn_loc_in_R_Denoise = pcdenoise(pointCloud(Ref_syn_loc_in_R), 'NumNeighbors',20, 'Threshold',2.5).Location;

        Ref_Table.syn_loc_L{i}         = Ref_syn_loc_in_L;
        Ref_Table.syn_loc_R{i}         = Ref_syn_loc_in_R;
        Ref_Table.syn_loc_L_Denoise{i} = Ref_syn_loc_in_L_Denoise;
        Ref_Table.syn_loc_R_Denoise{i} = Ref_syn_loc_in_R_Denoise;
    end

    % --- PCA per hemisphere ---
    Ref_mean_R = mean(Ref_Table.syn_loc_R_Denoise{1}, 1);
    [coeffs_R, transformedData_R] = pca(Ref_Table.syn_loc_R_Denoise{1});

    Ref_mean_L = mean(Ref_Table.syn_loc_L_Denoise{1}, 1);
    [coeffs_L, transformedData_L] = pca(Ref_Table.syn_loc_L_Denoise{1});

    if strcmpi(reversePC3{bi},'true')
        transformedData_L(:,3) = -transformedData_L(:,3);   coeffs_L(:,3) = -coeffs_L(:,3);
        transformedData_R(:,3) = -transformedData_R(:,3);   coeffs_R(:,3) = -coeffs_R(:,3);
    end

    % --- Parabolic surface fit (poly33 for LoP, poly22 for Me/Lo) ---
    [xData_R, yData_R, zData_R] = prepareSurfaceData(transformedData_R(:,1), transformedData_R(:,2), transformedData_R(:,3));
    [xData_L, yData_L, zData_L] = prepareSurfaceData(transformedData_L(:,1), transformedData_L(:,2), transformedData_L(:,3));

    if strcmpi(Ref_Where,'Lop')
        ft = fittype('poly33');
    elseif strcmpi(Ref_Where,'Lo') || strcmpi(Ref_Where,'Me')
        ft = fittype('poly22');
    else
        error('Lop or Lo or Me');
    end

    [fitresult_R, gof_R] = fit([xData_R, yData_R], zData_R, ft, 'Normalize','off');
    [fitresult_L, gof_L] = fit([xData_L, yData_L], zData_L, ft, 'Normalize','off');
    fprintf('R-squared_R: %f\n', gof_R.adjrsquare);
    fprintf('R-squared_L: %f\n', gof_L.adjrsquare);

    paraboloid_surface_R = make_surface(fitresult_R);
    paraboloid_surface_L = make_surface(fitresult_L);

    % --- Visualization (right then left): fitted surface + PCA-aligned neuropil mesh ---
    switch upper(Ref_Where)
        case 'ME',  Faces_R = ME_R_Faces;  Verts_R = ME_R_Vertices;  Faces_L = ME_L_Faces;  Verts_L = ME_L_Vertices;
        case 'LO',  Faces_R = LO_R_Faces;  Verts_R = LO_R_Vertices;  Faces_L = LO_L_Faces;  Verts_L = LO_L_Vertices;
        case 'LOP', Faces_R = LOP_R_Faces; Verts_R = LOP_R_Vertices; Faces_L = LOP_L_Faces; Verts_L = LOP_L_Vertices;
    end

    close all
    figure(1); set(gcf,'Color','w');
    plot_layer(transformedData_R, paraboloid_surface_R, Faces_R, Verts_R, Ref_mean_R, coeffs_R);

    figure(2); set(gcf,'Color','w');
    plot_layer(transformedData_L, paraboloid_surface_L, Faces_L, Verts_L, Ref_mean_L, coeffs_L);

    % --- Store basis (mean, coeffs, fit, gof) ---
    key = matlab.lang.makeValidName(sprintf('%s__%s__%s', upper(Ref_Where), strtrim(Ref_Neuron{1}), lower(strtrim(Ref_What))));
    PCA_Basis.(key).Ref_Neuron      = strtrim(Ref_Neuron{1});
    PCA_Basis.(key).Ref_What        = strtrim(Ref_What);
    PCA_Basis.(key).Ref_Where       = Ref_Where;
    PCA_Basis.(key).reverseLeftPC3  = strcmpi(reversePC3{bi},'true');   % both hemispheres flipped
    PCA_Basis.(key).reverseRightPC3 = strcmpi(reversePC3{bi},'true');

    PCA_Basis.(key).R.mean   = Ref_mean_R;
    PCA_Basis.(key).R.coeffs = coeffs_R;
    PCA_Basis.(key).L.mean   = Ref_mean_L;
    PCA_Basis.(key).L.coeffs = coeffs_L;

    PCA_Basis.(key).fit.R = fitresult_R;
    PCA_Basis.(key).fit.L = fitresult_L;
    PCA_Basis.(key).gof.R = gof_R;
    PCA_Basis.(key).gof.L = gof_L;

    fprintf('Saved basis in struct field: PCA_Basis.%s\n', key);
end

save(fullfile(baseDir, 'Processed_Data', 'neuropil_PCA_basis.mat'), 'PCA_Basis');
fprintf('\n>> All bases saved to neuropil_PCA_basis.mat\n');

%% ====================== Local functions ======================
function [F, V] = load_roi_mesh(meshDir, roi)
% Read a neuprint ROI mesh; faces are 0-indexed (corrected to 1-based) and
% vertices are in 8 nm voxels (converted to metres).
F = readmatrix(fullfile(meshDir, roi, 'faces.csv')) + 1;
V = readmatrix(fullfile(meshDir, roi, 'vertices.csv')) * 8e-9;
end

function s = make_surface(fitresult)
% Anonymous z(x,y) for a poly22 or poly33 surface fit (detected from coefficients).
if any(strcmp(coeffnames(fitresult), 'p30'))
    s = @(x,y) fitresult.p00 + fitresult.p10*x + fitresult.p01*y ...
        + fitresult.p20*x.^2 + fitresult.p11*x.*y + fitresult.p02*y.^2 ...
        + fitresult.p30*x.^3 + fitresult.p21*x.^2.*y + fitresult.p12*x.*y.^2 + fitresult.p03*y.^3;
else
    s = @(x,y) fitresult.p00 + fitresult.p10*x + fitresult.p01*y ...
        + fitresult.p20*x.^2 + fitresult.p11*x.*y + fitresult.p02*y.^2;
end
end

function plot_layer(transformedData, surface_fn, Faces, Vertices, Ref_mean, coeffs)
% Fitted surface + reference synapses + PCA-aligned neuropil mesh.
% Grid spacing is 500 voxels * 8 nm = 4 um (coordinates are in metres).
hold on;
spacing = 500*8e-9;
[X,Y] = meshgrid(min(transformedData(:,1)):spacing:max(transformedData(:,1)), ...
                 min(transformedData(:,2)):spacing:max(transformedData(:,2)));
surf(X, Y, surface_fn(X,Y), 'EdgeColor','none')
scatter3(transformedData(:,1), transformedData(:,2), transformedData(:,3), 20, [0.8500 0.3250 0.0980], 'filled');
grid on;

Vertices_PCA = (Vertices - Ref_mean) * coeffs;
trisurf(Faces, Vertices_PCA(:,1), Vertices_PCA(:,2), Vertices_PCA(:,3), 'FaceAlpha',0.2, 'EdgeColor','None');

axis('equal')
xlabel('PC1','Interpreter','none'); ylabel('PC2','Interpreter','none'); zlabel('PC3','Interpreter','none');
view([-31.459612469645986, 5.86855640905495])
end
