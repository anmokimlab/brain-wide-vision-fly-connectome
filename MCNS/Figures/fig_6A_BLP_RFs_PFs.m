%% ============================================
%  BLP neurons: weighted (theta, phi) projective/receptive fields  [MCNS]
%  MCNS analogue of the FAFB Figures/fig_6E_G_BLP_RFs_PFs.m.
%  - Tm3/T4a columns store theta/phi as N x 2: [value, weight], sum(weights)=1
%  - Mi1 columns keep a scalar (theta, phi) per column
%  - Each synapse maps to its nearest column; for Tm3/T4a the whole weighted
%    (theta, phi) distribution is added, for Mi1 a single (theta, phi) with weight 1.
%  - IN = orange (receptive field), OUT = blue (projective field); per-neuron heatmaps
%    and overlays plus type-level per-neuron contours.
%  Unlike the FAFB version (which assigns synapses to a neuropil via mesh
%  intriangulation), the MCNS synapse table carries a `neuropil` column, so synapses
%  are assigned to a region by that column.
%% ============================================

%% 1. Load data and initialize
clear all; close all; clc;

baseDir = fileparts(fileparts(mfilename('fullpath')));  % resolved from this script's location (portable; no absolute paths)

% Reference columns (Mi1 / Tm3 / T4a), from Data_Processing/s16_Mi1_Tm3_T4a_matching.m
load(fullfile(baseDir, 'Processed_Data', 'Mi1_Tm3_T4a_columns.mat'))   % Mi1_Columns, Tm3_Columns, T4a_Columns

% Reference-layer PCA bases, from Data_Processing/s17_making_layers.m
load(fullfile(baseDir, 'Processed_Data', 'neuropil_PCA_basis.mat'))    % PCA_Basis

% BLP neuron classification (provides BLP_R_NPIs, BLP_R_by_type)
load(fullfile(baseDir, 'Processed_Data', 'BLP_neurons_thr0.mat'))

% MCNS connectivity
opt = detectImportOptions(fullfile(baseDir, 'MCNS_Data', 'male-cns-v0.9-all-connections.csv'));
opt = setvartype(opt,'pre_root_id','int64');
opt = setvartype(opt,'post_root_id','int64');
MCNSConnections = readtable(fullfile(baseDir, 'MCNS_Data', 'male-cns-v0.9-all-connections.csv'),opt);

% Synapse coordinates (8 nm voxels -> metres; columns 1:3 = pre_x/pre_y/pre_z)
opt = detectImportOptions(fullfile(baseDir, 'MCNS_Data', 'male-cns-v0.9-synapse-coordinates.csv'));
opt = setvartype(opt,'pre_root_id','int64');
opt = setvartype(opt,'post_root_id','int64');
MCNS_synapse_coordinates = readtable(fullfile(baseDir, 'MCNS_Data', 'male-cns-v0.9-synapse-coordinates.csv'),opt);
MCNS_synapse_coordinates.pre_x = MCNS_synapse_coordinates.pre_x*8e-9;
MCNS_synapse_coordinates.pre_y = MCNS_synapse_coordinates.pre_y*8e-9;
MCNS_synapse_coordinates.pre_z = MCNS_synapse_coordinates.pre_z*8e-9;

%% Color setup
n = 256; startColor = [1,1,1];
endColorOut_orange = [0.8500, 0.3250, 0.0980]; % IN: orange
endColorIn_blue    = [0, 0.4470, 0.7410];      % OUT: blue
fillColor_in  = endColorOut_orange;
fillColor_out = endColorIn_blue;

%% Neuron types to analyze
neuron_types = unique(BLP_R_by_type.type);

% Prepare T4a columns (theta,phi as {i}=N×2 [value,weight])
T4a_R_Columns = T4a_Columns(strcmp(T4a_Columns.hemisphere,'R'),:);
T4a_L_Columns = T4a_Columns(strcmp(T4a_Columns.hemisphere,'L'),:);
T4a_R_out_syn_positions = cellfun(@table2array, T4a_R_Columns.out_syn_loc_denoise, 'UniformOutput', false); % OK
T4a_L_out_syn_positions = cellfun(@table2array, T4a_L_Columns.out_syn_loc_denoise, 'UniformOutput', false); % OK

% Prepare Tm3 columns (theta,phi as {i}=N×2 [value,weight])
Tm3_R_Columns = Tm3_Columns(strcmp(Tm3_Columns.hemisphere,'R'),:);
Tm3_L_Columns = Tm3_Columns(strcmp(Tm3_Columns.hemisphere,'L'),:);
Tm3_R_out_syn_positions = cellfun(@table2array, Tm3_R_Columns.out_syn_loc_denoise, 'UniformOutput', false); % OK
Tm3_L_out_syn_positions = cellfun(@table2array, Tm3_L_Columns.out_syn_loc_denoise, 'UniformOutput', false); % OK

% Prepare Mi1 columns (existing scalar form)
Mi1_R_Columns = Mi1_Columns(strcmp(Mi1_Columns.hemisphere,'R'),:);
Mi1_L_Columns = Mi1_Columns(strcmp(Mi1_Columns.hemisphere,'L'),:);
Mi1_R_out_syn_positions = cellfun(@table2array, Mi1_R_Columns.out_syn_loc_denoise, 'UniformOutput', false);
Mi1_L_out_syn_positions = cellfun(@table2array, Mi1_L_Columns.out_syn_loc_denoise, 'UniformOutput', false);
Mi1_theta = Mi1_Columns.theta;  % for background scatter
Mi1_phi   = Mi1_Columns.phi;    % for background scatter

%% Main loop
all_BLP_by_type = struct();  % stores current_neurons per type

for nt = 1:length(neuron_types)
    type_name = neuron_types{nt};
    fprintf('Processing: %s\n', type_name)

    this_type_idx   = strcmp(BLP_R_NPIs.type, type_name);
    current_neurons = table(BLP_R_NPIs.root_id(this_type_idx), 'VariableNames', {'root_id'});

    % --- Extract & aggregate synapse positions ---
    for i = 1:height(current_neurons)
        rid = current_neurons.root_id(i);

        InConnections_idx  = MCNSConnections.post_root_id==rid;
        InConnections      = MCNSConnections(InConnections_idx,:);
        OutConnections_idx = MCNSConnections.pre_root_id==rid;
        OutConnections     = MCNSConnections(OutConnections_idx,:);

        N_in_synapse_Me_R  = sum(InConnections.syn_count(strcmp(InConnections.neuropil,'ME(R)')));
        N_in_synapse_Me_L  = sum(InConnections.syn_count(strcmp(InConnections.neuropil,'ME(L)')));
        N_in_synapse_Lo_R  = sum(InConnections.syn_count(strcmp(InConnections.neuropil,'LO(R)')));
        N_in_synapse_Lo_L  = sum(InConnections.syn_count(strcmp(InConnections.neuropil,'LO(L)')));
        N_in_synapse_LoP_R = sum(InConnections.syn_count(strcmp(InConnections.neuropil,'LOP(R)')));
        N_in_synapse_LoP_L = sum(InConnections.syn_count(strcmp(InConnections.neuropil,'LOP(L)')));

        N_out_synapse_Me_R  = sum(OutConnections.syn_count(strcmp(OutConnections.neuropil,'ME(R)')));
        N_out_synapse_Me_L  = sum(OutConnections.syn_count(strcmp(OutConnections.neuropil,'ME(L)')));
        N_out_synapse_Lo_R  = sum(OutConnections.syn_count(strcmp(OutConnections.neuropil,'LO(R)')));
        N_out_synapse_Lo_L  = sum(OutConnections.syn_count(strcmp(OutConnections.neuropil,'LO(L)')));
        N_out_synapse_LoP_R = sum(OutConnections.syn_count(strcmp(OutConnections.neuropil,'LOP(R)')));
        N_out_synapse_LoP_L = sum(OutConnections.syn_count(strcmp(OutConnections.neuropil,'LOP(L)')));

        current_neurons.N_in_synapse_Me_R{i}  = N_in_synapse_Me_R;
        current_neurons.N_in_synapse_Me_L{i}  = N_in_synapse_Me_L;
        current_neurons.N_in_synapse_Lo_R{i}  = N_in_synapse_Lo_R;
        current_neurons.N_in_synapse_Lo_L{i}  = N_in_synapse_Lo_L;
        current_neurons.N_in_synapse_LoP_R{i} = N_in_synapse_LoP_R;
        current_neurons.N_in_synapse_LoP_L{i} = N_in_synapse_LoP_L;

        current_neurons.N_out_synapse_Me_R{i}  = N_out_synapse_Me_R;
        current_neurons.N_out_synapse_Me_L{i}  = N_out_synapse_Me_L;
        current_neurons.N_out_synapse_Lo_R{i}  = N_out_synapse_Lo_R;
        current_neurons.N_out_synapse_Lo_L{i}  = N_out_synapse_Lo_L;
        current_neurons.N_out_synapse_LoP_R{i} = N_out_synapse_LoP_R;
        current_neurons.N_out_synapse_LoP_L{i} = N_out_synapse_LoP_L;

        current_neurons.N_in_synpase_Optic{i}  = N_in_synapse_Me_R + N_in_synapse_Me_L + N_in_synapse_Lo_R + N_in_synapse_Lo_L + N_in_synapse_LoP_R + N_in_synapse_LoP_L;
        current_neurons.N_out_synpase_Optic{i} = N_out_synapse_Me_R + N_out_synapse_Me_L + N_out_synapse_Lo_R + N_out_synapse_Lo_L + N_out_synapse_LoP_R + N_out_synapse_LoP_L;
    end

    % --- Compute mean fractions ---
    in_Me_R   = cell2mat(current_neurons.N_in_synapse_Me_R);
    in_Me_L   = cell2mat(current_neurons.N_in_synapse_Me_L);
    in_Lo_R   = cell2mat(current_neurons.N_in_synapse_Lo_R);
    in_Lo_L   = cell2mat(current_neurons.N_in_synapse_Lo_L);
    in_LoP_R  = cell2mat(current_neurons.N_in_synapse_LoP_R);
    in_LoP_L  = cell2mat(current_neurons.N_in_synapse_LoP_L);

    out_Me_R  = cell2mat(current_neurons.N_out_synapse_Me_R);
    out_Me_L  = cell2mat(current_neurons.N_out_synapse_Me_L);
    out_Lo_R  = cell2mat(current_neurons.N_out_synapse_Lo_R);
    out_Lo_L  = cell2mat(current_neurons.N_out_synapse_Lo_L);
    out_LoP_R = cell2mat(current_neurons.N_out_synapse_LoP_R);
    out_LoP_L = cell2mat(current_neurons.N_out_synapse_LoP_L);

    Total_in_synases  = in_Me_R + in_Me_L + in_Lo_R + in_Lo_L + in_LoP_R + in_LoP_L;
    Total_out_synases = out_Me_R + out_Me_L + out_Lo_R + out_Lo_L + out_LoP_R + out_LoP_L;

    tin  = Total_in_synases;  tin(tin==0)   = NaN;
    tout = Total_out_synases; tout(tout==0) = NaN;

    Mean_N_in_synapse_Me_R   = mean(in_Me_R  ./ tin,  'omitnan');
    Mean_N_in_synapse_Me_L   = mean(in_Me_L  ./ tin,  'omitnan');
    Mean_N_in_synapse_Lo_R   = mean(in_Lo_R  ./ tin,  'omitnan');
    Mean_N_in_synapse_Lo_L   = mean(in_Lo_L  ./ tin,  'omitnan');
    Mean_N_in_synapse_LoP_R  = mean(in_LoP_R ./ tin,  'omitnan');
    Mean_N_in_synapse_LoP_L  = mean(in_LoP_L ./ tin,  'omitnan');

    Mean_N_out_synapse_Me_R  = mean(out_Me_R  ./ tout, 'omitnan');
    Mean_N_out_synapse_Me_L  = mean(out_Me_L  ./ tout, 'omitnan');
    Mean_N_out_synapse_Lo_R  = mean(out_Lo_R  ./ tout, 'omitnan');
    Mean_N_out_synapse_Lo_L  = mean(out_Lo_L  ./ tout, 'omitnan');
    Mean_N_out_synapse_LoP_R = mean(out_LoP_R ./ tout, 'omitnan');
    Mean_N_out_synapse_LoP_L = mean(out_LoP_L ./ tout, 'omitnan');

    % --- Filter and store synapse coordinates (by neuropil column) ---
    for i = 1:height(current_neurons)
        rid = current_neurons.root_id(i);

        out_syn_idx = find(ismember(MCNS_synapse_coordinates.pre_root_id, rid));
        out_locs    = MCNS_synapse_coordinates(out_syn_idx, :);
        in_syn_idx  = find(ismember(MCNS_synapse_coordinates.post_root_id, rid));
        in_locs     = MCNS_synapse_coordinates(in_syn_idx, :);

        current_neurons.in_syn_loc{i}  = in_locs{:,1:3};
        current_neurons.out_syn_loc{i} = out_locs{:,1:3};

        empty_loc_table = table([], [], [], 'VariableNames', {'x', 'y', 'z'});

        % Store input synapse positions (threshold 0.2)
        if Mean_N_in_synapse_Me_R >= 0.2
            in_locs_Me_R = in_locs(strcmpi(in_locs.neuropil,'ME(R)'),:);
            current_neurons.in_syn_loc_Me_R{i} = table(in_locs_Me_R{:,1},in_locs_Me_R{:,2},in_locs_Me_R{:,3},'VariableNames', {'x', 'y', 'z'});
        else
            current_neurons.in_syn_loc_Me_R{i} = empty_loc_table;
        end

        if Mean_N_in_synapse_Me_L >= 0.2
            in_locs_Me_L = in_locs(strcmpi(in_locs.neuropil,'ME(L)'),:);
            current_neurons.in_syn_loc_Me_L{i} = table(in_locs_Me_L{:,1},in_locs_Me_L{:,2},in_locs_Me_L{:,3},'VariableNames', {'x', 'y', 'z'});
        else
            current_neurons.in_syn_loc_Me_L{i} = empty_loc_table;
        end

        if Mean_N_in_synapse_Lo_R >= 0.2
            in_locs_Lo_R = in_locs(strcmpi(in_locs.neuropil,'LO(R)'),:);
            current_neurons.in_syn_loc_Lo_R{i} = table(in_locs_Lo_R{:,1},in_locs_Lo_R{:,2},in_locs_Lo_R{:,3},'VariableNames', {'x', 'y', 'z'});
        else
            current_neurons.in_syn_loc_Lo_R{i} = empty_loc_table;
        end

        if Mean_N_in_synapse_Lo_L >= 0.2
            in_locs_Lo_L = in_locs(strcmpi(in_locs.neuropil,'LO(L)'),:);
            current_neurons.in_syn_loc_Lo_L{i} = table(in_locs_Lo_L{:,1},in_locs_Lo_L{:,2},in_locs_Lo_L{:,3},'VariableNames', {'x', 'y', 'z'});
        else
            current_neurons.in_syn_loc_Lo_L{i} = empty_loc_table;
        end

        if Mean_N_in_synapse_LoP_R >= 0.2
            in_locs_LoP_R = in_locs(strcmpi(in_locs.neuropil,'LOP(R)'),:);
            current_neurons.in_syn_loc_LoP_R{i} = table(in_locs_LoP_R{:,1},in_locs_LoP_R{:,2},in_locs_LoP_R{:,3},'VariableNames', {'x', 'y', 'z'});
        else
            current_neurons.in_syn_loc_LoP_R{i} = empty_loc_table;
        end

        if Mean_N_in_synapse_LoP_L >= 0.2
            in_locs_LoP_L = in_locs(strcmpi(in_locs.neuropil,'LOP(L)'),:);
            current_neurons.in_syn_loc_LoP_L{i} = table(in_locs_LoP_L{:,1},in_locs_LoP_L{:,2},in_locs_LoP_L{:,3},'VariableNames', {'x', 'y', 'z'});
        else
            current_neurons.in_syn_loc_LoP_L{i} = empty_loc_table;
        end

        % Store output synapse positions
        if Mean_N_out_synapse_Me_R >= 0.2
            out_locs_Me_R = out_locs(strcmpi(out_locs.neuropil,'ME(R)'),:);
            current_neurons.out_syn_loc_Me_R{i} = table(out_locs_Me_R{:,1},out_locs_Me_R{:,2},out_locs_Me_R{:,3},'VariableNames', {'x', 'y', 'z'});
        else
            current_neurons.out_syn_loc_Me_R{i} = empty_loc_table;
        end

        if Mean_N_out_synapse_Me_L >= 0.2
            out_locs_Me_L = out_locs(strcmpi(out_locs.neuropil,'ME(L)'),:);
            current_neurons.out_syn_loc_Me_L{i} = table(out_locs_Me_L{:,1},out_locs_Me_L{:,2},out_locs_Me_L{:,3},'VariableNames', {'x', 'y', 'z'});
        else
            current_neurons.out_syn_loc_Me_L{i} = empty_loc_table;
        end

        if Mean_N_out_synapse_Lo_R >= 0.2
            out_locs_Lo_R = out_locs(strcmpi(out_locs.neuropil,'LO(R)'),:);
            current_neurons.out_syn_loc_Lo_R{i} = table(out_locs_Lo_R{:,1},out_locs_Lo_R{:,2},out_locs_Lo_R{:,3},'VariableNames', {'x', 'y', 'z'});
        else
            current_neurons.out_syn_loc_Lo_R{i} = empty_loc_table;
        end

        if Mean_N_out_synapse_Lo_L >= 0.2
            out_locs_Lo_L = out_locs(strcmpi(out_locs.neuropil,'LO(L)'),:);
            current_neurons.out_syn_loc_Lo_L{i} = table(out_locs_Lo_L{:,1},out_locs_Lo_L{:,2},out_locs_Lo_L{:,3},'VariableNames', {'x', 'y', 'z'});
        else
            current_neurons.out_syn_loc_Lo_L{i} = empty_loc_table;
        end

        if Mean_N_out_synapse_LoP_R >= 0.2
            out_locs_LoP_R = out_locs(strcmpi(out_locs.neuropil,'LOP(R)'),:);
            current_neurons.out_syn_loc_LoP_R{i} = table(out_locs_LoP_R{:,1},out_locs_LoP_R{:,2},out_locs_LoP_R{:,3},'VariableNames', {'x', 'y', 'z'});
        else
            current_neurons.out_syn_loc_LoP_R{i} = empty_loc_table;
        end

        if Mean_N_out_synapse_LoP_L >= 0.2
            out_locs_LoP_L = out_locs(strcmpi(out_locs.neuropil,'LOP(L)'),:);
            current_neurons.out_syn_loc_LoP_L{i} = table(out_locs_LoP_L{:,1},out_locs_LoP_L{:,2},out_locs_LoP_L{:,3},'VariableNames', {'x', 'y', 'z'});
        else
            current_neurons.out_syn_loc_LoP_L{i} = empty_loc_table;
        end
    end

    % Remove neurons with zero input/output
    idx_in_zero  = cell2mat(current_neurons.N_in_synpase_Optic)  == 0;
    idx_out_zero = cell2mat(current_neurons.N_out_synpase_Optic) == 0;
    current_neurons(idx_in_zero | idx_out_zero, :) = [];

    column_selection = 'all';

    % --- Map theta/phi (using weighted distribution) ---
    for i = 1:height(current_neurons)
        % ===== Input =====
        input_targets = {
            'Me_R',  Mi1_R_out_syn_positions, Mi1_R_Columns, 'unit';     % Mi1: existing approach (single)
            'Me_L',  Mi1_L_out_syn_positions, Mi1_L_Columns, 'unit';
            'Lo_R',  Tm3_R_out_syn_positions,  Tm3_R_Columns,  'weighted'; % Tm3: value+frequency
            'Lo_L',  Tm3_L_out_syn_positions,  Tm3_L_Columns,  'weighted';
            'LoP_R', T4a_R_out_syn_positions, T4a_R_Columns, 'weighted'; % T4a: value+frequency
            'LoP_L', T4a_L_out_syn_positions, T4a_L_Columns, 'weighted'};

        theta_in_all = []; phi_in_all = []; w_in_all = [];
        for k = 1:size(input_targets,1)
            region      = input_targets{k,1};
            ref_syns  = input_targets{k,2};
            ref_cols    = input_targets{k,3};
            mode_k      = input_targets{k,4};

            pos = table2array(current_neurons.(['in_syn_loc_' region]){i});
            if isempty(pos)
                theta = []; phi = []; w = []; root_id=[];
            else
                [theta, phi, w, root_id] = map_pos_to_weighted_angles(pos, ref_syns, ref_cols, mode_k, column_selection,k,PCA_Basis);
            end

            theta_in_all = [theta_in_all; theta];
            phi_in_all   = [phi_in_all;   phi];
            w_in_all     = [w_in_all;     w];

            % (optional) recording root id is skipped
            current_neurons.([region '_in']){i} = [root_id];
        end
        current_neurons.in_syn_theta{i} = theta_in_all;
        current_neurons.in_syn_phi{i}   = phi_in_all;
        current_neurons.in_weights{i}   = w_in_all;   % store weights

        % ===== Output =====
        output_targets = {
            'Me_R',  Mi1_R_out_syn_positions, Mi1_R_Columns, 'unit';
            'Me_L',  Mi1_L_out_syn_positions, Mi1_L_Columns, 'unit';
            'Lo_R',  Tm3_R_out_syn_positions,  Tm3_R_Columns,  'weighted';
            'Lo_L',  Tm3_L_out_syn_positions,  Tm3_L_Columns,  'weighted';
            'LoP_R', T4a_R_out_syn_positions, T4a_R_Columns, 'weighted';
            'LoP_L', T4a_L_out_syn_positions, T4a_L_Columns, 'weighted'};

        theta_out_all = []; phi_out_all = []; w_out_all = [];
        for k = 1:size(output_targets,1)
            region      = output_targets{k,1};
            ref_syns  = output_targets{k,2};
            ref_cols    = output_targets{k,3};
            mode_k      = output_targets{k,4};

            pos = table2array(current_neurons.(['out_syn_loc_' region]){i});
            if isempty(pos)
                theta = []; phi = []; w = []; root_id=[];
            else
                [theta, phi, w, root_id] = map_pos_to_weighted_angles(pos, ref_syns, ref_cols, mode_k,column_selection,k,PCA_Basis);
            end
            theta_out_all = [theta_out_all; theta];
            phi_out_all   = [phi_out_all;   phi];
            w_out_all     = [w_out_all;     w];

            current_neurons.([region '_out']){i} = [root_id];
        end
        current_neurons.out_syn_theta{i} = theta_out_all;
        current_neurons.out_syn_phi{i}   = phi_out_all;
        current_neurons.out_weights{i}   = w_out_all;   % store weights
    end

    %% ===== Per-neuron: IN(orange) + OUT(blue) heatmap overlay =====
    % - Accumulate Tin/Tout (phi, theta, count) points into a 2D weighted histogram
    % - Overlay two axes on the same plot, visualized with different colormaps and alpha
    % - Heatmap handles empty/sparse regions naturally (transparent where weight=0)

    % Heatmap resolution (spacing, unit: degrees)
    dphi   = 2;                       % horizontal (φ) spacing
    dtheta = 2;                       % vertical (θ) spacing
    phi_edges   = -180:dphi:180;      phi_centers   = (phi_edges(1:end-1) + phi_edges(2:end))/2;
    theta_edges =  -90:dtheta: 90;    theta_centers = (theta_edges(1:end-1) + theta_edges(2:end))/2;

    smoothing_on = true;     % set false to disable smoothing
    sigma_phi_deg   = 9.5/2.355;     % horizontal (φ) smoothing amount
    sigma_theta_deg = 9.5/2.355;     % vertical (θ) smoothing amount

    draw_contours   = true;
    contour_levels  = [0.1 0.1];

    % colormap (white→orange, white→blue)
    nC = 256;
    cmap_orange = [linspace(1,0.85,nC)', linspace(1,0.325,nC)', linspace(1,0.098,nC)'];
    cmap_blue   = [linspace(1,0.00,nC)', linspace(1,0.447,nC)', linspace(1,0.741,nC)'];

    % gamma (emphasize/soften high frequencies)
    gamma_in  = 1.0;
    gamma_out = 1.0;

    for i = 1:height(current_neurons)

        % ============================================================
        % (added) build in_tiles / out_tiles if missing (weighted sum)
        % ============================================================
        % Build IN tiles
        need_make_in = ~ismember('in_tiles', current_neurons.Properties.VariableNames) ...
            || isempty(current_neurons.in_tiles{i});
        if need_make_in
            phi_in   = current_neurons.in_syn_phi{i};
            theta_in = current_neurons.in_syn_theta{i};
            w_in     = current_neurons.in_weights{i};

            % convert to numeric vector if cell
            if iscell(phi_in),   phi = cell2mat(phi_in(:));   else, phi = phi_in(:);   end
            if iscell(theta_in), the = cell2mat(theta_in(:)); else, the = theta_in(:); end
            % tidy up weights
            if isempty(w_in)
                w = ones(min(numel(phi), numel(the)), 1);
            else
                if iscell(w_in), w = cell2mat(w_in(:)); else, w = w_in(:); end
            end

            % match lengths
            n = min([numel(phi), numel(the), numel(w)]);
            phi = phi(1:n); the = the(1:n); w = w(1:n);

            if isempty(phi)
                Tin_build = array2table(zeros(0,3),'VariableNames',{'phi','theta','count'});
            else
                pts = [phi(:) the(:)];
                [uniq, ~, ic] = unique(pts, 'rows', 'stable');
                counts = accumarray(ic, w, [], @sum);
                Tin_build = array2table([uniq counts], 'VariableNames', {'phi','theta','count'});
            end
            current_neurons.in_tiles{i} = Tin_build;
        end

        % (optional) build OUT tiles if also missing
        need_make_out = ~ismember('out_tiles', current_neurons.Properties.VariableNames) ...
            || isempty(current_neurons.out_tiles{i});
        if need_make_out
            phi_out   = current_neurons.out_syn_phi{i};
            theta_out = current_neurons.out_syn_theta{i};
            w_out     = current_neurons.out_weights{i};

            if iscell(phi_out),   phi = cell2mat(phi_out(:));   else, phi = phi_out(:);   end
            if iscell(theta_out), the = cell2mat(theta_out(:)); else, the = theta_out(:); end
            if isempty(w_out)
                w = ones(min(numel(phi), numel(the)), 1);
            else
                if iscell(w_out), w = cell2mat(w_out(:)); else, w = w_out(:); end
            end

            n = min([numel(phi), numel(the), numel(w)]);
            phi = phi(1:n); the = the(1:n); w = w(1:n);

            if isempty(phi)
                Tout_build = array2table(zeros(0,3),'VariableNames',{'phi','theta','count'});
            else
                pts = [phi(:) the(:)];
                [uniq, ~, ic] = unique(pts, 'rows', 'stable');
                counts = accumarray(ic, w, [], @sum);
                Tout_build = array2table([uniq counts], 'VariableNames', {'phi','theta','count'});
            end
            current_neurons.out_tiles{i} = Tout_build;
        end
        % ============================================================

        % --- Prepare data ---
        Tin  = current_neurons.in_tiles{i};   % columns: {'phi','theta','count'}
        Tout = current_neurons.out_tiles{i};

        % guard
        if isempty(Tin);  Tin = array2table(zeros(0,3),'VariableNames',{'phi','theta','count'}); end
        if isempty(Tout); Tout = array2table(zeros(0,3),'VariableNames',{'phi','theta','count'}); end

        % --- Weighted 2D histogram
        H_in  = local_weighted_hist2(Tin.theta,  Tin.phi,  Tin.count,  theta_edges, phi_edges);
        H_out = local_weighted_hist2(Tout.theta, Tout.phi, Tout.count, theta_edges, phi_edges);

        % ================== Gaussian smoothing option ==================


        if smoothing_on
            % bin spacing (deg) -> sigma in bin units


            % φ: circular padding, θ: edge-fixed
            H_in = local_gauss_smooth2_circ(H_in, [], [], struct( ...
                'mode','deg', ...
                'sigma_theta_deg', 9.5/2.355, ...  % based on Δρ=9.5°
                'sigma_phi_deg',   9.5/2.355, ...
                'bin_theta_deg',   2, ...
                'bin_phi_deg',     2));
            H_out = local_gauss_smooth2_circ(H_out, [], [], struct( ...
                'mode','deg', ...
                'sigma_theta_deg', 9.5/2.355, ...  % based on Δρ=9.5°
                'sigma_phi_deg',   9.5/2.355, ...
                'bin_theta_deg',   2, ...
                'bin_phi_deg',     2));

        end
        % ================================================================

        % Normalize (each independently; for a common scale use max_in_out = max([H_in(:);H_out(:)]))
        max_in  = max(H_in(:));  if max_in<=0,  max_in = 1; end
        max_out = max(H_out(:)); if max_out<=0, max_out = 1; end
        Hin  = (H_in  / max_in ).^gamma_in;
        Hout = (H_out / max_out).^gamma_out;

        % Clip very small values so they render transparent (optional)
        eps_clip = 1e-6;
        Hin(Hin < eps_clip)   = 0;
        Hout(Hout < eps_clip) = 0;

        % ========== [important] store results reused at the type level below ==========
        current_neurons.Hin_norm{i}  = Hin;     % after normalization + gamma + clipping
        current_neurons.Hout_norm{i} = Hout;
        % ===============================================================

        % --- Figure & two stacked axes (overlay) ---
        fig = figure('Name', sprintf('IN/OUT heatmap overlay: %s #%d', type_name, i)); clf
        set(fig,'Color','w');

        % Bottom axes: IN heatmap
        ax1 = axes; hold(ax1,'on');
        % Background Mi1 scatter
        scatter(ax1, Mi1_phi, Mi1_theta, 6, [0.85 0.85 0.85], 'filled', 'MarkerFaceAlpha', 0.20);
        % IN heatmap
        hi = imagesc(ax1, 'XData', [phi_edges(1) phi_edges(end)], ...
            'YData', [theta_edges(1) theta_edges(end)], ...
            'CData', Hin);
        set(ax1,'YDir','normal'); axis(ax1,'equal'); xlim(ax1,[-180 180]); ylim(ax1,[-90 90]);
        colormap(ax1, cmap_orange);
        % IN transparency: transparent where value is 0
        set(hi, 'AlphaData', Hin, 'AlphaDataMapping','none');

        % Top axes: OUT heatmap (transparent background)
        ax2 = axes('Position', get(ax1,'Position'), 'Color','none'); hold(ax2,'on');
        ho = imagesc(ax2, 'XData', [phi_edges(1) phi_edges(end)], ...
            'YData', [theta_edges(1) theta_edges(end)], ...
            'CData', Hout);
        set(ax2,'YDir','normal'); axis(ax2,'equal'); xlim(ax2,[-180 180]); ylim(ax2,[-90 90]);
        colormap(ax2, cmap_blue);
        % OUT transparency
        set(ho, 'AlphaData', Hout, 'AlphaDataMapping','none');

        % ================== IN/OUT contours (connect specific value) ==================


        if draw_contours
            col_in_1  = cmap_orange(end,:);
            col_out_1 = cmap_blue(end,:);

            % Draw contours on ax2 (the top axes)!
            if any(Hin(:)  >= min(contour_levels))
                contour(ax2, phi_centers, theta_centers, Hin,  contour_levels, ...
                    'LineColor', col_in_1,  'LineWidth', 1.8, 'LineStyle','-','Clipping','off');
            end
            if any(Hout(:) >= min(contour_levels))
                contour(ax2, phi_centers, theta_centers, Hout, contour_levels, ...
                    'LineColor', col_out_1, 'LineWidth', 1.8, 'LineStyle','-','Clipping','off');
            end
        end
        % ======================================================================

        % Grid/labels/title/axes setup (sync the two axes)
        linkaxes([ax1,ax2],'xy');
        set([ax1,ax2],'XTick',-180:45:180,'YTick',-90:45:90,'TickDir','out');
        xlabel(ax1,'\phi'); ylabel(ax1,'\theta');
        title(ax1, sprintf('IN (orange) + OUT (blue) heatmaps: %s (root id: %d)', ...
            type_name, current_neurons.root_id(i)));
        set(ax2,'XTick',[],'YTick',[]);  % hide top-axes ticks
    end

    %% ===== [type level] Draw all per-neuron contours at once (unique color per neuron) =====
    if height(current_neurons) > 0
        % Contour levels (relative to normalized heatmap) — same as above

        % (reuse the bin definitions from the per-neuron section)
        % dphi, dtheta, phi_edges, theta_edges, phi_centers, theta_centers are already defined above

        % Per-neuron colors (distinguishable_colors; fall back to lines if unavailable)
        try
            cmap = distinguishable_colors(height(current_neurons));
        catch
            warning('distinguishable_colors not on path; using default lines colors.');
            cmap = lines(height(current_neurons));
        end

        fig = figure('Name', sprintf('All neuron contours (per-neuron colors): %s', type_name));
        set(fig,'Color','w');
        ax = axes; hold(ax,'on');
        scatter(ax, Mi1_phi, Mi1_theta, 6, [0.85 0.85 0.85], 'filled', 'MarkerFaceAlpha', 0.20);

        % (optional) background ticks/axes
        axis(ax,'equal'); xlim(ax,[-180 180]); ylim(ax,[-90 90]);
        set(ax,'XTick',-180:45:180,'YTick',-90:45:90,'TickDir','out');
        xlabel(ax,'\phi'); ylabel(ax,'\theta');
        title(ax, sprintf('IN (solid) / OUT (dashed) contours by neuron: %s', type_name));

        for ii = 1:height(current_neurons)
            % ---------- reuse the (final) heatmap saved in the section above ----------
            Hin = []; Hout = [];
            if ismember('Hin_norm', current_neurons.Properties.VariableNames)
                Hin = current_neurons.Hin_norm{ii};
            end
            if ismember('Hout_norm', current_neurons.Properties.VariableNames)
                Hout = current_neurons.Hout_norm{ii};
            end

            % IN contour (solid)
            if ~isempty(Hin) && any(Hin(:) >= min(contour_levels))
                contour(ax, phi_centers, theta_centers, Hin,  contour_levels, ...
                    'LineColor', cmap(ii,:), 'LineStyle','-',  'LineWidth', 1.6, ...
                    'HandleVisibility','off', 'Clipping','off');
            end
            % OUT contour (dashed)
            if ~isempty(Hout) && any(Hout(:) >= min(contour_levels))
                contour(ax, phi_centers, theta_centers, Hout, contour_levels, ...
                    'LineColor', cmap(ii,:), 'LineStyle','--', 'LineWidth', 1.6, ...
                    'HandleVisibility','off', 'Clipping','off');
            end
        end

        title('Neuron (IN solid / OUT dashed)');
    end


    close all;
    safe_type_name = matlab.lang.makeValidName(type_name);
    all_BLP_by_type.(safe_type_name) = current_neurons;

end

% Save results
save(fullfile(baseDir, 'Processed_Data', 'BLP_RFs_PFs.mat'), 'all_BLP_by_type', '-v7.3');

%% ---------- Local functions ----------
function [theta, phi, w, src_root_id] = map_pos_to_weighted_angles(pos_array, ref_syn, ref_columns, mode, column_selection,k,PCA_Basis)
% mode: 'weighted' (Tm3/T4a) | 'unit' (Mi1)
% column_selection: 'all' | 'unique'
%   - 'all'    : keep every column each synapse matches, including duplicates
%   - 'unique' : keep a column only once even if selected multiple times
%
% Returns:
%   theta, phi : concatenated angle values (column vectors)
%   w          : per-angle weights (column vector)
%   src_root_id: identifier of which reference column each row came from (column vector)
%                priority: ref_columns.root_id > id > column_id > (else the column index)
%%%% keys
%ME__Dm6__dendrite

%LO__LT1a__dendrite

% LOP__T5a__axon

if k==1
    key='ME__Dm6__dendrite';
    PCA_Mean=PCA_Basis.(key).R.mean;
    PCA_coeffs=PCA_Basis.(key).R.coeffs;

    p00=PCA_Basis.(key).fit.R.p00;
    p10=PCA_Basis.(key).fit.R.p10;
    p01=PCA_Basis.(key).fit.R.p01;
    p20=PCA_Basis.(key).fit.R.p20;
    p11=PCA_Basis.(key).fit.R.p11;
    p02=PCA_Basis.(key).fit.R.p02;
    P=[p00 p10 p01 p20 p11 p02];

elseif k==2
    key='ME__Dm6__dendrite';
    PCA_Mean=PCA_Basis.(key).L.mean;
    PCA_coeffs=PCA_Basis.(key).L.coeffs;

    p00=PCA_Basis.(key).fit.L.p00;
    p10=PCA_Basis.(key).fit.L.p10;
    p01=PCA_Basis.(key).fit.L.p01;
    p20=PCA_Basis.(key).fit.L.p20;
    p11=PCA_Basis.(key).fit.L.p11;
    p02=PCA_Basis.(key).fit.L.p02;
    P=[p00 p10 p01 p20 p11 p02];

elseif k==3
    key='LO__LT1a__dendrite';
    PCA_Mean=PCA_Basis.(key).R.mean;
    PCA_coeffs=PCA_Basis.(key).R.coeffs;

    p00=PCA_Basis.(key).fit.R.p00;
    p10=PCA_Basis.(key).fit.R.p10;
    p01=PCA_Basis.(key).fit.R.p01;
    p20=PCA_Basis.(key).fit.R.p20;
    p11=PCA_Basis.(key).fit.R.p11;
    p02=PCA_Basis.(key).fit.R.p02;
    P=[p00 p10 p01 p20 p11 p02];

elseif k==4
    key='LO__LT1a__dendrite';
    PCA_Mean=PCA_Basis.(key).L.mean;
    PCA_coeffs=PCA_Basis.(key).L.coeffs;

    p00=PCA_Basis.(key).fit.L.p00;
    p10=PCA_Basis.(key).fit.L.p10;
    p01=PCA_Basis.(key).fit.L.p01;
    p20=PCA_Basis.(key).fit.L.p20;
    p11=PCA_Basis.(key).fit.L.p11;
    p02=PCA_Basis.(key).fit.L.p02;
    P=[p00 p10 p01 p20 p11 p02];

elseif k==5
    key='LOP__T5a__axon';
    PCA_Mean=PCA_Basis.(key).R.mean;
    PCA_coeffs=PCA_Basis.(key).R.coeffs;

    p00=PCA_Basis.(key).fit.R.p00;
    p10=PCA_Basis.(key).fit.R.p10;
    p01=PCA_Basis.(key).fit.R.p01;
    p20=PCA_Basis.(key).fit.R.p20;
    p11=PCA_Basis.(key).fit.R.p11;
    p02=PCA_Basis.(key).fit.R.p02;
    p30=PCA_Basis.(key).fit.R.p30;
    p21=PCA_Basis.(key).fit.R.p21;
    p12=PCA_Basis.(key).fit.R.p12;
    p03=PCA_Basis.(key).fit.R.p03;

    P=[p00 p10 p01 p20 p11 p02 p30 p21 p12 p03];


elseif k==6
    key='LOP__T5a__axon';
    PCA_Mean=PCA_Basis.(key).L.mean;
    PCA_coeffs=PCA_Basis.(key).L.coeffs;

    p00=PCA_Basis.(key).fit.L.p00;
    p10=PCA_Basis.(key).fit.L.p10;
    p01=PCA_Basis.(key).fit.L.p01;
    p20=PCA_Basis.(key).fit.L.p20;
    p11=PCA_Basis.(key).fit.L.p11;
    p02=PCA_Basis.(key).fit.L.p02;
    p30=PCA_Basis.(key).fit.L.p30;
    p21=PCA_Basis.(key).fit.L.p21;
    p12=PCA_Basis.(key).fit.L.p12;
    p03=PCA_Basis.(key).fit.L.p03;

    P=[p00 p10 p01 p20 p11 p02 p30 p21 p12 p03];
else
    msg = 'Check neuropil';
    error(msg)
end



%% ---- (NEW) Step 1: find the "nearest ref column" index by mean on-surface distance ----
% Assumes key/PCA_Mean/PCA_coeffs/P were already selected above from the k value
% Distinguish Poly22 (true closest point) vs Poly33 (vertical drop) by length of P
is_poly22 = numel(P) == 6;
is_poly33 = numel(P) == 10;
if ~(is_poly22 || is_poly33)
    error('P must be length 6 (poly22) or 10 (poly33).');
end

% poly33 z evaluation function
poly33_eval = @(PP, x, y) ( ...
    PP(1) + PP(2)*x + PP(3)*y + PP(4)*x.^2 + PP(5)*x.*y + PP(6)*y.^2 + ...
    PP(7)*x.^3 + PP(8)*x.^2.*y + PP(9)*x.*y.^2 + PP(10)*y.^3 );

% Numerical options (could be exposed as an opts argument if desired)
cp_opts.tol   = 1e-12;
cp_opts.maxit = 100;

% 0) Precompute on-surface points for each ref_syn neuron (in PCA coordinates)
n_ref = numel(ref_syn);
ref_onSurf = cell(n_ref,1);     % each cell: (m_i x 3) on-surface points
for i = 1:n_ref
    Xi = ref_syn{i};
    if isempty(Xi), ref_onSurf{i} = []; continue; end
    Xi = double(Xi);
    % into PCA coordinates
    Xi_pca = (Xi - PCA_Mean) * PCA_coeffs;  % (m_i x 3)
    if is_poly22
        % project each point to its true closest point
        Ri = zeros(size(Xi_pca));
        for j = 1:size(Xi_pca,1)
            [xs, ys, zs] = closest_point_poly22(P, Xi_pca(j,:), cp_opts);
            Ri(j,:) = [xs, ys, zs];
        end
    else
        % poly33: vertical drop on z only via f33(x,y)
        x = Xi_pca(:,1); y = Xi_pca(:,2);
        z = poly33_eval(P, x, y);
        Ri = [x, y, z];
    end
    ref_onSurf{i} = Ri;
end

% 1) Project each point of pos_array via PCA -> on-surface point
N = size(pos_array, 1);
nearest_idx = zeros(N,1);
for j = 1:N
    pj = double(pos_array(j,:));
    pj_pca = (pj - PCA_Mean) * PCA_coeffs;   % 1x3
    if is_poly22
        [qx, qy, qz] = closest_point_poly22(P, pj_pca, cp_opts);
        q = [qx, qy, qz];
    else
        % poly33: vertical drop on z only via f33(x,y)
        qx = pj_pca(1); qy = pj_pca(2);
        qz = poly33_eval(P, qx, qy);
        q  = [qx, qy, qz];
    end

    % 2) Compute "mean on-surface 3D distance" to each ref neuron i -> pick the minimum i
    best_i = 1; best_d = inf;
    for i = 1:n_ref
        Ri = ref_onSurf{i};
        if isempty(Ri), continue; end
        d = mean( sqrt(sum( (Ri - q).^2, 2 )) );  % mean 3D Euclidean distance in PCA coordinates
        if d < best_d
            best_d = d;
            best_i = i;
        end
    end
    nearest_idx(j) = best_i;
end

%% ---- (UNCHANGED) Step 2: derive idx_list according to column_selection ----
switch lower(column_selection)
    case 'unique'
        idx_list = unique(nearest_idx, 'stable');   % preserve order of appearance
    otherwise % 'all'
        idx_list = nearest_idx;
end

theta = []; phi = []; w = []; src_root_id = int64([]);


% (existing code retained) decide the id field name up front
vars     = ref_columns.Properties.VariableNames;


% (existing code retained) avoid loop-variable name clash when iterating idx_list
for ii = 1:numel(idx_list)
    min_idx = idx_list(ii);

    % Determine the reference ID
    rid = ref_columns.root_id(min_idx);


    switch mode
        case 'weighted'  % (existing) theta{idx}, phi{idx} = N×2 [value, weight]
            theta_mat = ref_columns.theta{min_idx};
            phi_mat   = ref_columns.phi{min_idx};
            if isempty(theta_mat) || isempty(phi_mat), continue; end

            n_ = min(size(theta_mat,1), size(phi_mat,1));
            tval = theta_mat(1:n_,1);
            pval = phi_mat(1:n_,1);
            wt   = theta_mat(1:n_,2);
            if any(~isfinite(wt)) || any(wt<0) || sum(wt)<=0
                wt = ones(n_,1)/n_;
            else
                wt = wt / sum(wt);
            end
            theta = [theta; tval];
            phi   = [phi;   pval];
            w     = [w;     wt];
            src_root_id = [src_root_id; repmat(rid, n_, 1)];

        case 'unit'      % (existing) scalar coordinate, weight=1
            tval = ref_columns.theta(min_idx);
            pval = ref_columns.phi(min_idx);
            theta = [theta; tval];
            phi   = [phi;   pval];
            w     = [w;     1];
            src_root_id = [src_root_id; rid];
    end
end

end

function H = local_weighted_hist2(theta_vals, phi_vals, weights, theta_edges, phi_edges)
% local_weighted_hist2
% - For MATLAB versions where histcounts2 does not support 'Weights'
% - Keeps the convention X-axis=theta (rows), Y-axis=phi (cols)
%
% Inputs: theta_vals, phi_vals, weights : all column vectors (double)
%         theta_edges, phi_edges        : bin edge vectors
% Output: H (numel(theta_edges)-1) x (numel(phi_edges)-1) weighted-sum matrix

% Safety: convert table-like inputs to arrays
if istable(theta_vals), theta_vals = theta_vals{:,:}; end
if istable(phi_vals),   phi_vals   = phi_vals{:,:};   end
if istable(weights),    weights    = weights{:,:};    end

theta_vals = theta_vals(:);
phi_vals   = phi_vals(:);
if isempty(weights)
    weights = ones(min(numel(theta_vals), numel(phi_vals)),1);
else
    weights = weights(:);
end

% match lengths
n = min([numel(theta_vals), numel(phi_vals), numel(weights)]);
theta_vals = theta_vals(1:n);
phi_vals   = phi_vals(1:n);
weights    = weights(1:n);

% keep only valid values (drop NaN/Inf)
valid = isfinite(theta_vals) & isfinite(phi_vals) & isfinite(weights);
theta_vals = theta_vals(valid);
phi_vals   = phi_vals(valid);
weights    = weights(valid);

% bin indices
[~,~,ib_theta] = histcounts(theta_vals, theta_edges); % 0 means outside
[~,~,ib_phi]   = histcounts(phi_vals,   phi_edges);

good = (ib_theta>0) & (ib_phi>0);
ib_theta = ib_theta(good);
ib_phi   = ib_phi(good);
weights  = weights(good);

nRow = numel(theta_edges)-1; % theta bins -> rows
nCol = numel(phi_edges)-1;   % phi bins   -> cols

% weighted sum
H = accumarray([ib_theta, ib_phi], weights, [nRow, nCol], @sum, 0);
end

function Hs = local_gauss_smooth2_circ(H, sigma_theta_bins, sigma_phi_bins, opts)
% Gaussian smoothing (separable kernel). φ is circular, θ is replicate.
% Compatible with the original signature; opts lets you flexibly specify units/parameters.
%
% Inputs:
%   H                 : (nTheta x nPhi) histogram (θ x φ)
%   sigma_theta_bins  : θ-direction σ (in bins)  [may be ignored depending on mode]
%   sigma_phi_bins    : φ-direction σ (in bins)  [may be ignored depending on mode]
%   opts (struct, optional) :
%       .mode            : 'bins' | 'deg' | 'accept' (default='bins')
%                          'bins'   : use sigma_*_bins directly
%                          'deg'    : convert to bins from sigma_*_deg and bin_*_deg
%                          'accept' : compute σ from acceptance_deg (=FWHM)
%       .sigma_theta_deg : θ-direction σ (deg, when mode='deg')
%       .sigma_phi_deg   : φ-direction σ (deg, when mode='deg')
%       .acceptance_deg  : Δρ = FWHM (deg, when mode='accept'; isotropic assumption)
%       .bin_theta_deg   : θ bin size (deg, required for mode='deg'/'accept')
%       .bin_phi_deg     : φ bin size (deg, required for mode='deg'/'accept')
%       .truncate        : kernel radius factor k (default=3; radius=ceil(k*σ))
%
% Output:
%   Hs : smoothed histogram

    if nargin < 2 || isempty(sigma_theta_bins), sigma_theta_bins = 1; end
    if nargin < 3 || isempty(sigma_phi_bins),   sigma_phi_bins   = 1; end
    if nargin < 4 || isempty(opts), opts = struct(); end
    if ~isfield(opts,'mode') || isempty(opts.mode), opts.mode = 'bins'; end
    if ~isfield(opts,'truncate') || isempty(opts.truncate), opts.truncate = 3; end

    [nRow, nCol] = size(H);

    % --- units/mode handling ---
    switch lower(opts.mode)
        case 'bins'
            sig_th_bins = max(eps, sigma_theta_bins);
            sig_ph_bins = max(eps, sigma_phi_bins);

        case 'deg'
            assert(isfield(opts,'sigma_theta_deg') && isfield(opts,'sigma_phi_deg') && ...
                   isfield(opts,'bin_theta_deg')   && isfield(opts,'bin_phi_deg'), ...
                   'mode="deg" requires sigma_*_deg and bin_*_deg.');
            sig_th_bins = max(eps, opts.sigma_theta_deg / opts.bin_theta_deg);
            sig_ph_bins = max(eps, opts.sigma_phi_deg   / opts.bin_phi_deg);

        case 'accept'
            % Δρ(=FWHM) -> σ(deg) -> σ(bins)
            assert(isfield(opts,'acceptance_deg') && isfield(opts,'bin_theta_deg') && isfield(opts,'bin_phi_deg'), ...
                   'mode="accept" requires acceptance_deg and bin_*_deg.');
            sigma_deg = opts.acceptance_deg / 2.355;   % σ = FWHM/2.355
            sig_th_bins = max(eps, sigma_deg / opts.bin_theta_deg);
            sig_ph_bins = max(eps, sigma_deg / opts.bin_phi_deg);

        otherwise
            error('Unknown mode: %s', opts.mode);
    end

    % If σ is very small (≈0 bins) return the original
    if sig_th_bins < 1e-6 && sig_ph_bins < 1e-6
        Hs = H; return;
    end

    % --- build kernel (normalized: energy-preserving) ---
    r_th = max(1, ceil(opts.truncate * sig_th_bins));
    r_ph = max(1, ceil(opts.truncate * sig_ph_bins));
    th   = -r_th:r_th;
    ph   = -r_ph:r_ph;

    g_th = exp(-(th.^2) / (2*sig_th_bins^2)); g_th = g_th / sum(g_th);
    g_ph = exp(-(ph.^2) / (2*sig_ph_bins^2)); g_ph = g_ph / sum(g_ph);

    % --- θ padding: replicate ---
    topPad    = repmat(H(1,:),  r_th, 1);
    bottomPad = repmat(H(end,:), r_th, 1);
    H_pad_th  = [topPad; H; bottomPad];

    % --- φ padding: circular ---
    leftPad   = H_pad_th(:, nCol-r_ph+1:nCol);
    rightPad  = H_pad_th(:, 1:r_ph);
    H_pad     = [leftPad, H_pad_th, rightPad];

    % --- separable convolution (vertical -> horizontal) ---
    H_tmp = conv2(g_th.', 1, H_pad, 'same');   % θ direction
    H_sm  = conv2(1, g_ph, H_tmp, 'same');     % φ direction

    % --- crop ---
    Hs = H_sm(r_th+1:r_th+nRow, r_ph+1:r_ph+nCol);
end

function [xstar, ystar, zstar, info] = closest_point_poly22(P, p, opts)
    if nargin<3, opts=struct; end
    if ~isfield(opts,'tol'),   opts.tol = 1e-16; end
    if ~isfield(opts,'maxit'), opts.maxit = 100;   end
    p00=P(1); p10=P(2); p01=P(3); p20=P(4); p11=P(5); p02=P(6);
    px=p(1); py=p(2); pz=p(3);
    f  = @(x,y) p00 + p10*x + p01*y + p20*x.^2 + p11*x.*y + p02*y.^2;
    fx = @(x,y) p10 + 2*p20*x + p11*y;
    fy = @(x,y) p01 + p11*x + 2*p02*y;
    fxx=2*p20; fxy=p11; fyy=2*p02;

    x=px; y=py;   % init: vertical drop
    for k=1:opts.maxit
        F=f(x,y); gx=fx(x,y); gy=fy(x,y); dz=F-pz;
        r1 = x - px + dz*gx;  r2 = y - py + dz*gy;
        if hypot(r1,r2) < opts.tol*(1+hypot(px,py)), break; end
        J11 = 1 + gx*gx + dz*fxx;   J12 = gx*gy + dz*fxy;
        J22 = 1 + gy*gy + dz*fyy;   J21 = J12;
        delta = -[J11 J12; J21 J22] \ [r1; r2];

        % simple backtracking
        step=1; old=hypot(r1,r2);
        for bt=1:6
            xn=x+step*delta(1); yn=y+step*delta(2);
            Fn=f(xn,yn); gxn=fx(xn,yn); gyn=fy(xn,yn); dzn=Fn-pz;
            new=hypot(xn-px + dzn*gxn, yn-py + dzn*gyn);
            if new <= 0.9*old || step < 1/64, x=xn; y=yn; break; end
            step=step*0.5;
        end
    end
    xstar=x; ystar=y; zstar=f(x,y);
    info.iter=k; info.dist=norm([x-px, y-py, f(x,y)-pz]);
end

function cmap = distinguishable_colors(n, bg)
    % Default background color: white
    if nargin < 2
        bg = [1, 1, 1];  % RGB
    end

    % Generate candidate RGB colors
    steps = 30;
    [R, G, B] = ndgrid(linspace(0, 1, steps));
    RGB = [R(:) G(:) B(:)];

    % Drop colors too close to the background
    lab_bg = rgb2lab(bg);
    lab_RGB = rgb2lab(RGB);
    dE_bg = sqrt(sum((lab_RGB - lab_bg).^2, 2));
    RGB = RGB(dE_bg > 10, :); % keep only colors with color difference >= 10

    % Pick the first color at random
    cmap = RGB(randi(size(RGB,1)), :);
    RGB_used = false(size(RGB,1),1);

    for i = 2:n
        % Distance to already-chosen colors (in Lab space)
        lab_cmap = rgb2lab(cmap);
        lab_all = rgb2lab(RGB);
        dists = pdist2(lab_all, lab_cmap);
        min_dist = min(dists, [], 2); % distance to the nearest chosen color
        min_dist(RGB_used) = -inf;    % exclude already-chosen colors

        [~, idx] = max(min_dist);
        cmap(end+1, :) = RGB(idx, :); %#ok<AGROW>
        RGB_used(idx) = true;
    end
end
