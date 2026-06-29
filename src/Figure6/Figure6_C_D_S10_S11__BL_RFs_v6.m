%% ============================================
%  Bilateral neurons: weighted (theta,phi) tiling & outlines
%  - Tm3/T4a columns now store theta/phi as N×2: [value, weight], sum(weights)=1
%  - Mi1 columns keep original (scalar per column)
%  - Each synapse maps to nearest column; for Tm3/T4a we add the whole (theta,phi) distribution
%    with its weights; for Mi1 we add a single (theta,phi) with weight=1.
%  - Tiles aggregate by SUM of weights (not raw counts).
%  - IN=orange, OUT=blue; hex tiling + per-neuron outlines + type-level outlines.
%  - Drop-in script: copy–paste and run (local helper funcs at bottom).
%% ============================================

%% 초기 설정 및 데이터 로드
clear all; close all; clc;

% Reference Column data
load('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\Figure5_BL\Tm3_voxels.mat')
load('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\Figure5_BL\Mi1_voxels.mat')
load('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\Figure5_BL\T4a_voxels.mat')

opt = detectImportOptions('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\connections_no_threshold.csv');
opt = setvartype(opt,'pre_root_id','int64');
opt = setvartype(opt,'post_root_id','int64');
FAFBConnections = readtable('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\connections_no_threshold.csv',opt);

opt = detectImportOptions('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\fafb_v783_princeton_synapse_table.csv');
opt = setvartype(opt,'pre_root_id_720575940','int64');
opt = setvartype(opt,'post_root_id_720575940','int64');
FAFB_synapse_coordinates = readtable('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\fafb_v783_princeton_synapse_table.csv',opt);
FAFB_synapse_coordinates.Properties.VariableNames(11) = "pre_root_id";
FAFB_synapse_coordinates.Properties.VariableNames(12) = "post_root_id";
FAFB_synapse_coordinates.pre_root_id  = FAFB_synapse_coordinates.pre_root_id  + int64(720575940000000000);
FAFB_synapse_coordinates.post_root_id = FAFB_synapse_coordinates.post_root_id + int64(720575940000000000);

Me_R_Faces    = readmatrix('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\Neuropil\Me_R_faces.csv') + 1;
Me_R_Vertices = readmatrix('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\Neuropil\Me_R_vertices.csv');
Me_L_Faces    = readmatrix('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\Neuropil\Me_L_faces.csv') + 1;
Me_L_Vertices = readmatrix('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\Neuropil\Me_L_vertices.csv');

Lo_R_Faces    = readmatrix('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\Neuropil\LO_R_faces.csv')+1;
Lo_R_Vertices = readmatrix('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\Neuropil\LO_R_vertices.csv');
Lo_L_Faces    = readmatrix('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\Neuropil\LO_L_faces.csv')+1;
Lo_L_Vertices = readmatrix('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\Neuropil\LO_L_vertices.csv');

LoP_R_Faces    = readmatrix('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\Neuropil\LoP_R_faces.csv')+1;
LoP_R_Vertices = readmatrix('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\Neuropil\LoP_R_vertices.csv');
LoP_L_Faces    = readmatrix('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\Neuropil\LoP_L_faces.csv')+1;
LoP_L_Vertices = readmatrix('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\Neuropil\LoP_L_vertices.csv');

load('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\Figure1_Overall\Bilateral_Neurons_Thr0.mat') % Bilateral_R_NPIs 포함

load('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\Figure5_BL\Neuropil_PCA_Basis.mat')


%% 색상 세팅
n = 256; startColor = [1,1,1];
endColorOut_orange = [0.8500, 0.3250, 0.0980]; % IN(이번 요구): 오렌지
endColorIn_blue    = [0, 0.4470, 0.7410];      % OUT(이번 요구): 파랑
fillColor_in  = endColorOut_orange;
fillColor_out = endColorIn_blue;

%% 분석할 뉴런 타입 리스트
neuron_types = unique(Bilateral_R_by_type.type);

% T4a column 준비 (theta,phi가 {i}=N×2 [value,weight])
T4a_R_Columns = T4a_Columns(strcmp(T4a_Columns.hemisphere,'R'),:);
T4a_L_Columns = T4a_Columns(strcmp(T4a_Columns.hemisphere,'L'),:);
T4a_R_out_syn_positions = cellfun(@table2array, T4a_R_Columns.out_syn_loc_denoise, 'UniformOutput', false); % OK
T4a_L_out_syn_positions = cellfun(@table2array, T4a_L_Columns.out_syn_loc_denoise, 'UniformOutput', false); % OK

% Tm3 column 준비 (theta,phi가 {i}=N×2 [value,weight])
Tm3_R_Columns = Tm3_Columns(strcmp(Tm3_Columns.hemisphere,'R'),:);
Tm3_L_Columns = Tm3_Columns(strcmp(Tm3_Columns.hemisphere,'L'),:);
Tm3_R_out_syn_positions = cellfun(@table2array, Tm3_R_Columns.out_syn_loc_denoise, 'UniformOutput', false); % OK
Tm3_L_out_syn_positions = cellfun(@table2array, Tm3_L_Columns.out_syn_loc_denoise, 'UniformOutput', false); % OK

% Mi1 column 준비 (기존 스칼라형)
Mi1_R_Columns = Mi1_Columns(strcmp(Mi1_Columns.hemisphere,'R'),:);
Mi1_L_Columns = Mi1_Columns(strcmp(Mi1_Columns.hemisphere,'L'),:);
Mi1_R_out_syn_positions = cellfun(@table2array, Mi1_R_Columns.out_syn_loc_denoise, 'UniformOutput', false);
Mi1_L_out_syn_positions = cellfun(@table2array, Mi1_L_Columns.out_syn_loc_denoise, 'UniformOutput', false);
Mi1_theta = Mi1_Columns.theta;  % 배경 스캐터용
Mi1_phi   = Mi1_Columns.phi;    % 배경 스캐터용

%% 반복 처리
all_neurons_by_type = struct();  % 타입별 current_neurons 저장용

for nt = 2:length(neuron_types)
    type_name = neuron_types{nt};
    fprintf('Processing: %s\n', type_name)

    this_type_idx   = strcmp(Bilateral_R_NPIs.type, type_name);
    current_neurons = table(Bilateral_R_NPIs.root_id(this_type_idx), 'VariableNames', {'root_id'});

    % --- 시냅스 위치 추출 & 집계 ---
    for i = 1:height(current_neurons)
        rid = current_neurons.root_id(i);

        InConnections_idx  = FAFBConnections.post_root_id==rid;
        InConnections      = FAFBConnections(InConnections_idx,:);
        OutConnections_idx = FAFBConnections.pre_root_id==rid;
        OutConnections     = FAFBConnections(OutConnections_idx,:);

        N_in_synapse_Me_R  = sum(InConnections.syn_count(strcmp(InConnections.neuropil,'ME_R')));
        N_in_synapse_Me_L  = sum(InConnections.syn_count(strcmp(InConnections.neuropil,'ME_L')));
        N_in_synapse_Lo_R  = sum(InConnections.syn_count(strcmp(InConnections.neuropil,'LO_R')));
        N_in_synapse_Lo_L  = sum(InConnections.syn_count(strcmp(InConnections.neuropil,'LO_L')));
        N_in_synapse_LoP_R = sum(InConnections.syn_count(strcmp(InConnections.neuropil,'LOP_R')));
        N_in_synapse_LoP_L = sum(InConnections.syn_count(strcmp(InConnections.neuropil,'LOP_L')));

        N_out_synapse_Me_R  = sum(OutConnections.syn_count(strcmp(OutConnections.neuropil,'ME_R')));
        N_out_synapse_Me_L  = sum(OutConnections.syn_count(strcmp(OutConnections.neuropil,'ME_L')));
        N_out_synapse_Lo_R  = sum(OutConnections.syn_count(strcmp(OutConnections.neuropil,'LO_R')));
        N_out_synapse_Lo_L  = sum(OutConnections.syn_count(strcmp(OutConnections.neuropil,'LO_L')));
        N_out_synapse_LoP_R = sum(OutConnections.syn_count(strcmp(OutConnections.neuropil,'LOP_R')));
        N_out_synapse_LoP_L = sum(OutConnections.syn_count(strcmp(OutConnections.neuropil,'LOP_L')));

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

    % --- 비율 평균 계산 ---
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

    % --- 시냅스 좌표 필터링 및 저장 ---
    for i = 1:height(current_neurons)
        rid = current_neurons.root_id(i);

        out_syn_idx = find(ismember(FAFB_synapse_coordinates.pre_root_id, rid));
        out_locs    = FAFB_synapse_coordinates(out_syn_idx, 1:3);
        in_syn_idx  = find(ismember(FAFB_synapse_coordinates.post_root_id, rid));
        in_locs     = FAFB_synapse_coordinates(in_syn_idx, 1:3);

        current_neurons.in_syn_loc{i}  = in_locs;
        current_neurons.out_syn_loc{i} = out_locs;

        in_array  = table2array(in_locs);
        out_array = table2array(out_locs);

        empty_loc_table = table([], [], [], 'VariableNames', {'x', 'y', 'z'});

        % 입력 시냅스 위치 저장 (threshold 0.2)
        if Mean_N_in_synapse_Me_R >= 0.2
            in_locs_Me_R = in_locs(intriangulation(Me_R_Vertices,Me_R_Faces,in_array),:);
            current_neurons.in_syn_loc_Me_R{i} = in_locs_Me_R;
        else
            current_neurons.in_syn_loc_Me_R{i} = empty_loc_table;
        end

        if Mean_N_in_synapse_Me_L >= 0.2
            in_locs_Me_L = in_locs(intriangulation(Me_L_Vertices,Me_L_Faces,in_array),:);
            current_neurons.in_syn_loc_Me_L{i} = in_locs_Me_L;
        else
            current_neurons.in_syn_loc_Me_L{i} = empty_loc_table;
        end

        if Mean_N_in_synapse_Lo_R >= 0.2
            in_locs_Lo_R = in_locs(intriangulation(Lo_R_Vertices,Lo_R_Faces,in_array),:);
            current_neurons.in_syn_loc_Lo_R{i} = in_locs_Lo_R;
        else
            current_neurons.in_syn_loc_Lo_R{i} = empty_loc_table;
        end

        if Mean_N_in_synapse_Lo_L >= 0.2
            in_locs_Lo_L = in_locs(intriangulation(Lo_L_Vertices,Lo_L_Faces,in_array),:);
            current_neurons.in_syn_loc_Lo_L{i} = in_locs_Lo_L;
        else
            current_neurons.in_syn_loc_Lo_L{i} = empty_loc_table;
        end

        if Mean_N_in_synapse_LoP_R >= 0.2
            in_locs_LoP_R = in_locs(intriangulation(LoP_R_Vertices,LoP_R_Faces,in_array),:);
            current_neurons.in_syn_loc_LoP_R{i} = in_locs_LoP_R;
        else
            current_neurons.in_syn_loc_LoP_R{i} = empty_loc_table;
        end

        if Mean_N_in_synapse_LoP_L >= 0.2
            in_locs_LoP_L = in_locs(intriangulation(LoP_L_Vertices,LoP_L_Faces,in_array),:);
            current_neurons.in_syn_loc_LoP_L{i} = in_locs_LoP_L;
        else
            current_neurons.in_syn_loc_LoP_L{i} = empty_loc_table;
        end

        % 출력 시냅스 위치 저장
        if Mean_N_out_synapse_Me_R >= 0.2
            out_locs_Me_R = out_locs(intriangulation(Me_R_Vertices,Me_R_Faces,out_array),:);
            current_neurons.out_syn_loc_Me_R{i} = out_locs_Me_R;
        else
            current_neurons.out_syn_loc_Me_R{i} = empty_loc_table;
        end

        if Mean_N_out_synapse_Me_L >= 0.2
            out_locs_Me_L = out_locs(intriangulation(Me_L_Vertices,Me_L_Faces,out_array),:);
            current_neurons.out_syn_loc_Me_L{i} = out_locs_Me_L;
        else
            current_neurons.out_syn_loc_Me_L{i} = empty_loc_table;
        end

        if Mean_N_out_synapse_Lo_R >= 0.2
            out_locs_Lo_R = out_locs(intriangulation(Lo_R_Vertices,Lo_R_Faces,out_array),:);
            current_neurons.out_syn_loc_Lo_R{i} = out_locs_Lo_R;
        else
            current_neurons.out_syn_loc_Lo_R{i} = empty_loc_table;
        end

        if Mean_N_out_synapse_Lo_L >= 0.2
            out_locs_Lo_L = out_locs(intriangulation(Lo_L_Vertices,Lo_L_Faces,out_array),:);
            current_neurons.out_syn_loc_Lo_L{i} = out_locs_Lo_L;
        else
            current_neurons.out_syn_loc_Lo_L{i} = empty_loc_table;
        end

        if Mean_N_out_synapse_LoP_R >= 0.2
            out_locs_LoP_R = out_locs(intriangulation(LoP_R_Vertices,LoP_R_Faces,out_array),:);
            current_neurons.out_syn_loc_LoP_R{i} = out_locs_LoP_R;
        else
            current_neurons.out_syn_loc_LoP_R{i} = empty_loc_table;
        end

        if Mean_N_out_synapse_LoP_L >= 0.2
            out_locs_LoP_L = out_locs(intriangulation(LoP_L_Vertices,LoP_L_Faces,out_array),:);
            current_neurons.out_syn_loc_LoP_L{i} = out_locs_LoP_L;
        else
            current_neurons.out_syn_loc_LoP_L{i} = empty_loc_table;
        end
    end

    % 입력/출력 0인 뉴런 제거
    idx_in_zero  = cell2mat(current_neurons.N_in_synpase_Optic)  == 0;
    idx_out_zero = cell2mat(current_neurons.N_out_synpase_Optic) == 0;
    current_neurons(idx_in_zero | idx_out_zero, :) = [];


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    column_selection='all';%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    % --- theta/phi 매핑 (가중 분포 반영) ---
    for i = 1:height(current_neurons)
        % ===== 입력 =====
        input_targets = {
            'Me_R',  Mi1_R_out_syn_positions, Mi1_R_Columns, 'unit';     % Mi1: 기존 방식 (단일)
            'Me_L',  Mi1_L_out_syn_positions, Mi1_L_Columns, 'unit';
            'Lo_R',  Tm3_R_out_syn_positions,  Tm3_R_Columns,  'weighted'; % Tm3: 값+빈도
            'Lo_L',  Tm3_L_out_syn_positions,  Tm3_L_Columns,  'weighted';
            'LoP_R', T4a_R_out_syn_positions, T4a_R_Columns, 'weighted'; % T4a: 값+빈도
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

            % (선택) root id 기록은 생략
            current_neurons.([region '_in']){i} = [root_id];
        end
        current_neurons.in_syn_theta{i} = theta_in_all;
        current_neurons.in_syn_phi{i}   = phi_in_all;
        current_neurons.in_weights{i}   = w_in_all;   % ★ 가중치 저장

        % ===== 출력 =====
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
        current_neurons.out_weights{i}   = w_out_all;   % ★ 가중치 저장
    end

    %% ===== 개별 뉴런: IN(오렌지) + OUT(파랑) 히트맵 오버레이 =====
    % - Tin/Tout의 (phi, theta, count) 점들을 2D 가중 히스토그램으로 누적
    % - 같은 축 위에 두 개의 axes를 겹쳐서, 서로 다른 colormap과 알파로 시각화
    % - 히트맵은 빈 값/희소 구역을 자연스럽게 처리 (가중치=0이면 투명)

    % 히트맵 해상도(간격, 단위: 도)
    dphi   = 2;                       % 가로(φ) 간격
    dtheta = 2;                       % 세로(θ) 간격
    phi_edges   = -180:dphi:180;      phi_centers   = (phi_edges(1:end-1) + phi_edges(2:end))/2;
    theta_edges =  -90:dtheta: 90;    theta_centers = (theta_edges(1:end-1) + theta_edges(2:end))/2;
    
    smoothing_on = true;     % 스무딩 끄려면 false
    sigma_phi_deg   = 9.5/2.355;     % 가로(φ) 방향 스무딩 정도
    sigma_theta_deg = 9.5/2.355;     % 세로(θ) 방향 스무딩 정도

    draw_contours   = true;
    contour_levels  = [0.1 0.1];

    % colormap (흰→오렌지, 흰→파랑)
    nC = 256;
    cmap_orange = [linspace(1,0.85,nC)', linspace(1,0.325,nC)', linspace(1,0.098,nC)'];
    cmap_blue   = [linspace(1,0.00,nC)', linspace(1,0.447,nC)', linspace(1,0.741,nC)'];

    % 감마(고빈도 강조/완화)
    gamma_in  = 1.0;
    gamma_out = 1.0;

    for i = 1:height(current_neurons)

        % ============================================================
        % (추가) in_tiles / out_tiles 없으면 만들어 넣기 (가중 합)
        % ============================================================
        % IN 타일 생성
        need_make_in = ~ismember('in_tiles', current_neurons.Properties.VariableNames) ...
            || isempty(current_neurons.in_tiles{i});
        if need_make_in
            phi_in   = current_neurons.in_syn_phi{i};
            theta_in = current_neurons.in_syn_theta{i};
            w_in     = current_neurons.in_weights{i};

            % cell이면 숫자 벡터로 변환
            if iscell(phi_in),   phi = cell2mat(phi_in(:));   else, phi = phi_in(:);   end
            if iscell(theta_in), the = cell2mat(theta_in(:)); else, the = theta_in(:); end
            % 가중치 정리
            if isempty(w_in)
                w = ones(min(numel(phi), numel(the)), 1);
            else
                if iscell(w_in), w = cell2mat(w_in(:)); else, w = w_in(:); end
            end

            % 길이 정합
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

        % (선택) OUT 타일도 없을 경우 생성
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

        % --- 데이터 준비 ---
        Tin  = current_neurons.in_tiles{i};   % columns: {'phi','theta','count'}
        Tout = current_neurons.out_tiles{i};

        % 방어
        if isempty(Tin);  Tin = array2table(zeros(0,3),'VariableNames',{'phi','theta','count'}); end
        if isempty(Tout); Tout = array2table(zeros(0,3),'VariableNames',{'phi','theta','count'}); end

        % --- 가중 2D 히스토그램
        H_in  = local_weighted_hist2(Tin.theta,  Tin.phi,  Tin.count,  theta_edges, phi_edges);
        H_out = local_weighted_hist2(Tout.theta, Tout.phi, Tout.count, theta_edges, phi_edges);

        % ================== 가우시안 스무딩 옵션 ==================


        if smoothing_on
            % bin 간격(도) -> bin단위 시그마로 변환


            % φ: 원형(circular) 패딩, θ: 가장자리 고정
            H_in = local_gauss_smooth2_circ(H_in, [], [], struct( ...
                'mode','deg', ...
                'sigma_theta_deg', 9.5/2.355, ...  % Δρ=9.5° 기준
                'sigma_phi_deg',   9.5/2.355, ...
                'bin_theta_deg',   2, ...
                'bin_phi_deg',     2));
            H_out = local_gauss_smooth2_circ(H_out, [], [], struct( ...
                'mode','deg', ...
                'sigma_theta_deg', 9.5/2.355, ...  % Δρ=9.5° 기준
                'sigma_phi_deg',   9.5/2.355, ...
                'bin_theta_deg',   2, ...
                'bin_phi_deg',     2));
            % H_in = local_gauss_smooth2_circ(H_in, [], [], struct( ...
            %     'mode','deg', ...
            %     'sigma_theta_deg', 10, ...  % Δρ=5° 기준
            %     'sigma_phi_deg',   10, ...
            %     'bin_theta_deg',   2, ...
            %     'bin_phi_deg',     2));
            % H_out = local_gauss_smooth2_circ(H_out, [], [], struct( ...
            %     'mode','deg', ...
            %     'sigma_theta_deg', 10, ...  % Δρ=5° 기준
            %     'sigma_phi_deg',   10, ...
            %     'bin_theta_deg',   2, ...
            %     'bin_phi_deg',     2));

        end
        % ================================================================

        % 정규화 (각각 독립 정규화; 공통스케일 원하면 max_in_out = max([H_in(:);H_out(:)]))
        max_in  = max(H_in(:));  if max_in<=0,  max_in = 1; end
        max_out = max(H_out(:)); if max_out<=0, max_out = 1; end
        Hin  = (H_in  / max_in ).^gamma_in;
        Hout = (H_out / max_out).^gamma_out;

        % 너무 작은 수치는 투명 처리되도록 클리핑(선택)
        eps_clip = 1e-6;
        Hin(Hin < eps_clip)   = 0;
        Hout(Hout < eps_clip) = 0;

        % ========== [중요] 아래 타입-수준에서 재사용할 결과 저장 ==========
        current_neurons.Hin_norm{i}  = Hin;     % 정규화 + 감마 + 클리핑 후 결과
        current_neurons.Hout_norm{i} = Hout;
        % ===============================================================

        % --- Figure & axes 두 겹(오버레이) ---
        fig = figure('Name', sprintf('IN/OUT heatmap overlay: %s #%d', type_name, i)); clf
        set(fig,'Color','w');

        % 아래 축: IN 히트맵
        ax1 = axes; hold(ax1,'on');
        % 배경 Mi1 스캐터
        scatter(ax1, Mi1_phi, Mi1_theta, 6, [0.85 0.85 0.85], 'filled', 'MarkerFaceAlpha', 0.20);
        % IN heatmap
        hi = imagesc(ax1, 'XData', [phi_edges(1) phi_edges(end)], ...
            'YData', [theta_edges(1) theta_edges(end)], ...
            'CData', Hin);
        set(ax1,'YDir','normal'); axis(ax1,'equal'); xlim(ax1,[-180 180]); ylim(ax1,[-90 90]);
        colormap(ax1, cmap_orange);
        % IN 투명도: 값이 0인 곳은 투명
        set(hi, 'AlphaData', Hin, 'AlphaDataMapping','none');

        % 위 축: OUT 히트맵 (투명 배경)
        ax2 = axes('Position', get(ax1,'Position'), 'Color','none'); hold(ax2,'on');
        ho = imagesc(ax2, 'XData', [phi_edges(1) phi_edges(end)], ...
            'YData', [theta_edges(1) theta_edges(end)], ...
            'CData', Hout);
        set(ax2,'YDir','normal'); axis(ax2,'equal'); xlim(ax2,[-180 180]); ylim(ax2,[-90 90]);
        colormap(ax2, cmap_blue);
        % OUT 투명도
        set(ho, 'AlphaData', Hout, 'AlphaDataMapping','none');

        % ================== IN/OUT 등치선(특정 값 연결) ==================


        if draw_contours
            col_in_1  = cmap_orange(end,:);
            col_out_1 = cmap_blue(end,:);

            % 등치선은 ax2(맨 위 축)에!
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

        % 격자/라벨/제목/축 설정 (두 축 동기화)
        linkaxes([ax1,ax2],'xy');
        set([ax1,ax2],'XTick',-180:45:180,'YTick',-90:45:90,'TickDir','out');
        xlabel(ax1,'\phi'); ylabel(ax1,'\theta');
        title(ax1, sprintf('IN (orange) + OUT (blue) heatmaps: %s (root id: %d)', ...
            type_name, current_neurons.root_id(i)));
        set(ax2,'XTick',[],'YTick',[]);  % 위축 눈금 숨김

        % 저장
        save_dir = fullfile(pwd, 'NeuronFigures', type_name);
        if ~exist(save_dir, 'dir'); mkdir(save_dir); end
        filename_base = sprintf('%s_INOUT_heatmap_%d', type_name, current_neurons.root_id(i));
        print(fig, '-dpng', '-r300', fullfile(save_dir, [filename_base '.png']));
        savefig(fig, fullfile(save_dir, [filename_base '.fig']));
        print(gcf,'-depsc2','-vector',fullfile(save_dir, [filename_base '.eps']))

        % close(fig);
    end

    %% ===== [타입 수준] 뉴런별 등치선 한 번에 그리기 (각 뉴런 고유 색) =====
    if height(current_neurons) > 0
        % 등치선 레벨(정규화된 히트맵 기준) — 위와 동일하게

        % (개별 섹션에서 쓴 bin 정의를 그대로 사용)
        % dphi, dtheta, phi_edges, theta_edges, phi_centers, theta_centers는 위에서 이미 정의됨

        % 뉴런별 색상 (distinguishable_colors; 없으면 lines로 대체)
        try
            cmap = distinguishable_colors(height(current_neurons));
        catch
            warning('distinguishable_colors가 경로에 없어 기본 lines 컬러를 사용합니다.');
            cmap = lines(height(current_neurons));
        end

        fig = figure('Name', sprintf('All neuron contours (per-neuron colors): %s', type_name));
        set(fig,'Color','w');
        ax = axes; hold(ax,'on');
        scatter(ax, Mi1_phi, Mi1_theta, 6, [0.85 0.85 0.85], 'filled', 'MarkerFaceAlpha', 0.20);

        % (선택) 배경 틱/축
        axis(ax,'equal'); xlim(ax,[-180 180]); ylim(ax,[-90 90]);
        set(ax,'XTick',-180:45:180,'YTick',-90:45:90,'TickDir','out');
        xlabel(ax,'\phi'); ylabel(ax,'\theta');
        title(ax, sprintf('IN (solid) / OUT (dashed) contours by neuron: %s', type_name));

        % 범례용 프록시(뉴런별 1개만 만들고 "OUT은 동일색 점선"이라고 표시)
        legend_entries = gobjects(height(current_neurons),1);
        legend_labels  = strings(height(current_neurons),1);

        for ii = 1:height(current_neurons)
            % ---------- 위 섹션에서 저장한 히트맵(최종본)을 그대로 사용 ----------
            Hin = []; Hout = [];
            if ismember('Hin_norm', current_neurons.Properties.VariableNames)
                Hin = current_neurons.Hin_norm{ii};
            end
            if ismember('Hout_norm', current_neurons.Properties.VariableNames)
                Hout = current_neurons.Hout_norm{ii};
            end

            % IN 등치선 (실선)
            if ~isempty(Hin) && any(Hin(:) >= min(contour_levels))
                contour(ax, phi_centers, theta_centers, Hin,  contour_levels, ...
                    'LineColor', cmap(ii,:), 'LineStyle','-',  'LineWidth', 1.6, ...
                    'HandleVisibility','off', 'Clipping','off');
            end
            % OUT 등치선 (점선)
            if ~isempty(Hout) && any(Hout(:) >= min(contour_levels))
                contour(ax, phi_centers, theta_centers, Hout, contour_levels, ...
                    'LineColor', cmap(ii,:), 'LineStyle','--', 'LineWidth', 1.6, ...
                    'HandleVisibility','off', 'Clipping','off');
            end

            % 범례 프록시(실선으로 1개만; 같은 색 점선=OUT임을 제목/축제목으로 명시)
            legend_entries(ii) = plot(ax, NaN,NaN,'-','Color',cmap(ii,:),'LineWidth',1.8);
            legend_labels(ii)  = sprintf('root %d', current_neurons.root_id(ii));
        end

        % % 범례: "IN=실선, OUT=같은색 점선"
        % lgd = legend(ax, legend_entries, cellstr(legend_labels), ...
        %     'Location','southoutside','Orientation','vertical','Box','off');
        title('Neuron (IN solid / OUT dashed)');

        % 저장
        save_dir = fullfile(pwd, 'NeuronFigures', type_name);
        if ~exist(save_dir, 'dir'); mkdir(save_dir); end
        filename_base = sprintf('%s_AllNeurons_Contours_PerNeuronColors_REUSED', type_name);
        print(fig, '-dpng', '-r300', fullfile(save_dir, [filename_base '.png']));
        savefig(fig, fullfile(save_dir, [filename_base '.fig']));
        print(gcf,'-depsc2','-vector',fullfile(save_dir, [filename_base '.eps']))

    end


    close all;
    safe_type_name = matlab.lang.makeValidName(type_name);
    all_neurons_by_type.(safe_type_name) = current_neurons;

end

% 결과 저장
save('NeuronAnalysis_byType.mat', 'all_neurons_by_type');

%% ---------- 헬퍼 함수들 (스크립트 끝에 배치) ----------
function [theta, phi, w, src_root_id] = map_pos_to_weighted_angles(pos_array, ref_syn, ref_columns, mode, column_selection,k,PCA_Basis)
% mode: 'weighted' (Tm3/T4a) | 'unit' (Mi1)
% column_selection: 'all' | 'unique'
%   - 'all'    : 각 시냅스가 매칭한 컬럼을 중복 포함하여 전부 저장
%   - 'unique' : 같은 컬럼이 여러 번 선택돼도 1회만 저장
%
% 반환:
%   theta, phi : 이어붙인 각도 값들 (열 벡터)
%   w          : 각도별 가중치 (열 벡터)
%   src_root_id: 각 행이 어느 참조 컬럼에서 왔는지 나타내는 식별자(열 벡터)
%                우선순위: ref_columns.root_id > id > column_id > (없으면 해당 컬럼 인덱스)
%%%% keys
%ME__Dm6__dendrite

%LO__LT1a__dendrite

% LOP__T5a__axon
% key = 'LO__LT1a__dendrite';  % 저장된 필드명 (위에서 출력된 key 참고)
% muR   = PCA_Basis.(key).R.mean;
% U_R   = PCA_Basis.(key).R.coeffs;
% P     = [x,y,z];                     % 1x3 또는 Nx3
% P_pca = (P - muR) * U_R;             % 해당 기준 좌표계로 변환

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



%% ---- (NEW) Step 1: 표면 위 평균거리 기반으로 "가장 가까운 ref 컬럼" 인덱스 계산 ----
% 전제: 위에서 k 값으로 key/PCA_Mean/PCA_coeffs/P를 이미 선택함
% P 길이에 따라 Poly22(최단거리) vs Poly33(수직드롭) 구분
is_poly22 = numel(P) == 6;
is_poly33 = numel(P) == 10;
if ~(is_poly22 || is_poly33)
    error('P must be length 6 (poly22) or 10 (poly33).');
end

% poly33 z 평가 함수
poly33_eval = @(PP, x, y) ( ...
    PP(1) + PP(2)*x + PP(3)*y + PP(4)*x.^2 + PP(5)*x.*y + PP(6)*y.^2 + ...
    PP(7)*x.^3 + PP(8)*x.^2.*y + PP(9)*x.*y.^2 + PP(10)*y.^3 );

% 수치옵션 (원하면 opts 인자로 노출 가능)
cp_opts.tol   = 1e-12;
cp_opts.maxit = 100;

% 0) ref_syn 각 뉴런별 표면 위 점들 미리 계산 (PCA좌표계에서)
n_ref = numel(ref_syn);
ref_onSurf = cell(n_ref,1);     % 각 셀: (m_i x 3) 표면 위 점들
for i = 1:n_ref
    Xi = ref_syn{i};
    if isempty(Xi), ref_onSurf{i} = []; continue; end
    Xi = double(Xi);
    % PCA 좌표계로
    Xi_pca = (Xi - PCA_Mean) * PCA_coeffs;  % (m_i x 3)
    if is_poly22
        % 각 점을 진짜 최단거리점으로 사영
        Ri = zeros(size(Xi_pca));
        for j = 1:size(Xi_pca,1)
            [xs, ys, zs] = closest_point_poly22(P, Xi_pca(j,:), cp_opts);
            Ri(j,:) = [xs, ys, zs];
        end
    else
        % poly33: z만 f33(x,y)로 수직 드롭
        x = Xi_pca(:,1); y = Xi_pca(:,2);
        z = poly33_eval(P, x, y);
        Ri = [x, y, z];
    end
    ref_onSurf{i} = Ri;
end

% 1) pos_array의 각 점도 PCA→표면 위 점으로 사영
N = size(pos_array, 1);
nearest_idx = zeros(N,1);
for j = 1:N
    pj = double(pos_array(j,:));
    pj_pca = (pj - PCA_Mean) * PCA_coeffs;   % 1x3
    if is_poly22
        [qx, qy, qz] = closest_point_poly22(P, pj_pca, cp_opts);
        q = [qx, qy, qz];
    else
        % poly33: z만 f33(x,y)로 수직 드롭
        qx = pj_pca(1); qy = pj_pca(2);
        qz = poly33_eval(P, qx, qy);
        q  = [qx, qy, qz];
    end

    % 2) 각 ref 뉴런 i와의 "표면상 평균 3D 거리" 계산 → 최소 i 선택
    best_i = 1; best_d = inf;
    for i = 1:n_ref
        Ri = ref_onSurf{i};
        if isempty(Ri), continue; end
        d = mean( sqrt(sum( (Ri - q).^2, 2 )) );  % PCA좌표계에서의 3D 유클리드 거리 평균
        if d < best_d
            best_d = d;
            best_i = i;
        end
    end
    nearest_idx(j) = best_i;
end

%% ---- (UNCHANGED) Step 2: column_selection에 따라 idx_list 산출 ----
switch lower(column_selection)
    case 'unique'
        idx_list = unique(nearest_idx, 'stable');   % 등장 순서 유지
    otherwise % 'all'
        idx_list = nearest_idx;
end

theta = []; phi = []; w = []; src_root_id = int64([]);


% (이하 기존 코드 유지) 미리 id 필드명 결정
vars     = ref_columns.Properties.VariableNames;


% (이하 기존 코드 유지) idx_list 순회 시 루프변수 명 충돌 방지
for ii = 1:numel(idx_list)
    min_idx = idx_list(ii);

    % 참조 ID 결정
    rid = ref_columns.root_id(min_idx);


    switch mode
        case 'weighted'  % (기존) theta{idx}, phi{idx} = N×2 [value, weight]
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

        case 'unit'      % (기존) 스칼라 좌표, 가중치=1
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
% - histcounts2에서 'Weights' 미지원 버전 대응용
% - X축=theta(행), Y축=phi(열) 규약 유지
%
% 입력: theta_vals, phi_vals, weights : 모두 컬럼 벡터(double)
%       theta_edges, phi_edges        : bin edge 벡터
% 출력: H (numel(theta_edges)-1) x (numel(phi_edges)-1) 가중 합 행렬

% 안전 처리: 테이블형 입력이면 배열로
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

% 길이 정합
n = min([numel(theta_vals), numel(phi_vals), numel(weights)]);
theta_vals = theta_vals(1:n);
phi_vals   = phi_vals(1:n);
weights    = weights(1:n);

% 유효치만 사용 (NaN/Inf 제거)
valid = isfinite(theta_vals) & isfinite(phi_vals) & isfinite(weights);
theta_vals = theta_vals(valid);
phi_vals   = phi_vals(valid);
weights    = weights(valid);

% bin 인덱스
[~,~,ib_theta] = histcounts(theta_vals, theta_edges); % 0은 바깥
[~,~,ib_phi]   = histcounts(phi_vals,   phi_edges);

good = (ib_theta>0) & (ib_phi>0);
ib_theta = ib_theta(good);
ib_phi   = ib_phi(good);
weights  = weights(good);

nRow = numel(theta_edges)-1; % theta bins -> 행
nCol = numel(phi_edges)-1;   % phi bins   -> 열

% 가중 합산
H = accumarray([ib_theta, ib_phi], weights, [nRow, nCol], @sum, 0);
end


% function Hs = local_gauss_smooth2_circ(H, sigma_theta_bins, sigma_phi_bins)
% % 가우시안 스무딩(분리가능 커널). φ방향은 원형(circular) 패딩, θ방향은 가장자리 고정(replicate)
% % 입력:
% %   H                 : (nTheta x nPhi) 히스토그램
% %   sigma_theta_bins  : θ방향 가우시안 표준편차 (bin 단위)
% %   sigma_phi_bins    : φ방향 가우시안 표준편차 (bin 단위)
% 
% if nargin < 2 || isempty(sigma_theta_bins), sigma_theta_bins = 1; end
% if nargin < 3 || isempty(sigma_phi_bins),   sigma_phi_bins   = 1; end
% 
% % 커널 반경 = ceil(3*sigma) (odd size 보장)
% r_th = max(1, ceil(3*sigma_theta_bins));
% r_ph = max(1, ceil(3*sigma_phi_bins));
% th   = -r_th:r_th;
% ph   = -r_ph:r_ph;
% 
% g_th = exp(-(th.^2) / (2*sigma_theta_bins^2)); g_th = g_th / sum(g_th);
% g_ph = exp(-(ph.^2) / (2*sigma_phi_bins^2));   g_ph = g_ph / sum(g_ph);
% 
% [nRow, nCol] = size(H);
% 
% % --- θ 방향 패딩: 가장자리 값 복제 (replicate) ---
% topPad    = repmat(H(1,:),  r_th, 1);
% bottomPad = repmat(H(end,:), r_th, 1);
% H_pad_th  = [topPad; H; bottomPad];
% 
% % --- φ 방향 패딩: 원형(circular) ---
% leftPad   = H_pad_th(:, nCol-r_ph+1:nCol);
% rightPad  = H_pad_th(:, 1:r_ph);
% H_pad     = [leftPad, H_pad_th, rightPad];
% 
% % --- 분리가능 컨볼루션: 먼저 θ(세로), 다음 φ(가로) ---
% % θ 방향(세로) 컨볼브
% H_tmp = conv2(g_th.', 1, H_pad, 'same');   % (길이 x 1) 커널
% % φ 방향(가로) 컨볼브
% H_sm  = conv2(1, g_ph, H_tmp, 'same');     % (1 x 길이) 커널
% 
% % --- 원래 영역으로 크롭 ---
% Hs = H_sm(r_th+1:r_th+nRow, r_ph+1:r_ph+nCol);
% end

function Hs = local_gauss_smooth2_circ(H, sigma_theta_bins, sigma_phi_bins, opts)
% 가우시안 스무딩(분리가능 커널). φ는 circular, θ는 replicate.
% 기존 시그니처와 호환되며, opts로 단위/파라미터를 유연하게 지정.
%
% 입력:
%   H                 : (nTheta x nPhi) 히스토그램 (θ x φ)
%   sigma_theta_bins  : θ방향 σ (bin 단위)  [모드에 따라 무시될 수 있음]
%   sigma_phi_bins    : φ방향 σ (bin 단위)  [모드에 따라 무시될 수 있음]
%   opts (struct, optional) :
%       .mode            : 'bins' | 'deg' | 'accept' (default='bins')
%                          'bins'   : sigma_*_bins 그대로 사용
%                          'deg'    : sigma_*_deg와 bin_*_deg로 bins 환산
%                          'accept' : acceptance_deg(=FWHM)로부터 σ 계산
%       .sigma_theta_deg : θ방향 σ (deg, mode='deg'일 때)
%       .sigma_phi_deg   : φ방향 σ (deg, mode='deg'일 때)
%       .acceptance_deg  : Δρ = FWHM (deg, mode='accept'일 때; 등방 가정)
%       .bin_theta_deg   : θ bin 크기 (deg, mode='deg'/'accept' 필수)
%       .bin_phi_deg     : φ bin 크기 (deg, mode='deg'/'accept' 필수)
%       .truncate        : 커널 반경 계수 k (default=3; 반경=ceil(k*σ))
%
% 출력:
%   Hs : 스무딩된 히스토그램

    if nargin < 2 || isempty(sigma_theta_bins), sigma_theta_bins = 1; end
    if nargin < 3 || isempty(sigma_phi_bins),   sigma_phi_bins   = 1; end
    if nargin < 4 || isempty(opts), opts = struct(); end
    if ~isfield(opts,'mode') || isempty(opts.mode), opts.mode = 'bins'; end
    if ~isfield(opts,'truncate') || isempty(opts.truncate), opts.truncate = 3; end

    [nRow, nCol] = size(H);

    % --- 단위/모드 처리 ---
    switch lower(opts.mode)
        case 'bins'
            sig_th_bins = max(eps, sigma_theta_bins);
            sig_ph_bins = max(eps, sigma_phi_bins);

        case 'deg'
            assert(isfield(opts,'sigma_theta_deg') && isfield(opts,'sigma_phi_deg') && ...
                   isfield(opts,'bin_theta_deg')   && isfield(opts,'bin_phi_deg'), ...
                   'mode="deg"일 때 sigma_*_deg와 bin_*_deg가 필요합니다.');
            sig_th_bins = max(eps, opts.sigma_theta_deg / opts.bin_theta_deg);
            sig_ph_bins = max(eps, opts.sigma_phi_deg   / opts.bin_phi_deg);

        case 'accept'
            % Δρ(=FWHM) -> σ(deg) -> σ(bins)
            assert(isfield(opts,'acceptance_deg') && isfield(opts,'bin_theta_deg') && isfield(opts,'bin_phi_deg'), ...
                   'mode="accept"일 때 acceptance_deg, bin_*_deg가 필요합니다.');
            sigma_deg = opts.acceptance_deg / 2.355;   % σ = FWHM/2.355
            sig_th_bins = max(eps, sigma_deg / opts.bin_theta_deg);
            sig_ph_bins = max(eps, sigma_deg / opts.bin_phi_deg);

        otherwise
            error('Unknown mode: %s', opts.mode);
    end

    % σ가 매우 작으면(≈0 bins) 원본 반환
    if sig_th_bins < 1e-6 && sig_ph_bins < 1e-6
        Hs = H; return;
    end

    % --- 커널 생성 (정규화: 에너지 보존) ---
    r_th = max(1, ceil(opts.truncate * sig_th_bins));
    r_ph = max(1, ceil(opts.truncate * sig_ph_bins));
    th   = -r_th:r_th;
    ph   = -r_ph:r_ph;

    g_th = exp(-(th.^2) / (2*sig_th_bins^2)); g_th = g_th / sum(g_th);
    g_ph = exp(-(ph.^2) / (2*sig_ph_bins^2)); g_ph = g_ph / sum(g_ph);

    % --- θ 패딩: replicate ---
    topPad    = repmat(H(1,:),  r_th, 1);
    bottomPad = repmat(H(end,:), r_th, 1);
    H_pad_th  = [topPad; H; bottomPad];

    % --- φ 패딩: circular ---
    leftPad   = H_pad_th(:, nCol-r_ph+1:nCol);
    rightPad  = H_pad_th(:, 1:r_ph);
    H_pad     = [leftPad, H_pad_th, rightPad];

    % --- 분리 가능 컨볼루션 (세로 -> 가로) ---
    H_tmp = conv2(g_th.', 1, H_pad, 'same');   % θ 방향
    H_sm  = conv2(1, g_ph, H_tmp, 'same');     % φ 방향

    % --- 크롭 ---
    Hs = H_sm(r_th+1:r_th+nRow, r_ph+1:r_ph+nCol);
end
