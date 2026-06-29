%% fig_S5D_E_F_G_H_MeMe_simulation
%  Rate-model simulation of the MeMe / Sm / Dm2 / MeTu network driven by a
%  rolling visual stimulus, comparing the intact network ("Normal") with one in
%  which the MeMe_e01 <-> MeMe_e02 interconnections are cut ("Ablated").
%
%  Dynamics:  tau dr/dt = -r + f( baseline + drive(t) + W r ),  f = tanh,
%             r clamped to >= 0. W is the neuron-level signed-synapse matrix
%             from s20_MeMe_connection.m; drive(t) is the per-type RF
%             (s21_MeMe_generate_RFs.m) dotted with the rolling stimulus
%             (s19_MeMe_roll_stimuli.m).
%
%  Figures (run once per stimulus; set stim_name below):
%   - Figure 1 (paper panel S5D) - RF-quadrant subgroup color key (stim-independent).
%   - Figures 2 & 3 (paper panels S5E [cw] / S5G [ccw]) - per-type mean response by
%     RF subgroup, Normal (Figure 2) and Ablated (Figure 3).
%   - Figure 4 (paper panel S5F [cw] / S5H [ccw]) - per-type Left vs Right mean final
%     activity (line slope = L-R asymmetry), Normal (blue) vs Ablated (red). Left/right
%     and upper/lower are split at each type's RF center (median phi / median theta).
%   - Figure 5 (Video S1 [cw] / S2 [ccw]) - per-type spatial response movie over time:
%     Normal vs Ablated vs their difference (Ablated - Normal). Saved as .mp4 in the
%     project root (set make_video=false to skip).
%
%  Panel mapping:
%     stim_name = 'Roll_60deg_cw'  ->  S5D (Fig 1), S5E (Fig 2,3), S5F (Fig 4), Video S1 (Fig 5)
%     stim_name = 'Roll_60deg_ccw' ->  S5G (Fig 2,3), S5H (Fig 4),  Video S2 (Fig 5)   (Fig 1 = S5D, identical)
%
%  Prerequisites in Processed_Data\ (regenerate with the s-scripts):
%     MeMe_AllSimNeurons_Graph.mat      (s20_MeMe_connection.m)
%     <type>_RFs.mat for each sim type  (s21_MeMe_generate_RFs.m)
%     Roll_60deg_cw.mat / _ccw.mat      (s19_MeMe_roll_stimuli.m)
clear; close all; clc

set(groot, 'defaultAxesTickDir', 'out');
set(groot, 'defaultAxesTickDirMode', 'manual');

baseDir  = fileparts(fileparts(mfilename('fullpath')));  % resolved from this script's location (portable; no absolute paths)
codexDir = fullfile(baseDir, 'Codex_Data');
procDir  = fullfile(baseDir, 'Processed_Data');

%% --- Settings ---
stim_name        = 'Roll_60deg_cw';   % stimulus .mat (without extension). ccw is 'Roll_60deg_ccw'
tau_ms           = 10;
dt_ms            = 1;
fps              = 50;

% Offset to display the stimulus (roll) onset as t=0 (pre 1.0 s -> onset=1000 ms)
stim_onset_ms    = 1000;

% --- Per-type sensory drive strength (max of RF · stim) ---
sensory_strength_by_type = struct();
sensory_strength_by_type.Dm2      = 0.5;
sensory_strength_by_type.MeMe_e01 = 0.0;
sensory_strength_by_type.MeMe_e02 = 0.0;
sensory_strength_by_type.Sm07     = 0.0;
sensory_strength_by_type.MeTu1    = 0.0;
sensory_strength_by_type.MeTu3c   = 1;

% --- Per-type baseline (tonic drive) ---
baseline_by_type = struct();
baseline_by_type.MeMe_e01 = 1;
baseline_by_type.MeMe_e02 = -1;
baseline_by_type.Sm07     = 0;
baseline_by_type.Dm2      = 0;
baseline_by_type.MeTu1    = -0.5;
baseline_by_type.MeTu3c   = 0;

% --- Activation function ---  'tanh' | 'tanh_pos' | 'leaky_relu' | 'relu' | 'linear'
activation       = 'tanh';
leaky_relu_alpha = 0.01;

% --- W matrix scaling : max|W| -> W_gain ---
W_gain = 1;

% --- Option to make edge weights uniform ---
uniform_weight_by_type = false;

% --- Option to include same-type (intra-type) connections ---
include_same_type    = true;
connections_csv_path = fullfile(codexDir, 'connections_no_threshold.csv');
syn_count_min        = 5;
exclude_neuropils          = {'UNASGD'};
inhib_nts            = {'GABA','GLUT'};

% --- Ablation mode option ---  'all_output' | 'between' | 'directed'
ablation_mode = 'between';
ablate_types  = {'MeMe_e01','MeMe_e02'};
abl_from_type = 'MeMe_e02';
abl_to_type   = 'MeMe_e01';

% --- Video (Figure 5 -> Video S1 / S2) options ---
make_video        = true;                     % generate the Normal/Ablated/Diff movie (.mp4)
vid_fps           = 20;                        % playback fps
vid_nframe        = 200;                       % number of rendered frames (capped by n_steps)
vid_exclude_types = {'Dm2','Sm07','MeTu1'};    % types (rows) to exclude from the movie
resp_cmap         = 'hot';                     % colormap for the response (Normal/Ablated) columns
diff_clim         = [-0.3 0.3];                % color limit for the diff column ([] = auto, 0-symmetric)

% --- Connectivity .mat (output of s20_MeMe_connection.m) ---
graph_mat_path = fullfile(procDir, 'MeMe_AllSimNeurons_Graph.mat');

%% --- Load connectivity ---
fprintf('Loading connectivity: %s\n', graph_mat_path);
CD = load(graph_mat_path);
pair_table   = CD.pair_table;
neuron_table = CD.neuron_table;
sim_types    = CD.sim_types;
edge_table   = CD.edge_table;
N = height(neuron_table);
fprintf('  N neurons = %d  across types: %s\n', N, strjoin(sim_types, ', '));

rid2idx = containers.Map(num2cell(neuron_table.root_id), num2cell((1:N)'));

%% --- Build W (sparse, N x N) ---
fprintf('Building W matrix...\n');
nP = height(pair_table);
pre_idx_v  = zeros(nP, 1);
post_idx_v = zeros(nP, 1);
keep_p     = false(nP, 1);
for k = 1:nP
    pre_rid  = pair_table.pre_root_id(k);
    post_rid = pair_table.post_root_id(k);
    if isKey(rid2idx, pre_rid) && isKey(rid2idx, post_rid)
        pre_idx_v(k)  = rid2idx(pre_rid);
        post_idx_v(k) = rid2idx(post_rid);
        keep_p(k)     = true;
    end
end

w_pair = double(pair_table.signed_syn);

if uniform_weight_by_type
    EN = edge_table.EndNodes;
    if ~iscell(EN), EN = cellstr(EN); end
    tp_w = containers.Map('KeyType','char','ValueType','double');
    for e = 1:height(edge_table)
        tp_w([EN{e,1} '__' EN{e,2}]) = edge_table.Weight(e);
    end
    for k = 1:nP
        key_k = [pair_table.src_type{k} '__' pair_table.tgt_type{k}];
        if isKey(tp_w, key_k), w_pair(k) = tp_w(key_k); end
    end
end

W = sparse(post_idx_v(keep_p), pre_idx_v(keep_p), w_pair(keep_p), N, N);

%% --- Add same-type (intra-type) connections (optional) ---
if include_same_type
    fprintf('Adding same-type (intra-type) connections from raw CSV...\n');
    optc = detectImportOptions(connections_csv_path);
    optc = setvartype(optc, 'pre_root_id',  'int64');
    optc = setvartype(optc, 'post_root_id', 'int64');
    Cc = readtable(connections_csv_path, optc);
    Cc(Cc.syn_count < syn_count_min, :) = [];
    Cc = Cc(~ismember(Cc.neuropil, exclude_neuropils), :);

    rid_set = neuron_table.root_id;
    in_both = ismember(Cc.pre_root_id, rid_set) & ismember(Cc.post_root_id, rid_set);
    Cc = Cc(in_both, :);

    if isempty(Cc)
        fprintf('  No connections among sim neurons in raw CSV -> not added\n');
    else
        pre_i_all  = cell2mat(values(rid2idx, num2cell(Cc.pre_root_id)));
        post_i_all = cell2mat(values(rid2idx, num2cell(Cc.post_root_id)));
        same_mask  = strcmp(neuron_table.primary_type(pre_i_all), ...
                            neuron_table.primary_type(post_i_all));
        Cc = Cc(same_mask, :);

        if isempty(Cc)
            fprintf('  No same-type (pre,post) connections -> not added\n');
        else
            [uPair, ~, ic] = unique([Cc.pre_root_id Cc.post_root_id], 'rows');
            syn_sum = accumarray(ic, Cc.syn_count, [], @sum);

            [uPre, iaP] = unique(Cc.pre_root_id);
            ntPre   = Cc.nt_type(iaP);
            signPre = ones(numel(uPre), 1);
            signPre(ismember(ntPre, inhib_nts)) = -1;
            pre2sign = containers.Map(num2cell(uPre), num2cell(signPre));
            sgn_pair = cell2mat(values(pre2sign, num2cell(uPair(:,1))));

            w_same = sgn_pair .* double(syn_sum);
            pi = cell2mat(values(rid2idx, num2cell(uPair(:,1))));
            qi = cell2mat(values(rid2idx, num2cell(uPair(:,2))));

            if uniform_weight_by_type
                ptype = neuron_table.primary_type(pi);
                utp   = unique(ptype);
                for u = 1:numel(utp)
                    m = strcmp(ptype, utp{u});
                    w_same(m) = mean(w_same(m));
                end
            end

            W_same = sparse(qi, pi, w_same, N, N);
            W = W + W_same;
            fprintf('  Same-type connections added: %d pairs\n', numel(syn_sum));
        end
    end
else
    fprintf('  Same-type connections not included (include_same_type=false)\n');
end

% Normalize: max|W| -> W_gain
maxAbsW = full(max(abs(nonzeros(W))));
if isempty(maxAbsW) || maxAbsW == 0, maxAbsW = 1; end
W = W / maxAbsW * W_gain;
nz_vals = nonzeros(W);
fprintf('  W : nnz=%d, max|W|=%.3f, n_pos=%d, n_neg=%d\n', ...
    nnz(W), full(max(abs(nz_vals))), sum(nz_vals>0), sum(nz_vals<0));

%% --- Load per-type RFs (for sensory drive) ---
type_rf = struct();
sens_types = fieldnames(sensory_strength_by_type);
phi_centers_rf = []; theta_centers_rf = [];
PHI_rf = []; THE_rf = [];
nR_grid = 0; nC_grid = 0;

for it = 1:numel(sens_types)
    tp = sens_types{it};
    strength = sensory_strength_by_type.(tp);

    msk_tp = strcmp(neuron_table.primary_type, tp);
    idx_tp = find(msk_tp);
    n_tp   = numel(idx_tp);
    if n_tp == 0
        fprintf('  [%s] not in neuron_table -> skip\n', tp);
        continue;
    end

    rf_path = fullfile(procDir, sprintf('%s_RFs.mat', tp));
    if ~exist(rf_path, 'file')
        warning('[%s] RF file not found (%s) -> drive=0', tp, rf_path);
        continue;
    end

    fprintf('Loading %s RFs : %s (strength=%.3f, n=%d)\n', tp, rf_path, strength, n_tp);
    RFD = load(rf_path);
    var_name = sprintf('%s_neurons', tp);
    if ~isfield(RFD, var_name)
        warning('[%s] variable %s not in RF file -> drive=0', tp, var_name);
        continue;
    end
    tp_neurons = RFD.(var_name);

    if isempty(phi_centers_rf)
        phi_centers_rf   = RFD.phi_centers;
        theta_centers_rf = RFD.theta_centers;
        nR_grid = numel(theta_centers_rf);
        nC_grid = numel(phi_centers_rf);
        [PHI_rf, THE_rf] = meshgrid(phi_centers_rf, theta_centers_rf);
    elseif ~isequal(RFD.phi_centers, phi_centers_rf) ...
        || ~isequal(RFD.theta_centers, theta_centers_rf)
        warning('[%s] RF grid differs from other types -> skip', tp);
        continue;
    end

    rfd_rid_map = containers.Map(num2cell(tp_neurons.root_id), ...
                                  num2cell((1:height(tp_neurons))'));

    M_RF       = zeros(n_tp, nR_grid * nC_grid);
    rf_size_t  = zeros(n_tp, 1);
    rf_phi_t   = nan(n_tp,   1);
    rf_the_t   = nan(n_tp,   1);
    n_matched  = 0;  n_unmatched = 0;  n_empty = 0;

    for k = 1:n_tp
        rid = neuron_table.root_id(idx_tp(k));
        if isKey(rfd_rid_map, rid)
            i_rfd = rfd_rid_map(rid);
            H = tp_neurons.Hin_norm{i_rfd};
            if ~isempty(H)
                M_RF(k, :) = H(:)';
                sH = sum(H(:));
                rf_size_t(k) = sH;
                if sH > 0
                    rf_phi_t(k) = sum(H(:) .* PHI_rf(:)) / sH;
                    rf_the_t(k) = sum(H(:) .* THE_rf(:)) / sH;
                end
            else
                n_empty = n_empty + 1;
            end
            n_matched = n_matched + 1;
        else
            n_unmatched = n_unmatched + 1;
        end
    end
    rf_sum_safe = rf_size_t;
    rf_sum_safe(rf_sum_safe == 0) = 1;

    fprintf('  [%s] matched %d / %d  (unmatched: %d, empty RF: %d)\n', ...
        tp, n_matched, n_tp, n_unmatched, n_empty);

    type_rf.(tp).idx_in_nt   = idx_tp;
    type_rf.(tp).M_RF        = M_RF;
    type_rf.(tp).rf_size     = rf_size_t;
    type_rf.(tp).rf_phi_ctr  = rf_phi_t;
    type_rf.(tp).rf_the_ctr  = rf_the_t;
    type_rf.(tp).rf_sum_safe = rf_sum_safe;
end

if isempty(phi_centers_rf)
    error('No RF loaded for any type. Check <type>_RFs.mat.');
end

%% --- Load stimulus ---
stim_path = fullfile(procDir, [stim_name '.mat']);
fprintf('Loading stimulus: %s\n', stim_path);
SD  = load(stim_path);
pat = SD.pat;
[H_stim, W_stim, nFrames] = size(pat);

dt_frame_ms = 1000 / fps;
T_total_ms  = (nFrames - 1) * dt_frame_ms;
deg_per_px_phi   = 360 / W_stim;
deg_per_px_theta = 180 / H_stim;
phi_stim   = (-180 + deg_per_px_phi/2)   : deg_per_px_phi   : ( 180 - deg_per_px_phi/2);
theta_stim = ( -90 + deg_per_px_theta/2) : deg_per_px_theta : (  90 - deg_per_px_theta/2);
fprintf('  stim : %d x %d x %d frames,  T_total=%.0f ms\n', H_stim, W_stim, nFrames, T_total_ms);

% --- Detect the time window where the stimulus moves (frames change) (for the gray box) ---
frame_motion  = [0, squeeze(sum(sum(abs(diff(pat,1,3)),1),2))'];
moving_frames = find(frame_motion > 0);
if isempty(moving_frames)
    stim_move_t0 = NaN;  stim_move_t1 = NaN;
else
    f_start = max(1, moving_frames(1) - 1);
    f_end   = moving_frames(end);
    stim_move_t0 = (f_start - 1) * dt_frame_ms;
    stim_move_t1 = (f_end   - 1) * dt_frame_ms;
end

% Downsample stim to RF grid
pat_rf = zeros(nR_grid, nC_grid, nFrames);
for k = 1:nFrames
    pat_rf(:,:,k) = interp2(phi_stim, theta_stim, double(pat(:,:,k)), PHI_rf, THE_rf, 'linear', 0);
end

%% --- Drive(t) for each sensory-driven type ---
t_frame = (0:nFrames-1) * dt_frame_ms;
t_axis  = 0:dt_ms:T_total_ms;
n_steps = numel(t_axis);
drive_t = zeros(N, n_steps);

loaded_types = fieldnames(type_rf);
for it = 1:numel(loaded_types)
    tp = loaded_types{it};
    strength = sensory_strength_by_type.(tp);
    if strength == 0, continue; end
    R = type_rf.(tp);
    drive_tp_frame = zeros(numel(R.idx_in_nt), nFrames);
    for k = 1:nFrames
        s = pat_rf(:, :, k);
        drive_tp_frame(:, k) = strength * (R.M_RF * s(:)) ./ R.rf_sum_safe;
    end
    drive_tp_t = interp1(t_frame, drive_tp_frame', t_axis, 'linear', 0)';
    drive_t(R.idx_in_nt, :) = drive_t(R.idx_in_nt, :) + drive_tp_t;
end

%% --- Baseline (per-type) ---
baseline = zeros(N, 1);
fn = fieldnames(baseline_by_type);
for k = 1:numel(fn)
    msk = strcmp(neuron_table.primary_type, fn{k});
    baseline(msk) = baseline_by_type.(fn{k});
end
miss_types = setdiff(sim_types, fn);
if ~isempty(miss_types)
    warning('Types missing from baseline_by_type (baseline=0): %s', strjoin(miss_types, ', '));
end

%% --- Activation + dynamics (Normal) ---
act_fn = build_activation(activation, leaky_relu_alpha);
fprintf('Running dynamics: tau=%d ms, dt=%d ms, %d steps, activation=%s\n', tau_ms, dt_ms, n_steps, activation);
r = zeros(N, n_steps);
r(:, 1) = max(act_fn(baseline), 0);
for t = 1:n_steps-1
    I = act_fn(baseline + drive_t(:, t) + W * r(:, t));
    r(:, t+1) = max(r(:, t) + (dt_ms/tau_ms) * (-r(:, t) + I), 0);
end
final_r = r(:, end);

%% --- Ablation run ---
abl_type_mask = ismember(neuron_table.primary_type, ablate_types);
W_abl = W;
switch lower(ablation_mode)
    case 'all_output'
        W_abl(:, abl_type_mask) = 0;
    case 'between'
        W_abl(abl_type_mask, abl_type_mask) = 0;
    case 'directed'
        pre_mask  = strcmp(neuron_table.primary_type, abl_from_type);
        post_mask = strcmp(neuron_table.primary_type, abl_to_type);
        W_abl(post_mask, pre_mask) = 0;
    otherwise
        error('Unknown ablation_mode: %s', ablation_mode);
end
switch lower(ablation_mode)
    case 'all_output', abl_label = sprintf('no %s output', strjoin(ablate_types, '/'));
    case 'between',    abl_label = sprintf('%s inter-connections cut', strjoin(ablate_types, '<->'));
    case 'directed',   abl_label = sprintf('%s -> %s ablated', abl_from_type, abl_to_type);
end

r_abl = zeros(N, n_steps);
r_abl(:, 1) = max(act_fn(baseline), 0);
for t = 1:n_steps-1
    I = act_fn(baseline + drive_t(:, t) + W_abl * r_abl(:, t));
    r_abl(:, t+1) = max(r_abl(:, t) + (dt_ms/tau_ms) * (-r_abl(:, t) + I), 0);
end
final_r_abl = r_abl(:, end);

nT = numel(sim_types);

%% --- RF-quadrant subgroup definition (shared by Figures 1/2/3 / color key) ---
xonly_types = {'MeMe_e01','MeMe_e02'};
grp_label = {};  grp_idx = {};  grp_type = {};  grp_x = [];  grp_color = {};
xc = 0;
for tt = 1:nT
    tp = sim_types{tt};
    if isfield(type_rf, tp)
        R     = type_rf.(tp);
        idx_n = R.idx_in_nt(:);
        phi_n = R.rf_phi_ctr(:);
        the_n = R.rf_the_ctr(:);
    else
        idx_n = find(strcmp(neuron_table.primary_type, tp));
        phi_n = nan(numel(idx_n),1);
        the_n = nan(numel(idx_n),1);
    end
    valid  = ~isnan(phi_n) & ~isnan(the_n);

    % Split at each type's RF center (median phi = azimuth, median theta = elevation)
    % instead of absolute 0.  L/R = phi below/above the center; U/D = theta below/above.
    phi_ctr = median(phi_n(valid));   the_ctr = median(the_n(valid));

    if ismember(tp, xonly_types)
        defs = { sprintf('%s (L)', tp),   phi_n<phi_ctr & valid,                 [-1  0]; ...
                 sprintf('%s (R)', tp),   phi_n>phi_ctr & valid,                 [ 1  0] };
    else
        defs = { sprintf('%s (LU)', tp),  phi_n<phi_ctr & the_n>the_ctr & valid, [-1  1]; ...
                 sprintf('%s (LD)', tp),  phi_n<phi_ctr & the_n<the_ctr & valid, [-1 -1]; ...
                 sprintf('%s (RU)', tp),  phi_n>phi_ctr & the_n>the_ctr & valid, [ 1  1]; ...
                 sprintf('%s (RD)', tp),  phi_n>phi_ctr & the_n<the_ctr & valid, [ 1 -1] };
    end
    for gg = 1:size(defs,1)
        xc = xc + 1;
        grp_label{end+1} = defs{gg,1};        %#ok<SAGROW>
        grp_idx{end+1}   = idx_n(defs{gg,2}); %#ok<SAGROW>
        grp_type{end+1}  = tp;                %#ok<SAGROW>
        grp_x(end+1)     = xc;                %#ok<SAGROW>
        grp_color{end+1} = rf_quadrant_color(defs{gg,3}(1), defs{gg,3}(2)); %#ok<SAGROW>
    end
    xc = xc + 0.6;
end
nG = numel(grp_label);

%% --- Figure 1 (panel S5D): RF-quadrant color legend (color key) ---
fig_key = figure(1);
set(fig_key, 'Color','w','Name','RF-quadrant color key', 'Position',[250 250 820 420]);
tl_key = tiledlayout(1, 2, 'TileSpacing','compact','Padding','compact');
title(tl_key, {'Subgroup naming by receptive-field (RF) center position', ...
    ['\rm\fontsize{9}azimuth (\phi): hue (left = blue, right = red)' ...
     '   ·   elevation (\theta): brightness (upper = light, lower = dark)']}, ...
    'Interpreter','tex');

ax_k1 = nexttile; hold(ax_k1,'on'); axis(ax_k1,'equal');
quad_list = {[-1  1],'LU'; [ 1  1],'RU'; [-1 -1],'LD'; [ 1 -1],'RD'};
for q = 1:size(quad_list,1)
    sx = quad_list{q,1}(1);  sy = quad_list{q,1}(2);
    x0 = (sx<0)*-1;  y0 = (sy<0)*-1;
    col = rf_quadrant_color(sx, sy);
    rectangle(ax_k1,'Position',[x0 y0 1 1], 'FaceColor',col, 'EdgeColor','w','LineWidth',2);
    tcol = 'k';  if mean(col) < 0.5, tcol = 'w'; end
    text(ax_k1, sx*0.5, sy*0.5, quad_list{q,2}, 'Color',tcol, ...
        'HorizontalAlignment','center','FontWeight','bold','FontSize',15);
end
plot(ax_k1, [-1 1],[0 0], 'k-','LineWidth',1.2);
plot(ax_k1, [0 0],[-1 1], 'k-','LineWidth',1.2);
xlim(ax_k1,[-1 1]); ylim(ax_k1,[-1 1]);
xlabel(ax_k1,'x = \phi (azimuth)','Interpreter','tex');
ylabel(ax_k1,'y = \theta (elevation)','Interpreter','tex');
title(ax_k1,'general types  (LU / RU / LD / RD)');
set(ax_k1,'XTick',[-0.5 0.5],'XTickLabel',{'L','R'},'YTick',[-0.5 0.5],'YTickLabel',{'D','U'});

ax_k2 = nexttile; hold(ax_k2,'on'); axis(ax_k2,'equal');
for sx = [-1 1]
    x0  = (sx<0)*-1;
    col = rf_quadrant_color(sx, 0);
    rectangle(ax_k2,'Position',[x0 -1 1 2], 'FaceColor',col, 'EdgeColor','w','LineWidth',2);
    tcol = 'k';  if mean(col) < 0.5, tcol = 'w'; end
    if sx < 0, lab = 'L'; else, lab = 'R'; end
    text(ax_k2, sx*0.5, 0, lab, 'Color',tcol, ...
        'HorizontalAlignment','center','FontWeight','bold','FontSize',15);
end
plot(ax_k2, [0 0],[-1 1], 'k-','LineWidth',1.2);
xlim(ax_k2,[-1 1]); ylim(ax_k2,[-1 1]);
xlabel(ax_k2,'x = \phi (azimuth)','Interpreter','tex');
title(ax_k2, sprintf('%s  (L / R only)', strjoin(xonly_types,' / ')), 'Interpreter','none');
set(ax_k2,'XTick',[-0.5 0.5],'XTickLabel',{'L','R'},'YTick',[]);

%% --- Figures 2 & 3 (panels S5E / S5G): per-type mean response by RF subgroup ---
fig2_ylim = zeros(nT, 2);
for tt = 1:nT
    g_this = find(strcmp(grp_type, sim_types{tt}));
    vmax = eps;
    for jj = 1:numel(g_this)
        idx = grp_idx{g_this(jj)};
        if isempty(idx), continue; end
        vmax = max([vmax, max(mean(r(idx,:),1)), max(mean(r_abl(idx,:),1))]);
    end
    fig2_ylim(tt,:) = [0, vmax*1.05];
end

fig2_specs = { r,     'Normal'; ...
               r_abl, sprintf('Ablated (%s)', abl_label) };
for fp = 1:size(fig2_specs,1)
    ruse = fig2_specs{fp,1};
    tagn = fig2_specs{fp,2};
    fh = figure(fp+1);   % Figure 2 = Normal, Figure 3 = Ablated
    set(fh, 'Color','w', ...
        'Name', sprintf('Per-type mean response by RF subgroup (%s)', tagn), ...
        'Position', [80+30*(fp-1) 60-20*(fp-1) 1150 950]);
    tl = tiledlayout(nT, 1, 'TileSpacing','compact','Padding','compact');
    title(tl, sprintf('Per-type mean response by RF subgroup \\bf[%s]\\rm   (stim: %s)', ...
        strrep(tagn,'_','\_'), strrep(stim_name,'_','\_')), 'Interpreter','tex');
    for tt = 1:nT
        tp     = sim_types{tt};
        g_this = find(strcmp(grp_type, tp));
        ax = nexttile;  hold(ax, 'on');
        any_line = false;
        for jj = 1:numel(g_this)
            gg  = g_this(jj);
            idx = grp_idx{gg};
            if isempty(idx), continue; end
            m_ = mean(ruse(idx, :), 1);
            paren = regexprep(grp_label{gg}, '^[^(]*\(', '(');
            lab   = sprintf('%s  n=%d', paren, numel(idx));
            plot(ax, t_axis - stim_onset_ms, m_, 'Color', grp_color{gg}, 'LineWidth', 1.6, ...
                'DisplayName', lab);
            any_line = true;
        end
        grid(ax, 'on');  xlim(ax, [-500 2500]);  ylim(ax, fig2_ylim(tt,:));
        ylabel(ax, strrep(tp,'_','\_'), 'Interpreter', 'tex');
        yl = fig2_ylim(tt,:);
        hbox = patch(ax, [0 2000 2000 0], [yl(1) yl(1) yl(2) yl(2)], [0.85 0.85 0.85], ...
            'EdgeColor','none', 'FaceAlpha', 0.5, 'HandleVisibility','off');
        uistack(hbox, 'bottom');
        if any_line
            legend(ax, 'Location','eastoutside', 'Box','off', 'Interpreter','none');
        end
        if tt == nT, xlabel(ax, 'time from stim onset (ms)'); end
    end
end

%% --- Figure 4 (panel S5F / S5H): per-type Left vs Right final activity (Normal vs Ablated) ---
%   One panel per type x elevation.  x-axis = {Left, Right}, y = mean final r.
%   normal (blue solid) / ablated (red dashed); the line slope = left-right (L-R)
%   difference, so the change in slope = how ablation shifts the L-R asymmetry.
%     - MeTu3c / Sm07 / MeTu1   : upper and lower panels both shown.
%     - MeMe_e01 / MeMe_e02     : azimuth-only split -> upper panel only.
%   Left/right(phi) and upper/lower(theta) are split at each type's RF center
%   (median phi / median theta).  Annotation 'L-R: dn -> da (p%)' : dn,da = normal /
%   ablated left-right difference, p = percent change in |L-R|.  Dm2 excluded.
types5 = {'MeMe_e01','MeMe_e02','MeTu3c','Sm07','MeTu1'};
keep5  = ismember(types5, sim_types) & cellfun(@(t) isfield(type_rf,t), types5);
types5 = types5(keep5);
nType5 = numel(types5);
is_xo  = ismember(types5, xonly_types);    % azimuth-only -> upper panel only
colN   = [0.10 0.45 0.85];  colA = [0.85 0.15 0.15];   % normal / ablated

% Panel list (left to right) : {type, elevation 'U'/'D'}
plist = cell(0,2);
for c = 1:nType5
    plist(end+1,:) = {types5{c}, 'U'};                       %#ok<SAGROW>
    if ~is_xo(c),  plist(end+1,:) = {types5{c}, 'D'}; end     %#ok<SAGROW>
end
nPanel = size(plist,1);

fig_lr = figure(4);
set(fig_lr, 'Color','w', 'Name','Left vs Right slope (Normal vs Ablated)', ...
    'Position',[40 150 max(210*nPanel, 1100) 400]);
tl4 = tiledlayout(fig_lr, 1, nPanel, 'TileSpacing','compact', 'Padding','compact');
title(tl4, {sprintf('Left vs Right activity (mean): Normal vs Ablated   |   stim: %s', ...
    strrep(stim_name,'_','\_')), ...
    ['\rm\fontsize{9}line slope = L-R asymmetry;  normal (blue solid) vs ablated (red dashed);  ' ...
     'annotation L\rightarrowR: normal \rightarrow ablated  (% = change in |L-R|)']}, ...
    'Interpreter','tex');

hN4 = []; hA4 = [];  sj = @(n) 0.05*(rand(n,1)-0.5);
for p = 1:nPanel
    tp = plist{p,1};  ev = plist{p,2};
    R = type_rf.(tp);  idx = R.idx_in_nt(:);  phi = R.rf_phi_ctr(:);  the = R.rf_the_ctr(:);
    vmask = ~isnan(phi) & ~isnan(the);
    phi_ctr = median(phi(vmask));   the_ctr = median(the(vmask));   % type RF center
    if ev == 'U', em = the > the_ctr; elab = 'upper'; else, em = the < the_ctr; elab = 'lower'; end
    Lm = phi < phi_ctr & em & vmask;
    Rm = phi > phi_ctr & em & vmask;
    vLn = final_r(idx(Lm));      vRn = final_r(idx(Rm));        % normal
    vLa = final_r_abl(idx(Lm));  vRa = final_r_abl(idx(Rm));    % ablated

    ax = nexttile(tl4);  hold(ax,'on');
    scatter(ax, 1-0.13+sj(numel(vLn)), vLn, 7, colN, 'filled', 'MarkerFaceAlpha',0.22, 'HandleVisibility','off');
    scatter(ax, 2-0.13+sj(numel(vRn)), vRn, 7, colN, 'filled', 'MarkerFaceAlpha',0.22, 'HandleVisibility','off');
    scatter(ax, 1+0.13+sj(numel(vLa)), vLa, 7, colA, 'filled', 'MarkerFaceAlpha',0.22, 'HandleVisibility','off');
    scatter(ax, 2+0.13+sj(numel(vRa)), vRa, 7, colA, 'filled', 'MarkerFaceAlpha',0.22, 'HandleVisibility','off');

    mLn = mean(vLn); mRn = mean(vRn); mLa = mean(vLa); mRa = mean(vRa);   % means
    h1 = plot(ax, [1 2]-0.13, [mLn mRn], '-o',  'Color',colN, 'MarkerFaceColor',colN, 'LineWidth',2.2, 'MarkerSize',6);
    h2 = plot(ax, [1 2]+0.13, [mLa mRa], '--o', 'Color',colA, 'MarkerFaceColor',colA, 'LineWidth',2.2, 'MarkerSize',6);
    if isempty(hN4), hN4 = h1; hA4 = h2; end

    dn = round(mLn - mRn, 2);  da = round(mLa - mRa, 2);   % L - R, rounded to displayed digits
    if abs(dn) > 1e-6                                       % so the printed value and % stay consistent
        if abs(da) < abs(dn), arr = '\downarrow'; else, arr = '\uparrow'; end   % |L-R| shrinks / grows
        pstr = sprintf('%s%.0f%%', arr, 100*abs(abs(da)-abs(dn))/abs(dn));
    else
        pstr = 'n/a';
    end

    ylim(ax, [0 1]);  xlim(ax, [0.55 2.45]);
    set(ax, 'XTick', [1 2], 'XTickLabel', {'Left','Right'});  grid(ax,'on');
    if p == 1, ylabel(ax, 'final r (mean)'); end
    title(ax, sprintf('%s (%s)  n_L=%d n_R=%d', strrep(tp,'_','\_'), elab, numel(vLn), numel(vRn)), ...
        'Interpreter','tex');
    text(ax, 1.5, 0.95, sprintf('L-R: %+.2f \\rightarrow %+.2f  (%s)', dn, da, pstr), ...
        'HorizontalAlignment','center', 'FontSize',8.5, 'Interpreter','tex');
end
if ~isempty(hN4)
    legend([hN4 hA4], {'Normal','Ablated'}, 'Box','off', 'Location','northeast');
end

%% --- Figure 5 (Video S1 / S2): per-type response movie (Normal vs Ablated vs Diff) ---
%   Column 1 = Normal, column 2 = Ablated, column 3 = Ablated - Normal (diff).
%   Rows = types with RF centroids, excluding vid_exclude_types.
%   Each frame = one time point; dot = neuron (RF centroid), color = its value at t.
%     col 1/2 color = response r (resp_cmap, clim_vec);  col 3 = diff (diff_cmap, clim_diff).
%   Dot stacking order = |final response change| (largest on top).

% Response / diff colormaps for the movie.  r is clamped to >= 0, so the response scale starts at 0.
clim_vec  = [0, max([1, max(r(:)), eps])];
cmap_resp = feval(resp_cmap, 256);
% diff colormap: redblue diverging centered at 0 (blue #0b318f .. white .. red #e60012)
diff_blue_end = [11 49 143]/255;   diff_red_end = [230 0 18]/255;   nhalf = 128;
blue_ramp = [linspace(diff_blue_end(1),1,nhalf)', linspace(diff_blue_end(2),1,nhalf)', linspace(diff_blue_end(3),1,nhalf)'];
red_ramp  = [linspace(1,diff_red_end(1),256-nhalf)', linspace(1,diff_red_end(2),256-nhalf)', linspace(1,diff_red_end(3),256-nhalf)'];
cmap_diff = [blue_ramp; red_ramp];

% Types to draw : those with at least one RF centroid (NaN-centroid neurons are skipped)
loaded_types = fieldnames(type_rf);
draw_types = {};
for it = 1:numel(loaded_types)
    tp = loaded_types{it};
    if sum(~isnan(type_rf.(tp).rf_phi_ctr)) > 0
        draw_types{end+1} = tp;                       %#ok<SAGROW>
    end
end

vid_nframe_use = min(n_steps, vid_nframe);
draw_types_vid = draw_types(~ismember(draw_types, vid_exclude_types));
n_draw_vid     = numel(draw_types_vid);
if ~isempty(vid_exclude_types)
    fprintf('  Video-excluded types: %s\n', strjoin(vid_exclude_types, ', '));
end
if make_video && n_draw_vid > 0
    vid_path = fullfile(baseDir, ['MeMe_sim_video__' stim_name '.mp4']);
    fprintf('Generating video: %s (%d frames @ %d fps)\n', vid_path, vid_nframe_use, vid_fps);

    vid_t     = linspace(500, 3500, vid_nframe_use);   % skip initial transient: 500-3500 ms
    vid_step  = arrayfun(@(t) max(1,min(n_steps, round(t/dt_ms)+1)),       vid_t);
    vid_frame = arrayfun(@(t) max(1,min(nFrames, round(t/dt_frame_ms)+1)), vid_t);

    % diff column color range : diff_clim if given, else auto (0-symmetric, over drawn neurons)
    if ~isempty(diff_clim)
        clim_diff = diff_clim;
    else
        vid_idx_all = cell2mat(cellfun(@(tp) type_rf.(tp).idx_in_nt(:), ...
            draw_types_vid(:), 'UniformOutput', false));
        dvals     = r_abl(vid_idx_all, vid_step) - r(vid_idx_all, vid_step);
        dmax      = max([max(abs(dvals(:))), eps]);
        clim_diff = [-dmax, dmax];
    end

    figv = figure(5);
    %   Capture the InnerPosition region -> fixed 1280x720 (HD).
    set(figv, 'Color','w','Name','Normal vs Ablated vs Diff (video)', ...
        'Resize','off', 'Units','pixels', 'InnerPosition',[408 175 1280 720]);
    tlv = tiledlayout(figv, n_draw_vid, 3, 'TileSpacing','compact','Padding','compact');
    ttl = title(tlv, '', 'Interpreter','tex');

    cond_nm   = {'Normal','Ablated','Ablated - Normal'};
    cond_cmap = {cmap_resp, cmap_resp, cmap_diff};
    cond_clim = {clim_vec,  clim_vec,  clim_diff};
    nColV     = 3;
    img_h     = gobjects(n_draw_vid, nColV);
    sc_h      = gobjects(n_draw_vid, nColV);
    idx_o_all = cell(n_draw_vid, 1);       % per-type neuron order (by |change|)
    ax_resp = gobjects(1);  ax_diff = gobjects(1);
    bg0 = repmat(pat(:,:,vid_frame(1))*0.6 + 0.2, [1 1 3]);
    for it = 1:n_draw_vid
        tp  = draw_types_vid{it};
        R   = type_rf.(tp);
        idx = R.idx_in_nt(:);
        chg_f = final_r_abl(idx) - final_r(idx);     % stack order by |final change|
        [~, ord] = sort(abs(chg_f), 'ascend');
        idx_o_all{it} = idx(ord);
        phi_o = R.rf_phi_ctr(ord);  the_o = R.rf_the_ctr(ord);
        tp_esc  = strrep(tp,'_','\_');
        n_drawn = sum(~isnan(R.rf_phi_ctr));
        for c = 1:nColV
            ax = nexttile(tlv, (it-1)*nColV + c);
            img_h(it,c) = image(ax, phi_stim, theta_stim, bg0);
            set(ax,'YDir','normal'); axis(ax,'equal','tight');
            xlim(ax,[-180 180]); ylim(ax,[-90 90]); hold(ax,'on');
            rr0 = vid_col_vals(c, r, r_abl, idx_o_all{it}, vid_step(1));
            sc_h(it,c) = scatter(ax, phi_o, the_o, 26, rr0, 'filled', ...
                'MarkerEdgeColor','k','LineWidth',0.2,'MarkerFaceAlpha',0.9);
            colormap(ax, cond_cmap{c}); caxis(ax, cond_clim{c});
            set(ax,'XTick',-180:90:180,'YTick',-90:45:90);
            if it == 1,           title(ax, cond_nm{c}); end
            if c  == 1,           ylabel(ax, sprintf('%s (n=%d)', tp_esc, n_drawn), 'Interpreter','tex'); end
            if it == n_draw_vid,  xlabel(ax, '\phi'); end
            if c == 1, ax_resp = ax; elseif c == 3, ax_diff = ax; end
        end
    end
    cb_r = colorbar(ax_resp); cb_r.Layout.Tile = 'east'; cb_r.Label.String = 'Response';
    cb_d = colorbar(ax_diff); cb_d.Layout.Tile = 'east'; cb_d.Label.String = 'Response change (Ablated - Normal)';

    vw = VideoWriter(vid_path, 'MPEG-4');
    vw.FrameRate = vid_fps;
    open(vw);
    for f = 1:vid_nframe_use
        bgf = repmat(pat(:,:,vid_frame(f))*0.6 + 0.2, [1 1 3]);
        for it = 1:n_draw_vid
            idx_o = idx_o_all{it};
            for c = 1:nColV
                set(img_h(it,c), 'CData', bgf);
                set(sc_h(it,c),  'CData', vid_col_vals(c, r, r_abl, idx_o, vid_step(f)));
            end
        end
        ttl.String = sprintf('t = %d ms   |   Stim: %s', ...
            round(vid_t(f) - stim_onset_ms), strrep(stim_name,'_','\_'));
        drawnow;
        writeVideo(vw, getframe(figv));
    end
    close(vw);
    fprintf('  Video saved: %s\n', vid_path);
end


%% =========================================================
% Local functions
%% =========================================================
function f = build_activation(name, alpha)
    switch lower(name)
        case 'tanh',       f = @(x) tanh(x);
        case 'tanh_pos',   f = @(x) max(tanh(x), 0);
        case 'leaky_relu', a = alpha;  f = @(x) max(x, 0) + a * min(x, 0);
        case 'relu',       f = @(x) max(x, 0);
        case 'linear',     f = @(x) x;
        otherwise
            error('Unknown activation: %s', name);
    end
end

function c = rf_quadrant_color(sx, sy)
% Color for the RF center position (quadrant).
%   sx : sign of x(phi) (-1 = left, +1 = right)
%   sy : sign of y(theta) (-1 = down, +1 = up, 0 = no y split)
%   left/right(x) = hue family (left = blue, right = red); up/down(y) = brightness.
    if sx < 0, hue = 0.58; else, hue = 0.00; end
    if     sy > 0, sat = 0.45; val = 0.96;
    elseif sy < 0, sat = 0.95; val = 0.62;
    else,          sat = 0.75; val = 0.82;
    end
    c = hsv2rgb([hue, sat, val]);
end

function v = vid_col_vals(c, r, r_abl, rows, step)
% Per-column dot color values for the Figure 5 video.
%   c==1 normal (r), c==2 ablated (r_abl), c==3 diff (ablated - normal)
    if c == 3
        v = r_abl(rows, step) - r(rows, step);
    elseif c == 2
        v = r_abl(rows, step);
    else
        v = r(rows, step);
    end
end
