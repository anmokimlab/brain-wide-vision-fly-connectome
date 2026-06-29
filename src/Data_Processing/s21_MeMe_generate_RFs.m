%% s21_MeMe_generate_RFs
% Build the per-type input receptive fields (RFs) used by the MeMe network
% simulation (Figure S5D-H). Running this script writes one <type>_RFs.mat into
% Processed_Data\ for each simulation cell type:
%   Dm2_RFs.mat, MeMe_e01_RFs.mat, MeMe_e02_RFs.mat,
%   Sm07_RFs.mat, MeTu1_RFs.mat, MeTu3c_RFs.mat
%
% For every neuron of a type, its input (postsynaptic) synapses inside the
% medulla are mapped to visual-space (theta, phi) via the nearest Mi1 reference
% column (projected onto the PCA + parabolic medulla layer), accumulated into a
% synapse-weighted 2D histogram, Gaussian-smoothed (Drho = 9.5 deg) and
% normalized to give Hin_norm. Neurons with no medulla input synapse get
% Hin_norm = [] (the simulation then leaves their RF centroid as NaN, so they
% are not plotted). No figure is drawn here.
%
% Prerequisites (produced by the main pipeline / provided inputs):
%   - Processed_Data\Mi1_Tm3_T4a_columns.mat   (Mi1_Columns; from s16_Mi1_Tm3_T4a_matching.m)
%   - Processed_Data\neuropil_PCA_basis.mat    (PCA_Basis;  from s17_making_layers.m)
%   - Processed_Data\optic_lobe_neuropil_mesh\Me_{R,L}_{vertices,faces}.csv  (from s01)
%   - Codex_Data\fafb_v783_princeton_synapse_table.csv
%   - Codex_Data\consolidated_cell_types.csv, Codex_Data\classification.csv
% The intriangulation helper (Helper_Function\) must be on the MATLAB path.
%
% (Korean inline comments are kept from the original analysis script.)
clear all; close all; clc

baseDir  = fileparts(fileparts(mfilename('fullpath')));  % resolved from this script's location (portable; no absolute paths)
codexDir = fullfile(baseDir, 'Codex_Data');
procDir  = fullfile(baseDir, 'Processed_Data');
meshDir  = fullfile(procDir, 'optic_lobe_neuropil_mesh');
if ~exist(procDir, 'dir'), mkdir(procDir); end
addpath(fullfile(baseDir, 'Helper_Function'));

% Simulation cell types whose input RFs are needed by fig_S5D_E_F_G_H.
sim_types = {'Dm2','MeMe_e01','MeMe_e02','Sm07','MeTu1','MeTu3c'};

%% ---------------- Load reference data (once) ----------------
fprintf('Loading Mi1 reference columns / PCA basis...\n');
S = load(fullfile(procDir, 'Mi1_Tm3_T4a_columns.mat'), 'Mi1_Columns');
Mi1_Columns = S.Mi1_Columns;
S = load(fullfile(procDir, 'neuropil_PCA_basis.mat'), 'PCA_Basis');
PCA_Basis = S.PCA_Basis;

fprintf('Loading medulla mesh...\n');
Me_R_Faces    = readmatrix(fullfile(meshDir, 'Me_R_faces.csv'))   + 1;
Me_R_Vertices = readmatrix(fullfile(meshDir, 'Me_R_vertices.csv'));
Me_L_Faces    = readmatrix(fullfile(meshDir, 'Me_L_faces.csv'))   + 1;
Me_L_Vertices = readmatrix(fullfile(meshDir, 'Me_L_vertices.csv'));
meshes = struct('Me_R_F', Me_R_Faces, 'Me_R_V', Me_R_Vertices, ...
                'Me_L_F', Me_L_Faces, 'Me_L_V', Me_L_Vertices);

fprintf('Loading consolidated_cell_types.csv / classification.csv...\n');
opt_t = detectImportOptions(fullfile(codexDir, 'consolidated_cell_types.csv'));
opt_t = setvartype(opt_t, 'root_id', 'int64');
cct = readtable(fullfile(codexDir, 'consolidated_cell_types.csv'), opt_t);
opt_t = detectImportOptions(fullfile(codexDir, 'classification.csv'));
opt_t = setvartype(opt_t, 'root_id', 'int64');
cls = readtable(fullfile(codexDir, 'classification.csv'), opt_t);

fprintf('Loading synapse coordinates (large CSV)...\n');
syn_csv = fullfile(codexDir, 'fafb_v783_princeton_synapse_table.csv');
opt_t = detectImportOptions(syn_csv);
opt_t = setvartype(opt_t, 'pre_root_id_720575940',  'int64');
opt_t = setvartype(opt_t, 'post_root_id_720575940', 'int64');
FAFB_synapse_coordinates = readtable(syn_csv, opt_t);
FAFB_synapse_coordinates.Properties.VariableNames(11) = "pre_root_id";
FAFB_synapse_coordinates.Properties.VariableNames(12) = "post_root_id";
FAFB_synapse_coordinates.pre_root_id  = FAFB_synapse_coordinates.pre_root_id  + int64(720575940000000000);
FAFB_synapse_coordinates.post_root_id = FAFB_synapse_coordinates.post_root_id + int64(720575940000000000);

%% ---------------- Mi1 reference setup ----------------
Mi1_R_Columns = Mi1_Columns(strcmp(Mi1_Columns.hemisphere, 'R'), :);
Mi1_L_Columns = Mi1_Columns(strcmp(Mi1_Columns.hemisphere, 'L'), :);
ref.R_syn  = cellfun(@table2array, Mi1_R_Columns.out_syn_loc_denoise, 'UniformOutput', false);
ref.L_syn  = cellfun(@table2array, Mi1_L_Columns.out_syn_loc_denoise, 'UniformOutput', false);
ref.R_cols = Mi1_R_Columns;
ref.L_cols = Mi1_L_Columns;

%% ---------------- Generate RFs per type ----------------
for ti = 1:numel(sim_types)
    generate_input_RFs(sim_types{ti}, cct, cls, FAFB_synapse_coordinates, ...
        meshes, ref, PCA_Basis, procDir);
end
fprintf('All RFs done.\n');


%% =========================================================
% Local functions
%% =========================================================
function generate_input_RFs(primary_type, cct, cls, syn_all, meshes, ref, PCA_Basis, outDir)
% Build and save <primary_type>_RFs.mat (input RFs only) into outDir.

% Heatmap binning / smoothing
dphi   = 2;
dtheta = 2;
sigma_deg = 9.5/2.355;          % Drho = 9.5 deg

type_mask = strcmp(cct.primary_type, primary_type);
neurons   = cct(type_mask, {'root_id'});
N = height(neurons);
fprintf('[%s] %d neurons.\n', primary_type, N);
if N == 0, warning('[%s] no neurons; skip.', primary_type); return; end

% --- side from classification.csv ---
neurons.side = repmat({''}, N, 1);
[is_in_cls, cls_loc] = ismember(neurons.root_id, cls.root_id);
for i = 1:N
    if is_in_cls(i)
        s = cls.side{cls_loc(i)};
        if ischar(s) || isstring(s), s = char(s); end
        neurons.side{i} = lower(strtrim(s));
    end
end

% --- input synapses for these neurons ---
target_ids = neurons.root_id;
syn_in_all = syn_all(ismember(syn_all.post_root_id, target_ids), :);

% --- heatmap grid ---
phi_edges   = -180:dphi:180;     phi_centers   = (phi_edges(1:end-1)   + phi_edges(2:end))/2;
theta_edges =  -90:dtheta:90;    theta_centers = (theta_edges(1:end-1) + theta_edges(2:end))/2;
smopts = struct('mode','deg', ...
    'sigma_theta_deg', sigma_deg, 'sigma_phi_deg', sigma_deg, ...
    'bin_theta_deg',   dtheta,    'bin_phi_deg',   dphi);

neurons.k_region     = zeros(N, 1);              % 1 = ME_R, 2 = ME_L
neurons.n_in_syn_Me  = zeros(N, 1);
neurons.in_syn_theta = cell(N, 1);
neurons.in_syn_phi   = cell(N, 1);
neurons.in_weights   = cell(N, 1);
neurons.Hin_norm     = cell(N, 1);               % empty [] -> NaN centroid in sim

n_empty = 0;
for i = 1:N
    rid  = neurons.root_id(i);
    side = neurons.side{i};

    in_rows  = syn_in_all.post_root_id == rid;
    in_array = table2array(syn_in_all(in_rows, 1:3));

    % --- side fallback (majority of in-syns in ME_R vs ME_L) ---
    if ~ismember(side, {'right', 'left'})
        n_R = 0; n_L = 0;
        if ~isempty(in_array)
            n_R = sum(intriangulation(meshes.Me_R_V, meshes.Me_R_F, in_array) == 1);
            n_L = sum(intriangulation(meshes.Me_L_V, meshes.Me_L_F, in_array) == 1);
        end
        if n_R == 0 && n_L == 0,      side = 'right';
        elseif n_R >= n_L,            side = 'right';
        else,                         side = 'left';
        end
        neurons.side{i} = side;
    end

    if strcmpi(side, 'right')
        Me_F = meshes.Me_R_F;  Me_V = meshes.Me_R_V;
        ref_syns = ref.R_syn;  ref_cols = ref.R_cols;  k_region = 1;
    else
        Me_F = meshes.Me_L_F;  Me_V = meshes.Me_L_V;
        ref_syns = ref.L_syn;  ref_cols = ref.L_cols;  k_region = 2;
    end
    neurons.k_region(i) = k_region;

    if isempty(in_array)
        in_locs_Me = [];
    else
        in_inMe = intriangulation(Me_V, Me_F, in_array);
        in_locs_Me = in_array(in_inMe == 1, :);
    end
    neurons.n_in_syn_Me(i) = size(in_locs_Me, 1);

    if isempty(in_locs_Me)
        theta_in = []; phi_in = []; w_in = [];
    else
        [theta_in, phi_in, w_in, ~] = map_pos_to_weighted_angles( ...
            in_locs_Me, ref_syns, ref_cols, 'unit', 'all', k_region, PCA_Basis);
    end
    neurons.in_syn_theta{i} = theta_in;
    neurons.in_syn_phi{i}   = phi_in;
    neurons.in_weights{i}   = w_in;

    if isempty(in_locs_Me)
        neurons.Hin_norm{i} = [];
        n_empty = n_empty + 1;
    else
        if isempty(w_in), w_in_use = ones(numel(theta_in), 1); else, w_in_use = w_in; end
        H_in = local_weighted_hist2(theta_in, phi_in, w_in_use, theta_edges, phi_edges);
        H_in = local_gauss_smooth2_circ(H_in, [], [], smopts);
        max_in = max(H_in(:));   if max_in <= 0, max_in = 1; end
        Hin    = H_in / max_in;
        Hin(Hin < 1e-6) = 0;
        if all(Hin(:) == 0)
            neurons.Hin_norm{i} = [];  n_empty = n_empty + 1;
        else
            neurons.Hin_norm{i} = Hin;
        end
    end
end
fprintf('[%s] RFs built: %d / %d  (empty / no medulla syn: %d)\n', ...
    primary_type, N - n_empty, N, n_empty);

% --- save ---
data_path = fullfile(outDir, sprintf('%s_RFs.mat', primary_type));
out_struct = struct();
out_struct.(sprintf('%s_neurons', primary_type)) = neurons;
out_struct.phi_centers   = phi_centers;
out_struct.theta_centers = theta_centers;
out_struct.phi_edges     = phi_edges;
out_struct.theta_edges   = theta_edges;
out_struct.primary_type  = primary_type;
save(data_path, '-struct', 'out_struct', '-v7.3');
fprintf('[%s] Saved: %s\n', primary_type, data_path);
end


function [theta, phi, w, src_root_id] = map_pos_to_weighted_angles(pos_array, ref_syn, ref_columns, mode, column_selection, k, PCA_Basis)
% k = 1 (ME_R) or 2 (ME_L). mode = 'unit' (Mi1 reference; scalar per column).
if k == 1
    key       = 'ME__Dm6__dendrite';
    PCA_Mean  = PCA_Basis.(key).R.mean;
    PCA_coeffs= PCA_Basis.(key).R.coeffs;
    p00=PCA_Basis.(key).fit.R.p00; p10=PCA_Basis.(key).fit.R.p10; p01=PCA_Basis.(key).fit.R.p01;
    p20=PCA_Basis.(key).fit.R.p20; p11=PCA_Basis.(key).fit.R.p11; p02=PCA_Basis.(key).fit.R.p02;
    P=[p00 p10 p01 p20 p11 p02];
elseif k == 2
    key       = 'ME__Dm6__dendrite';
    PCA_Mean  = PCA_Basis.(key).L.mean;
    PCA_coeffs= PCA_Basis.(key).L.coeffs;
    p00=PCA_Basis.(key).fit.L.p00; p10=PCA_Basis.(key).fit.L.p10; p01=PCA_Basis.(key).fit.L.p01;
    p20=PCA_Basis.(key).fit.L.p20; p11=PCA_Basis.(key).fit.L.p11; p02=PCA_Basis.(key).fit.L.p02;
    P=[p00 p10 p01 p20 p11 p02];
else
    error('Only k=1 (ME_R) or k=2 (ME_L) supported.');
end

is_poly22 = numel(P) == 6;
cp_opts.tol   = 1e-12;
cp_opts.maxit = 100;

n_ref = numel(ref_syn);
ref_onSurf = cell(n_ref, 1);
for i = 1:n_ref
    Xi = ref_syn{i};
    if isempty(Xi), ref_onSurf{i} = []; continue; end
    Xi = double(Xi);
    Xi_pca = (Xi - PCA_Mean) * PCA_coeffs;
    if is_poly22
        Ri = zeros(size(Xi_pca));
        for j = 1:size(Xi_pca, 1)
            [xs, ys, zs] = closest_point_poly22(P, Xi_pca(j,:), cp_opts);
            Ri(j,:) = [xs, ys, zs];
        end
    else
        error('Medulla expects poly22 fit (length 6).');
    end
    ref_onSurf{i} = Ri;
end

N_pos = size(pos_array, 1);
nearest_idx = zeros(N_pos, 1);
for j = 1:N_pos
    pj_pca = (double(pos_array(j,:)) - PCA_Mean) * PCA_coeffs;
    [qx, qy, qz] = closest_point_poly22(P, pj_pca, cp_opts);
    q = [qx, qy, qz];
    best_i = 1; best_d = inf;
    for i = 1:n_ref
        Ri = ref_onSurf{i};
        if isempty(Ri), continue; end
        d = mean(sqrt(sum((Ri - q).^2, 2)));
        if d < best_d, best_d = d; best_i = i; end
    end
    nearest_idx(j) = best_i;
end

switch lower(column_selection)
    case 'unique'
        idx_list = unique(nearest_idx, 'stable');
    otherwise
        idx_list = nearest_idx;
end

theta = []; phi = []; w = []; src_root_id = int64([]);
for ii = 1:numel(idx_list)
    min_idx = idx_list(ii);
    rid = ref_columns.root_id(min_idx);
    switch mode
        case 'unit'
            theta = [theta; ref_columns.theta(min_idx)]; %#ok<AGROW>
            phi   = [phi;   ref_columns.phi(min_idx)];   %#ok<AGROW>
            w     = [w;     1];                          %#ok<AGROW>
            src_root_id = [src_root_id; rid];            %#ok<AGROW>
        otherwise
            error('Expect mode=''unit''.');
    end
end
end


function H = local_weighted_hist2(theta_vals, phi_vals, weights, theta_edges, phi_edges)
if istable(theta_vals), theta_vals = theta_vals{:,:}; end
if istable(phi_vals),   phi_vals   = phi_vals{:,:};   end
if istable(weights),    weights    = weights{:,:};    end
theta_vals = theta_vals(:);  phi_vals = phi_vals(:);
if isempty(weights), weights = ones(min(numel(theta_vals), numel(phi_vals)), 1); else, weights = weights(:); end
n = min([numel(theta_vals), numel(phi_vals), numel(weights)]);
theta_vals = theta_vals(1:n);  phi_vals = phi_vals(1:n);  weights = weights(1:n);
valid = isfinite(theta_vals) & isfinite(phi_vals) & isfinite(weights);
theta_vals = theta_vals(valid);  phi_vals = phi_vals(valid);  weights = weights(valid);
[~,~,ib_theta] = histcounts(theta_vals, theta_edges);
[~,~,ib_phi]   = histcounts(phi_vals,   phi_edges);
good = (ib_theta > 0) & (ib_phi > 0);
ib_theta = ib_theta(good);  ib_phi = ib_phi(good);  weights = weights(good);
nRow = numel(theta_edges) - 1;  nCol = numel(phi_edges) - 1;
if isempty(ib_theta)
    H = zeros(nRow, nCol);
else
    H = accumarray([ib_theta, ib_phi], weights, [nRow, nCol], @sum, 0);
end
end


function Hs = local_gauss_smooth2_circ(H, sigma_theta_bins, sigma_phi_bins, opts)
% phi = circular, theta = replicate.
if nargin < 2 || isempty(sigma_theta_bins), sigma_theta_bins = 1; end
if nargin < 3 || isempty(sigma_phi_bins),   sigma_phi_bins   = 1; end
if nargin < 4 || isempty(opts), opts = struct(); end
if ~isfield(opts,'mode')     || isempty(opts.mode),     opts.mode = 'bins'; end
if ~isfield(opts,'truncate') || isempty(opts.truncate), opts.truncate = 3;  end
[nRow, nCol] = size(H);
switch lower(opts.mode)
    case 'bins'
        sig_th_bins = max(eps, sigma_theta_bins);
        sig_ph_bins = max(eps, sigma_phi_bins);
    case 'deg'
        sig_th_bins = max(eps, opts.sigma_theta_deg / opts.bin_theta_deg);
        sig_ph_bins = max(eps, opts.sigma_phi_deg   / opts.bin_phi_deg);
    otherwise
        error('Unknown mode: %s', opts.mode);
end
if sig_th_bins < 1e-6 && sig_ph_bins < 1e-6, Hs = H; return; end
r_th = max(1, ceil(opts.truncate * sig_th_bins));
r_ph = max(1, ceil(opts.truncate * sig_ph_bins));
th   = -r_th:r_th;  ph = -r_ph:r_ph;
g_th = exp(-(th.^2) / (2*sig_th_bins^2)); g_th = g_th / sum(g_th);
g_ph = exp(-(ph.^2) / (2*sig_ph_bins^2)); g_ph = g_ph / sum(g_ph);
topPad    = repmat(H(1,:),   r_th, 1);
bottomPad = repmat(H(end,:), r_th, 1);
H_pad_th  = [topPad; H; bottomPad];
leftPad   = H_pad_th(:, nCol-r_ph+1:nCol);
rightPad  = H_pad_th(:, 1:r_ph);
H_pad     = [leftPad, H_pad_th, rightPad];
H_tmp = conv2(g_th.', 1, H_pad, 'same');
H_sm  = conv2(1, g_ph, H_tmp, 'same');
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
    x=px; y=py;
    for k=1:opts.maxit
        F=f(x,y); gx=fx(x,y); gy=fy(x,y); dz=F-pz;
        r1 = x - px + dz*gx;  r2 = y - py + dz*gy;
        if hypot(r1,r2) < opts.tol*(1+hypot(px,py)), break; end
        J11 = 1 + gx*gx + dz*fxx;   J12 = gx*gy + dz*fxy;
        J22 = 1 + gy*gy + dz*fyy;   J21 = J12;
        delta = -[J11 J12; J21 J22] \ [r1; r2];
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
