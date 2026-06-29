%% s18 — RF/PF area & center metrics from the BLP RF/PF heatmaps
% Reads the per-type RF/PF heatmaps (all_BLP_by_type) from
% Processed_Data/BLP_RFs_PFs.mat (Figures/fig_6E_G_BLP_RFs_PFs.m) and, per neuron,
% computes from the normalized heatmaps:
%   - RF/PF size (deg^2) at a contour threshold
%   - weighted centers (phi circular, theta linear)
%   - RF->PF center distance, plus signed phi/theta components
%   - the same for the phi-mirrored RF center (phi -> -phi)
% The metric columns are appended to each type's table and saved to
% Processed_Data/BLP_RFs_PFs_area_center.mat, read by
% Figures/fig_6F_H_J_BLP_area_center_plot.m.

clear; clc;

baseDir = fileparts(fileparts(mfilename('fullpath')));  % resolved from this script's location (portable; no absolute paths)

matfile_in  = fullfile(baseDir, 'Processed_Data', 'BLP_RFs_PFs.mat');             % from fig_6E_G
matfile_out = fullfile(baseDir, 'Processed_Data', 'BLP_RFs_PFs_area_center.mat');

S = load(matfile_in, 'all_BLP_by_type');
all_BLP_by_type = S.all_BLP_by_type;

% ---------- Grid / contour params (must match fig_6E_G) ----------
dphi   = 2;                         % deg (phi bin size)
dtheta = 2;                         % deg (theta bin size)
phi_edges     = -180:dphi:180;
theta_edges   =  -90:dtheta: 90;
phi_centers   = (phi_edges(1:end-1)   + phi_edges(2:end))/2;
theta_centers = (theta_edges(1:end-1) + theta_edges(2:end))/2;

contour_level = 0.1;   % Hin_norm / Hout_norm >= 0.1 defines the inside of the contour

% ---------- utils ----------
if exist('wrapTo180','file')~=2
    wrapTo180 = @(ang) mod(ang+180,360)-180;  % [-180,180)
end
circ_mean_deg = @(ang_deg, w) atan2d( sum(w.*sind(ang_deg)), sum(w.*cosd(ang_deg)) );
circ_diff_deg = @(a_deg, b_deg) wrapTo180(a_deg - b_deg);  % a-b (deg, circular)

compute_field_metrics = @(Hnorm) ...
    local_field_metrics(Hnorm, phi_centers, theta_centers, contour_level, dphi, dtheta, circ_mean_deg);

% ---------- iterate over types ----------
type_names = fieldnames(all_BLP_by_type);
n_missing_Hin = 0; n_missing_Hout = 0;

for tt = 1:numel(type_names)
    T = all_BLP_by_type.(type_names{tt});
    if height(T)==0
        all_BLP_by_type.(type_names{tt}) = T;
        continue;
    end

    % Output columns (create if missing)
    needcols = { ...
        'rf_area_deg2','pf_area_deg2', ...
        'rf_center_phi','rf_center_theta', ...
        'pf_center_phi','pf_center_theta', ...
        'rf_pf_center_dist_deg', ...
        'rf_mirror_pf_center_dist_deg', ...
        'rf_pf_center_dx_phi_deg', ...
        'rf_pf_center_dy_theta_deg', ...
        'rf_mirror_pf_center_dx_phi_deg', ...
        'rf_mirror_pf_center_dy_theta_deg' };
    for c = 1:numel(needcols)
        if ~ismember(needcols{c}, T.Properties.VariableNames)
            T.(needcols{c}) = nan(height(T),1);
        end
    end

    for i = 1:height(T)
        % ---- Get Hin/Hout (rebuild from in/out_tiles if missing) ----
        Hin = []; Hout = [];

        if ismember('Hin_norm',T.Properties.VariableNames) && ~isempty(T.Hin_norm{i})
            Hin = T.Hin_norm{i};
            mx  = max(Hin(:)); if mx>0 && mx~=1, Hin = Hin./mx; end
        else
            n_missing_Hin = n_missing_Hin + 1;
            if ismember('in_tiles',T.Properties.VariableNames) && ~isempty(T.in_tiles{i})
                Tin = T.in_tiles{i};  % columns: {'phi','theta','count'}
                Hin = local_weighted_hist2(Tin.theta, Tin.phi, Tin.count, theta_edges, phi_edges);
                mx  = max(Hin(:)); if mx>0, Hin = Hin./mx; end
            end
        end

        if ismember('Hout_norm',T.Properties.VariableNames) && ~isempty(T.Hout_norm{i})
            Hout = T.Hout_norm{i};
            mx   = max(Hout(:)); if mx>0 && mx~=1, Hout = Hout./mx; end
        else
            n_missing_Hout = n_missing_Hout + 1;
            if ismember('out_tiles',T.Properties.VariableNames) && ~isempty(T.out_tiles{i})
                Tout = T.out_tiles{i};
                Hout = local_weighted_hist2(Tout.theta, Tout.phi, Tout.count, theta_edges, phi_edges);
                mx   = max(Hout(:)); if mx>0, Hout = Hout./mx; end
            end
        end

        % Skip if no usable heatmap
        if isempty(Hin)  || ~any(Hin(:)>0),  continue; end
        if isempty(Hout) || ~any(Hout(:)>0), continue; end

        % ---- RF/PF metrics ----
        rf = compute_field_metrics(Hin);
        pf = compute_field_metrics(Hout);

        % Thresholded centers
        phi_rf = rf.center_phi_thr;    th_rf = rf.center_theta_thr;
        phi_pf = pf.center_phi_thr;    th_pf = pf.center_theta_thr;

        % Component & total distances
        dx  = circ_diff_deg(phi_rf, phi_pf);  % phi circular diff (deg)
        dy  = th_rf - th_pf;                  % theta linear diff (deg)
        dist_centers = hypot(dx, dy);

        % phi-mirrored RF center
        phi_rf_mirror = wrapTo180(-phi_rf);
        dx_m = circ_diff_deg(phi_rf_mirror, phi_pf);
        dy_m = dy;                            % theta unchanged
        dist_mirror = hypot(dx_m, dy_m);

        % ---- Store ----
        T.rf_area_deg2(i)  = rf.area_deg2_thr;
        T.pf_area_deg2(i)  = pf.area_deg2_thr;

        T.rf_center_phi(i)   = phi_rf;
        T.rf_center_theta(i) = th_rf;
        T.pf_center_phi(i)   = phi_pf;
        T.pf_center_theta(i) = th_pf;

        T.rf_pf_center_dist_deg(i)        = dist_centers;
        T.rf_mirror_pf_center_dist_deg(i) = dist_mirror;

        T.rf_pf_center_dx_phi_deg(i)          = dx;
        T.rf_pf_center_dy_theta_deg(i)        = dy;
        T.rf_mirror_pf_center_dx_phi_deg(i)   = dx_m;
        T.rf_mirror_pf_center_dy_theta_deg(i) = dy_m;
    end

    all_BLP_by_type.(type_names{tt}) = T;
end

% Metadata
all_BLP_by_type_meta.bin_phi_deg   = dphi;
all_BLP_by_type_meta.bin_theta_deg = dtheta;
all_BLP_by_type_meta.contour_level = contour_level;

save(matfile_out, 'all_BLP_by_type', 'all_BLP_by_type_meta', '-v7.3');

fprintf('Done.\nHin missing attempts (reconstructed from in_tiles): %d\nHout missing attempts (reconstructed from out_tiles): %d\n', ...
    n_missing_Hin, n_missing_Hout);

%% ===================== Local functions =====================
function out = local_field_metrics(Hnorm, phi_centers, theta_centers, level, dphi, dtheta, circ_mean_deg)
% Contour area & centers from a normalized heatmap (0..1).
% - area: #pixels(H>=level) * dphi*dtheta (deg^2)
% - thresholded center: H-weighted within the mask
% - all-bin center: H-weighted over all bins (fallback)
    out = struct('area_deg2_thr',0, ...
                 'center_phi_thr',nan,'center_theta_thr',nan, ...
                 'center_phi_all',nan,'center_theta_all',nan);

    if isempty(Hnorm) || ~any(Hnorm(:)>0)
        return;
    end

    [TH, PH] = ndgrid(theta_centers, phi_centers);  % TH (rows=theta), PH (cols=phi)

    % All-bin weighted center (reference / fallback)
    W_all = Hnorm;
    wsum_all = sum(W_all(:));
    if wsum_all > 0
        out.center_theta_all = sum(W_all(:).*TH(:)) / wsum_all;
        out.center_phi_all   = circ_mean_deg(PH(:), W_all(:));
    end

    % Threshold mask
    mask = Hnorm >= level;
    if ~any(mask(:))
        % No region above threshold -> fall back to the all-bin center
        out.area_deg2_thr    = 0;
        out.center_phi_thr   = out.center_phi_all;
        out.center_theta_thr = out.center_theta_all;
        return;
    end

    % Area (deg^2)
    out.area_deg2_thr = nnz(mask) * (dphi * dtheta);

    % Thresholded weighted center
    W = Hnorm;
    W(~mask) = 0;
    wsum = sum(W(:));
    if wsum > 0
        out.center_theta_thr = sum(W(:).*TH(:)) / wsum;
        out.center_phi_thr   = circ_mean_deg(PH(:), W(:));
    else
        out.center_theta_thr = out.center_theta_all;
        out.center_phi_thr   = out.center_phi_all;
    end
end

function H = local_weighted_hist2(theta_vals, phi_vals, weights, theta_edges, phi_edges)
% Weighted 2D histogram (rows = theta bins, cols = phi bins). Drops NaN/Inf;
% samples outside the bin range are ignored.
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

    n = min([numel(theta_vals), numel(phi_vals), numel(weights)]);
    theta_vals = theta_vals(1:n);
    phi_vals   = phi_vals(1:n);
    weights    = weights(1:n);

    valid = isfinite(theta_vals) & isfinite(phi_vals) & isfinite(weights);
    theta_vals = theta_vals(valid);
    phi_vals   = phi_vals(valid);
    weights    = weights(valid);

    [~,~,ib_theta] = histcounts(theta_vals, theta_edges);
    [~,~,ib_phi]   = histcounts(phi_vals,   phi_edges);

    good = (ib_theta>0) & (ib_phi>0);
    ib_theta = ib_theta(good);
    ib_phi   = ib_phi(good);
    weights  = weights(good);

    nRow = numel(theta_edges)-1; % theta bins -> rows
    nCol = numel(phi_edges)-1;   % phi bins   -> cols
    H = accumarray([ib_theta, ib_phi], weights, [nRow, nCol], @sum, 0);
end
