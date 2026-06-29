%% 1. Load data and initialize
% MCNS analogue of the FAFB Figures/fig_6F_H_J_BLP_area_center_plot.m.
% Per-type RF/PF area & center metrics from Data_Processing/s18_BLP_RFs_PFs_area_center.m
close all; clc; clearvars;

baseDir = fileparts(fileparts(mfilename('fullpath')));  % resolved from this script's location (portable; no absolute paths)
load(fullfile(baseDir, 'Processed_Data', 'BLP_RFs_PFs_area_center.mat'))   % all_BLP_by_type
neuron_types = fieldnames(all_BLP_by_type);

%% 2. Per-type grouped-bar metrics (figures 1 / 2 / 3)
% Each config compares two metrics per neuron type, sorted by the first metric (desc).
configs = {
    struct('name','Mi1 Area (In vs Out)', ...
           'fn_pair', {{ @(T) T.rf_area_deg2, @(T) T.pf_area_deg2 }}, ...
           'labels',  {{ 'Mi1 in area', 'Mi1 out area' }}, ...
           'ylabel','Area', 'ylim',[], 'sort_by','first', ...
           'savename','fig_6_BLP_area_bar'), ...

    struct('name','Center diff (In - Out)', ...
           'fn_pair', {{ @(T) T.rf_pf_center_dx_phi_deg, @(T) T.rf_pf_center_dy_theta_deg }}, ...
           'labels',  {{ 'dx raw', 'dy raw' }}, ...
           'ylabel','Degrees', 'ylim',[], 'sort_by','first', ...
           'savename','fig_6_BLP_center_diff'), ...

    struct('name','Center diff flip (In - Out)', ...
           'fn_pair', {{ @(T) T.rf_mirror_pf_center_dx_phi_deg, @(T) T.rf_mirror_pf_center_dy_theta_deg }}, ...
           'labels',  {{ 'dx flip', 'dy flip' }}, ...
           'ylabel','Degrees', 'ylim',[], 'sort_by','first', ...
           'savename','fig_6_BLP_center_diff_flip')
};

for c = 1:numel(configs)
    cfg = configs{c};

    fnA = cfg.fn_pair{1};
    fnB = cfg.fn_pair{2};
    [types, muA, sdA] = compute_stats(all_BLP_by_type, neuron_types, fnA);
    [~,     muB, sdB] = compute_stats(all_BLP_by_type, neuron_types, fnB);

    switch cfg.sort_by
        case 'avg', base = (muA + muB) / 2;
        otherwise,  base = muA;   % 'first'
    end
    [types_s, ~, sdA_s, order] = sort_desc(types, base, sdA);
    muA_s = muA(order); muB_s = muB(order); sdB_s = sdB(order);

    M = [muA_s muB_s]; Sd = [sdA_s sdB_s];

    figure('Name', cfg.name, 'Color', 'w');
    bar(M, 'grouped'); hold on;
    ngroups = numel(types_s); nbars = 2;
    groupwidth = min(0.8, nbars/(nbars + 1.5));
    for iBar = 1:nbars
        x = (1:ngroups) - groupwidth/2 + (2*iBar-1)*groupwidth/(2*nbars);
        errorbar(x, M(:,iBar), Sd(:,iBar), '.k', 'LineWidth', 1.2);
    end
    hold off;

    set(gca, 'XTick', 1:ngroups, 'XTickLabel', sanitize_labels(types_s), ...
             'XTickLabelRotation', 45, 'Box','off', 'TickDir','out');
    legend(cfg.labels, 'Location','best');
    ylabel(cfg.ylabel);
    title(['Mean \pm STD of ', cfg.name, ' per Neuron Type']);
    if ~isempty(cfg.ylim), ylim(cfg.ylim); end

end

%% 3. Figure 4: RF area vs PF area scatter (mean per type)
figure('Name','Mi1 Area: RF (In) vs PF (Out) (Mean per Type)', 'Color','w');
ax = gca; hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on'); set(ax, 'TickDir', 'out');

all_X_data = []; all_Y_data = [];        % all individual neurons
type_means_X = []; type_means_Y = [];    % per-type means
plot_handles = []; plot_labels = {};

for i = 1:numel(neuron_types)
    T = all_BLP_by_type.(neuron_types{i});
    valsX = T.rf_area_deg2(:);
    valsY = T.pf_area_deg2(:);
    valid = ~isnan(valsX) & ~isnan(valsY);
    valsX = valsX(valid); valsY = valsY(valid);
    if isempty(valsX), continue; end

    all_X_data = [all_X_data; valsX];
    all_Y_data = [all_Y_data; valsY];

    mX = mean(valsX); mY = mean(valsY);
    type_means_X = [type_means_X; mX];
    type_means_Y = [type_means_Y; mY];

    h = scatter(ax, mX, mY, 36, 'filled', 'MarkerFaceAlpha', 0.6);
    plot_handles(end+1) = h;
    plot_labels{end+1}  = strrep(neuron_types{i}, '_', '\_');
end

xlabel(ax, 'RF Area (deg^2) - In');
ylabel(ax, 'PF Area (deg^2) - Out');
if ~isempty(plot_handles)
    legend(plot_handles, plot_labels, 'Location', 'eastoutside', 'Interpreter', 'tex');
end

% 1:1 line and equal axes
axis(ax, 'equal');
xlims = xlim(ax); ylims = ylim(ax);
max_lim = max(xlims(2), ylims(2));
line([0 max_lim], [0 max_lim], 'Color', 'k', 'LineStyle', '--', 'Parent', ax);
xlim(ax, [0 max_lim]); ylim(ax, [0 max_lim]);

% Correlation (individual neurons & per-type means)
[rho_p_all,  p_p_all]  = corr(all_X_data, all_Y_data, 'Type', 'Pearson');
[rho_s_all,  p_s_all]  = corr(all_X_data, all_Y_data, 'Type', 'Spearman');
[rho_p_mean, p_p_mean] = corr(type_means_X, type_means_Y, 'Type', 'Pearson');
[rho_s_mean, p_s_mean] = corr(type_means_X, type_means_Y, 'Type', 'Spearman');

fprintf('\n--- Correlation Analysis Results ---\n');
fprintf('[Individual Neurons (N=%d)]\n', numel(all_X_data));
fprintf('  Pearson:  r = %.3f, p = %.4f\n', rho_p_all, p_p_all);
fprintf('  Spearman: rho = %.3f, p = %.4f\n', rho_s_all, p_s_all);
fprintf('\n[Type Means (N=%d)]\n', numel(type_means_X));
fprintf('  Pearson:  r = %.3f, p = %.4f\n', rho_p_mean, p_p_mean);
fprintf('  Spearman: rho = %.3f, p = %.4f\n', rho_s_mean, p_s_mean);
fprintf('------------------------------------\n');

title({'Mi1 Area: RF (In) vs PF (Out) (Mean per Type)', ...
       sprintf('Pearson r=%.2f, Spearman \\rho=%.2f', rho_p_all, rho_s_all)}, 'FontSize', 10);

%% 4. Per-type mean RF / PF area table (console)
n_types = numel(neuron_types);
mean_rf_areas = nan(n_types, 1);
mean_pf_areas = nan(n_types, 1);
for i = 1:n_types
    T = all_BLP_by_type.(neuron_types{i});
    mean_rf_areas(i) = mean(T.rf_area_deg2, 'omitnan');
    mean_pf_areas(i) = mean(T.pf_area_deg2, 'omitnan');
end
mean_table = table(string(neuron_types), mean_rf_areas, mean_pf_areas, ...
    'VariableNames', {'Type', 'Mean_RF_Area_deg2', 'Mean_PF_Area_deg2'});
disp(mean_table);

%% ===================== Local functions =====================
function labels = sanitize_labels(types)
labels = cellfun(@(s) strrep(s,'_','\_'), types, 'UniformOutput', false);
end

function [types, mu, sd] = compute_stats(all_BLP_by_type, neuron_types, fn_handle)
n = numel(neuron_types);
mu = nan(n,1); sd = nan(n,1);
for i = 1:n
    T = all_BLP_by_type.(neuron_types{i});
    vals = fn_handle(T);
    vals = vals(~isnan(vals));
    if ~isempty(vals)
        mu(i) = mean(vals); sd(i) = std(vals);
    end
end
types = neuron_types(:);
end

function [types_s, base_s, sd_s, order] = sort_desc(types, base, sd)
base_for_sort = base;
base_for_sort(isnan(base_for_sort)) = -inf;   % push NaNs to the end
[~, order] = sort(base_for_sort, 'descend');
types_s = types(order);
base_s  = base(order);
sd_s    = sd(order);
end
