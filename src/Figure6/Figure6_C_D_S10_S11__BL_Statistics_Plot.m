%% ========================================================================
%  Neuron-type metrics plotting (sorted by descending mean)
%  - 단일 지표: fn
%  - 두 지표 비교(그룹 막대): fn_pair  (⚠️ 중괄호 두 겹으로 감싸기!!)
% ========================================================================
close all; clc; clearvars;

% 데이터 로드 (변수 all_neurons_by_type 포함)
load('NeuronAnalysis_byType_withMetrics.mat');

% ------------------------------------------------------------------------
% [설정] 그리고 싶은 지표들 정의
% ------------------------------------------------------------------------
configs = {
    % ⚠️ fn_pair, labels 는 중괄호 두 겹으로!
    struct('name','Mi1 Area (In vs Out)', ...
           'fn_pair', {{ @(T) T.rf_area_deg2, @(T) T.pf_area_deg2 }}, ...
           'labels',  {{ 'Mi1 in area', 'Mi1 out area' }}, ...
           'ylabel','Area', 'ylim',[], 'sort_by','first'), ...  % 콤마!

    struct('name','Center diff (In - Out)', ...
           'fn_pair', {{ @(T) T.rf_pf_center_dx_phi_deg, @(T) T.rf_pf_center_dy_theta_deg }}, ...
           'labels',  {{ 'dx raw', 'dy raw' }}, ...
           'ylabel','Degrees', 'ylim',[], 'sort_by','first'), ...

    struct('name','Center diff flip(In - Out)', ...
           'fn_pair', {{ @(T) T.rf_mirror_pf_center_dx_phi_deg, @(T) T.rf_mirror_pf_center_dy_theta_deg }}, ...
           'labels',  {{ 'dx flip', 'dy flip' }}, ...
           'ylabel','Degrees', 'ylim',[], 'sort_by','first')
};

% ------------------------------------------------------------------------
% 실행
% ------------------------------------------------------------------------
neuron_types = fieldnames(all_neurons_by_type);

% 방어: configs 는 cell 이어야 함
if ~iscell(configs)
    error('configs must be a cell array where each cell contains a 1x1 struct.');
end

for c = 1:numel(configs)
    cfg = configs{c};

    % 혹시 struct 배열이면 첫 요소만 사용 (스칼라화)
    if numel(cfg) ~= 1
        warning('configs{%d} is a struct array (%dx%d). Using the first element only.', ...
                 c, size(cfg,1), size(cfg,2));
        cfg = cfg(1);
    end

    % ------------------------ 단일 지표 ------------------------
    if isfield(cfg,'fn') && ~isempty(cfg.fn)
        [types, mu, sd] = compute_stats(all_neurons_by_type, neuron_types, cfg.fn);
        [types_s, mu_s, sd_s] = sort_desc(types, mu, sd);

        figure('Name', cfg.name);
        bar(mu_s, 'FaceColor', [0.2 0.6 0.8]); hold on;
        errorbar(1:numel(mu_s), mu_s, sd_s, '.k', 'LineWidth', 1.2);
        hold off;

        set(gca, 'XTick', 1:numel(types_s), ...
                 'XTickLabel', sanitize_labels(types_s), ...
                 'XTickLabelRotation', 45, 'Box','off', 'TickDir','out');
        ylabel(cfg.ylabel);
        title(['Mean \pm STD of ', cfg.name, ' per Neuron Type']);
        yl = getfield_def(cfg,'ylim',[]);
        if ~isempty(yl), ylim(yl); end

    % ---------------------- 두 지표 비교 -----------------------
    elseif isfield(cfg,'fn_pair') && ~isempty(cfg.fn_pair)
        fp = cfg.fn_pair;
        if ~iscell(fp) || numel(fp) < 2
            error('configs{%d}.fn_pair must be a cell with at least two function handles.', c);
        end
        fnA = fp{1};
        fnB = fp{2};

        [types, muA, sdA] = compute_stats(all_neurons_by_type, neuron_types, fnA);
        [~,     muB, sdB] = compute_stats(all_neurons_by_type, neuron_types, fnB);

        sort_by = getfield_def(cfg,'sort_by','first');
        switch sort_by
            case 'avg'
                base = (muA + muB) / 2;
            otherwise % 'first'
                base = muA;
        end

        [types_s, ~, sdA_s, order] = sort_desc(types, base, sdA);
        muA_s = muA(order); muB_s = muB(order);
        sdB_s = sdB(order);

        M = [muA_s muB_s]; S = [sdA_s sdB_s];

        figure('Name', cfg.name);
        bar(M, 'grouped'); hold on;

        ngroups = numel(types_s); nbars = 2;
        groupwidth = min(0.8, nbars/(nbars + 1.5));
        for iBar = 1:nbars
            x = (1:ngroups) - groupwidth/2 + (2*iBar-1)*groupwidth/(2*nbars);
            errorbar(x, M(:,iBar), S(:,iBar), '.k', 'LineWidth', 1.2);
        end
        hold off;

        set(gca, 'XTick', 1:ngroups, ...
                 'XTickLabel', sanitize_labels(types_s), ...
                 'XTickLabelRotation', 45, 'Box','off', 'TickDir','out');
        legends = getfield_def(cfg,'labels',{'Series A','Series B'});
        legend(legends, 'Location','best');
        ylabel(cfg.ylabel);
        title(['Mean \pm STD of ', cfg.name, ' per Neuron Type']);
        yl = getfield_def(cfg,'ylim',[]);
        if ~isempty(yl), ylim(yl); end

    else
        warning('Config %d (%s): neither fn nor valid fn_pair found. Skipped.', c, cfg.name);
    end
end

% ========================================================================
% 헬퍼 함수들
% ========================================================================
function v = getfield_def(S, field, default)
    if isfield(S, field) && ~isempty(S.(field))
        v = S.(field);
    else
        v = default;
    end
end

function labels = sanitize_labels(types)
    labels = cellfun(@(s) strrep(s,'_','\_'), types, 'UniformOutput', false);
end

function [types, mu, sd] = compute_stats(all_neurons_by_type, neuron_types, fn_handle)
    n = numel(neuron_types);
    mu = nan(n,1); sd = nan(n,1);
    for i = 1:n
        T = all_neurons_by_type.(neuron_types{i});
        vals = fn_handle(T);
        vals = vals(~isnan(vals));
        if isempty(vals)
            mu(i) = NaN; sd(i) = NaN;
        else
            mu(i) = mean(vals); sd(i) = std(vals);
        end
    end
    types = neuron_types(:);
end

function [types_s, base_s, sd_s, order] = sort_desc(types, base, sd)
    base_for_sort = base;
    base_for_sort(isnan(base_for_sort)) = -inf;   % NaN은 뒤로
    [~, order] = sort(base_for_sort, 'descend');
    types_s = types(order);
    base_s  = base(order);
    sd_s    = sd(order);
end
