clear all; close all; clc;

baseDir = fileparts(fileparts(mfilename('fullpath')));  % resolved from this script's location (portable; no absolute paths)

% FBP output-synapse depth profiles (FBP_Synapse_*), saved by
% Data_Processing/s11_FBP_output_layer.m
load(fullfile(baseDir, 'Processed_Data', 'FBP_output_synapse_layer.mat'));

% FBP upstream-synapse depth profiles (FBP_Upstream_Synapse_*), saved by
% Data_Processing/s12_FBP_upstream_layer.m
load(fullfile(baseDir, 'Processed_Data', 'FBP_upstream_synapse_layer.mat'));

% FFP / FBP / BDP neuron classification (provides RightFBP_NPIs), saved by
% Figures/fig_1D_E_postPI_prePI_right_neurons.m
load(fullfile(baseDir, 'Processed_Data', 'right_neurons_thr0.mat'),'RightFBP_NPIs')

%% Group the FBP neurons by cell type
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

%% Fill missing
% NaN means there is no visual input at all; 0 means one of the two sides has input.
FBP_Synapse_ME_R = fillmissing(FBP_Synapse_ME_R,'constant',0);
FBP_Synapse_ME_L = fillmissing(FBP_Synapse_ME_L ,'constant',0);
FBP_Upstream_Synapse_ME_R = fillmissing(FBP_Upstream_Synapse_ME_R,'constant',0);
FBP_Upstream_Synapse_ME_L = fillmissing(FBP_Upstream_Synapse_ME_L,'constant',0);

FBP_Synapse_LO_R = fillmissing(FBP_Synapse_LO_R,'constant',0);
FBP_Synapse_LO_L = fillmissing(FBP_Synapse_LO_L ,'constant',0);
FBP_Upstream_Synapse_LO_R = fillmissing(FBP_Upstream_Synapse_LO_R,'constant',0);
FBP_Upstream_Synapse_LO_L = fillmissing(FBP_Upstream_Synapse_LO_L,'constant',0);

FBP_Synapse_LOP_R = fillmissing(FBP_Synapse_LOP_R,'constant',0);
FBP_Synapse_LOP_L = fillmissing(FBP_Synapse_LOP_L ,'constant',0);
FBP_Upstream_Synapse_LOP_R = fillmissing(FBP_Upstream_Synapse_LOP_R,'constant',0);
FBP_Upstream_Synapse_LOP_L = fillmissing(FBP_Upstream_Synapse_LOP_L,'constant',0);

%% Gaussian kernel
sigma = 2;
x = -round(3*sigma):round(3*sigma);
gaussKernel = exp(-(x.^2)/(2*sigma^2));
gaussKernel = gaussKernel / sum(gaussKernel);  % normalize
% sigma = 2 -> kernel length = 13

%% Gaussian smoothing + overlap computation
% X = smoothed FBP output depth profile, Y_L / Y_R = smoothed upstream input depth
% profiles (contralateral L / ipsilateral R), each scaled by its share of the total
% upstream synapses. overlap = sum of the per-depth minimum of the two distributions.

% LOP
for i = 1:size(FBP_Synapse_LOP_R,2)
    X = FBP_Synapse_LOP_R(:,i);
    X = conv(X, gaussKernel, 'same');
    X = X ./ sum(X);
    X=fillmissing(X,'constant',0);

    Y_L =  FBP_Upstream_Synapse_LOP_L(:,i);
    Y_L = conv(Y_L, gaussKernel, 'same');
    Y_L = Y_L ./ sum(Y_L)*(sum(FBP_Upstream_Synapse_LOP_L(:,i)))/(sum(FBP_Upstream_Synapse_LOP_L(:,i))+sum(FBP_Upstream_Synapse_LOP_R(:,i)));
    Y_L=fillmissing(Y_L,'constant',0);
    Y_R =  FBP_Upstream_Synapse_LOP_R(:,i);
    Y_R = conv(Y_R, gaussKernel, 'same');
    Y_R = Y_R ./ sum(Y_R)*(sum(FBP_Upstream_Synapse_LOP_R(:,i)))/(sum(FBP_Upstream_Synapse_LOP_L(:,i))+sum(FBP_Upstream_Synapse_LOP_R(:,i)));
    Y_R=fillmissing(Y_R,'constant',0);

    overlap_idx_LOP(i) = sum(min(X, (Y_R+Y_L)));
    overlap_idx_LOP_L(i) = sum(min(X, (Y_L)));
    overlap_idx_LOP_R(i) = sum(min(X, (Y_R)));
end

% LO
for i = 1:size(FBP_Synapse_LO_R,2)
    X = FBP_Synapse_LO_R(:,i);
    X = conv(X, gaussKernel, 'same');
    X = X ./ sum(X);
    X=fillmissing(X,'constant',0);

    Y_L =  FBP_Upstream_Synapse_LO_L(:,i);
    Y_L = conv(Y_L, gaussKernel, 'same');
    Y_L = Y_L ./ sum(Y_L)*(sum(FBP_Upstream_Synapse_LO_L(:,i)))/(sum(FBP_Upstream_Synapse_LO_L(:,i))+sum(FBP_Upstream_Synapse_LO_R(:,i)));
    Y_L=fillmissing(Y_L,'constant',0);

    Y_R =  FBP_Upstream_Synapse_LO_R(:,i);
    Y_R = conv(Y_R, gaussKernel, 'same');
    Y_R = Y_R ./ sum(Y_R)*(sum(FBP_Upstream_Synapse_LO_R(:,i)))/(sum(FBP_Upstream_Synapse_LO_L(:,i))+sum(FBP_Upstream_Synapse_LO_R(:,i)));
    Y_R=fillmissing(Y_R,'constant',0);

    overlap_idx_LO(i) = sum(min(X, (Y_R+Y_L)));
    overlap_idx_LO_L(i) = sum(min(X, (Y_L)));
    overlap_idx_LO_R(i) = sum(min(X, (Y_R)));
end

% ME
for i = 1:size(FBP_Synapse_ME_R,2)
    X = FBP_Synapse_ME_R(:,i);
    X = conv(X, gaussKernel, 'same');
    X = X ./ sum(X);
    X=fillmissing(X,'constant',0);

    Y_L =  FBP_Upstream_Synapse_ME_L(:,i);
    Y_L = conv(Y_L, gaussKernel, 'same');
    Y_L = Y_L ./ sum(Y_L)*(sum(FBP_Upstream_Synapse_ME_L(:,i)))/(sum(FBP_Upstream_Synapse_ME_L(:,i))+sum(FBP_Upstream_Synapse_ME_R(:,i)));
    Y_L=fillmissing(Y_L,'constant',0);

    Y_R =  FBP_Upstream_Synapse_ME_R(:,i);
    Y_R = conv(Y_R, gaussKernel, 'same');
    Y_R = Y_R ./ sum(Y_R)*(sum(FBP_Upstream_Synapse_ME_R(:,i)))/(sum(FBP_Upstream_Synapse_ME_L(:,i))+sum(FBP_Upstream_Synapse_ME_R(:,i)));
    Y_R=fillmissing(Y_R,'constant',0);

    overlap_idx_ME(i) = sum(min(X, (Y_R+Y_L)));
    overlap_idx_ME_L(i) = sum(min(X, (Y_L)));
    overlap_idx_ME_R(i) = sum(min(X, (Y_R)));
end

%% Figure 1 (panel S2E): per-neuron overlap, ME
figure(1); set(gcf,'Color','w'); hold on;
y = [overlap_idx_ME_L' overlap_idx_ME_R'];

b = bar(FBP_Me.type, y);
b(1).FaceColor = [0.4660 0.6740 0.1880];  % column 1: upstream input contralateral (L)
b(2).FaceColor = [0 0.4470 0.7410];       % column 2: upstream input ipsilateral (R)

set(gca, 'TickDir', 'out', 'Box', 'off');
ylim([0 1]);
ylabel('Overlap coefficient (Gaussian smoothed)');
title('Overlap of input and output depth distributions (with Gaussian smoothing)');

%% Figure 2 (panel S2E): per-neuron overlap, LO
figure(2); set(gcf,'Color','w'); hold on;
y = [overlap_idx_LO_L' overlap_idx_LO_R'];

b = bar(FBP_Lo.type, y);
b(1).FaceColor = [0.4660 0.6740 0.1880];  % column 1: upstream input contralateral (L)
b(2).FaceColor = [0 0.4470 0.7410];       % column 2: upstream input ipsilateral (R)

set(gca, 'TickDir', 'out', 'Box', 'off');
ylim([0 1]);
ylabel('Overlap coefficient (Gaussian smoothed)');
title('Overlap of input and output depth distributions (with Gaussian smoothing)');

%% Figure 3 (panel S2E): per-neuron overlap, LOP
figure(3); set(gcf,'Color','w'); hold on;
y = [overlap_idx_LOP_L' overlap_idx_LOP_R'];

b = bar(FBP_Lop.type, y);
b(1).FaceColor = [0.4660 0.6740 0.1880];  % column 1: upstream input contralateral (L)
b(2).FaceColor = [0 0.4470 0.7410];       % column 2: upstream input ipsilateral (R)

set(gca, 'TickDir', 'out', 'Box', 'off');
ylim([0 1]);
ylabel('Overlap coefficient (Gaussian smoothed)');
title('Overlap of input and output depth distributions (with Gaussian smoothing)');

%% Overall mean and std across all neurons
all_overlap = [overlap_idx_ME, overlap_idx_LO, overlap_idx_LOP];
all_overlap_L = [overlap_idx_ME_L, overlap_idx_LO_L, overlap_idx_LOP_L];
all_overlap_R = [overlap_idx_ME_R, overlap_idx_LO_R, overlap_idx_LOP_R];

mean_overlap_all = mean(all_overlap, 'omitnan');
mean_overlap_all_L = mean(all_overlap_L, 'omitnan');
mean_overlap_all_R = mean(all_overlap_R, 'omitnan');

std_overlap_all = std(all_overlap, 0, 'omitnan');
std_overlap_all_L = std(all_overlap_L, 0, 'omitnan');
std_overlap_all_R = std(all_overlap_R, 0, 'omitnan');

fprintf('Overlap mean and std across all neurons:\n');
fprintf('- All  : mean = %.4f, std = %.4f\n', mean_overlap_all, std_overlap_all);
fprintf('- Left : mean = %.4f, std = %.4f\n', mean_overlap_all_L, std_overlap_all_L);
fprintf('- Right: mean = %.4f, std = %.4f\n', mean_overlap_all_R, std_overlap_all_R);

%% Figure 4 (panel 4G): overlap distribution across all neurons, per region/side
data = {overlap_idx_ME_L(:), overlap_idx_ME_R(:), ...
        overlap_idx_LO_L(:), overlap_idx_LO_R(:), ...
        overlap_idx_LOP_L(:), overlap_idx_LOP_R(:)};
group = repelem(1:length(data), cellfun(@length, data));   % group labels
combinedData = vertcat(data{:});                            % concatenate into one array

figure(4);set(gcf,'Color','w')
boxplot(combinedData, group, 'Labels', {'Me L', 'Me R', 'Lo L','Lo R','LoP L','LoP R'},'Notch','on','Symbol','');
set(gca,'TickDir','out','box','off')
ylim([0 0.8])
hold on;

%% Paired comparison: Wilcoxon signed-rank test (signrank)
% Within each target neuropil, contralateral (L) vs ipsilateral (R) overlap is
% paired by neuron (same column index), so a paired test is used.
% Boxplot group order: Me L/R, Lo L/R, LoP L/R.
np_names = {'ME', 'LO', 'LOP'};
np_L = {overlap_idx_ME_L(:), overlap_idx_LO_L(:), overlap_idx_LOP_L(:)};
np_R = {overlap_idx_ME_R(:), overlap_idx_LO_R(:), overlap_idx_LOP_R(:)};
nComp = numel(np_names);            % number of comparisons (Bonferroni, n = 3)
alpha = 0.05 / nComp;               % Bonferroni-corrected significance level alpha = 0.05/3 ~ 0.0167
% Asterisk thresholds (Bonferroni-corrected p):
%   *   p < alpha   = 0.05/n  ~ 0.0167
%   **  p < 0.01/n           ~ 0.0033
%   *** p < 0.001/n          ~ 0.00033
p_signrank = nan(1, nComp);         % raw p-values

fprintf('\n[Wilcoxon signed-rank test (signrank, paired): Left vs Right, Bonferroni n=%d, alpha=%.4f]\n', nComp, alpha);
for k = 1:numel(np_names)
    L = np_L{k};
    R = np_R{k};
    valid = ~isnan(L) & ~isnan(R);   % drop NaNs pairwise
    L = L(valid);
    R = R(valid);
    [p, ~, stats] = signrank(L, R);
    p_signrank(k) = p;
    if isfield(stats, 'zval')
        fprintf('- %-3s : p = %.4g (n = %d, signedrank = %g, z = %.3f)\n', ...
            np_names{k}, p, numel(L), stats.signedrank, stats.zval);
    else
        fprintf('- %-3s : p = %.4g (n = %d, signedrank = %g)\n', ...
            np_names{k}, p, numel(L), stats.signedrank);
    end

    % Significance marker (raw p vs Bonferroni-corrected thresholds), above the two boxes
    x1 = 2*k - 1;   % Left
    x2 = 2*k;       % Right
    yMax = max([L; R]);
    yBar = yMax + 0.05;
    plot([x1 x1 x2 x2], [yBar-0.01 yBar yBar yBar-0.01], 'k', 'LineWidth', 1);
    if p < 0.001/nComp
        star = '***';
    elseif p < 0.01/nComp
        star = '**';
    elseif p < alpha          % alpha = 0.05/nComp
        star = '*';
    else
        star = 'n.s.';
    end
    text((x1+x2)/2, yBar + 0.01, star, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom');
end
