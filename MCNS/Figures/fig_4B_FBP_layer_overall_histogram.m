%% fig_4B_FBP_layer_overall_histogram (MCNS)
% MCNS analogue of the FAFB Figures/fig_4F_FBP_layer_overall_histogram.m.
% Loads the FBP output- and upstream-synapse depth profiles saved by s11/s12 and,
% for each optic lobe (ME / LO / LOP), overlays the mean +/- 1 SD depth profile of
% the FBP output (orange), upstream input ipsilateral (blue), and upstream input
% contralateral (green).
clear all; close all; clc;

%% ===== Data load =====
baseDir = fileparts(fileparts(mfilename('fullpath')));  % resolved from this script's location (portable; no absolute paths)

% FBP output-synapse depth profiles (FBP_Synapse_*), saved by
% Data_Processing/s11_FBP_output_layer.m
load(fullfile(baseDir, 'Processed_Data', 'FBP_output_synapse_layer.mat'));

% FBP upstream-synapse depth profiles (FBP_Upstream_Synapse_*), saved by
% Data_Processing/s12_FBP_upstream_layer.m
load(fullfile(baseDir, 'Processed_Data', 'FBP_upstream_synapse_layer.mat'));

%% ===== Axis labels (depth ticks) =====
edges = -7e-5:1e-6:5e-5;

%% ===== Fill missing =====
% NaN = no visual input; replace with 0 as instructed.
FBP_Synapse_ME_R           = fillmissing(FBP_Synapse_ME_R,'constant',0);
FBP_Synapse_ME_L           = fillmissing(FBP_Synapse_ME_L,'constant',0);
FBP_Upstream_Synapse_ME_R  = fillmissing(FBP_Upstream_Synapse_ME_R,'constant',0);
FBP_Upstream_Synapse_ME_L  = fillmissing(FBP_Upstream_Synapse_ME_L,'constant',0);

FBP_Synapse_LO_R           = fillmissing(FBP_Synapse_LO_R,'constant',0);
FBP_Synapse_LO_L           = fillmissing(FBP_Synapse_LO_L,'constant',0);
FBP_Upstream_Synapse_LO_R  = fillmissing(FBP_Upstream_Synapse_LO_R,'constant',0);
FBP_Upstream_Synapse_LO_L  = fillmissing(FBP_Upstream_Synapse_LO_L,'constant',0);

FBP_Synapse_LOP_R          = fillmissing(FBP_Synapse_LOP_R,'constant',0);
FBP_Synapse_LOP_L          = fillmissing(FBP_Synapse_LOP_L,'constant',0);
FBP_Upstream_Synapse_LOP_R = fillmissing(FBP_Upstream_Synapse_LOP_R,'constant',0);
FBP_Upstream_Synapse_LOP_L = fillmissing(FBP_Upstream_Synapse_LOP_L,'constant',0);

%% ===== Colors =====
col_orange = [0.8500 0.3250 0.0980];   % FBP output (R)
col_blue   = [0      0.4470 0.7410];   % upstream input, ipsilateral (R)
col_green  = [0.4660 0.6740 0.1880];   % upstream input, contralateral (L)

%% ======================= Figure 1: ME =======================
figure(1); clf; set(gcf,'Color','w'); hold on;

y_ME = (1:size(FBP_Synapse_ME_R,1))';

% Means
m_syn_ME_R  = mean(FBP_Synapse_ME_R,  2);                            % orange
m_up_ME_R   = mean(FBP_Upstream_Synapse_ME_R, 2, 'omitnan');         % blue
m_up_ME_L   = mean(FBP_Upstream_Synapse_ME_L, 2, 'omitnan');         % green

% SDs
s_syn_ME_R  = std(FBP_Synapse_ME_R,  0, 2);                          % orange
s_up_ME_R   = std(FBP_Upstream_Synapse_ME_R, 0, 2, 'omitnan');       % blue
s_up_ME_L   = std(FBP_Upstream_Synapse_ME_L, 0, 2, 'omitnan');       % green

% Shaded ±1SD bands
shaded_band(m_syn_ME_R, s_syn_ME_R, y_ME, col_orange);
shaded_band(m_up_ME_R,  s_up_ME_R,  y_ME, col_blue);
shaded_band(m_up_ME_L,  s_up_ME_L,  y_ME, col_green);

% Mean lines
plot(m_syn_ME_R, y_ME, 'Color', col_orange, 'LineWidth', 2);
plot(m_up_ME_R,  y_ME, 'Color', col_blue,   'LineWidth', 2);
plot(m_up_ME_L,  y_ME, 'Color', col_green,  'LineWidth', 2);

xlabel('Normalized value');
ylabel('Innvervation depth');
title('ME innervation depth total');
set(gca,'TickDir','out','Box','off');
ylim([0.5 90.5]);
xlim([0 1]);
set(gca,'YTick',0.5:10:90+0.5,'YTickLabel',edges(1:10:end));

%% ======================= Figure 2: LO =======================
figure(2); clf; set(gcf,'Color','w'); hold on;

y_LO = (1:size(FBP_Synapse_LO_R,1))';

% Means
m_syn_LO_R  = mean(FBP_Synapse_LO_R,  2);                            % orange
m_up_LO_R   = mean(FBP_Upstream_Synapse_LO_R, 2, 'omitnan');         % blue
m_up_LO_L   = mean(FBP_Upstream_Synapse_LO_L, 2, 'omitnan');         % green

% SDs
s_syn_LO_R  = std(FBP_Synapse_LO_R,  0, 2);                          % orange
s_up_LO_R   = std(FBP_Upstream_Synapse_LO_R, 0, 2, 'omitnan');       % blue
s_up_LO_L   = std(FBP_Upstream_Synapse_LO_L, 0, 2, 'omitnan');       % green

% Shaded ±1SD bands
shaded_band(m_syn_LO_R, s_syn_LO_R, y_LO, col_orange);
shaded_band(m_up_LO_R,  s_up_LO_R,  y_LO, col_blue);
shaded_band(m_up_LO_L,  s_up_LO_L,  y_LO, col_green);

% Mean lines
plot(m_syn_LO_R, y_LO, 'Color', col_orange, 'LineWidth', 2);
plot(m_up_LO_R,  y_LO, 'Color', col_blue,   'LineWidth', 2);
plot(m_up_LO_L,  y_LO, 'Color', col_green,  'LineWidth', 2);

xlabel('Normalized value');
ylabel('Innvervation depth');
title('LO innervation depth total');
set(gca,'TickDir','out','Box','off');
set(gca,'YTick',0.5:10:100+0.5,'YTickLabel',edges(1:10:end));
ylim([15.5 90.5]);
xlim([0 1]);

%% ======================= Figure 3: LOP =======================
figure(3); clf; set(gcf,'Color','w'); hold on;

y_LOP = (1:size(FBP_Synapse_LOP_R,1))';

% Means
m_syn_LOP_R = mean(FBP_Synapse_LOP_R, 2);                            % orange
m_up_LOP_R  = mean(FBP_Upstream_Synapse_LOP_R, 2, 'omitnan');        % blue
m_up_LOP_L  = mean(FBP_Upstream_Synapse_LOP_L, 2, 'omitnan');        % green

% SDs
s_syn_LOP_R = std(FBP_Synapse_LOP_R, 0, 2);                          % orange
s_up_LOP_R  = std(FBP_Upstream_Synapse_LOP_R, 0, 2, 'omitnan');      % blue
s_up_LOP_L  = std(FBP_Upstream_Synapse_LOP_L, 0, 2, 'omitnan');      % green

% Shaded ±1SD bands
shaded_band(m_syn_LOP_R, s_syn_LOP_R, y_LOP, col_orange);
shaded_band(m_up_LOP_R,  s_up_LOP_R,  y_LOP, col_blue);
shaded_band(m_up_LOP_L,  s_up_LOP_L,  y_LOP, col_green);

% Mean lines
plot(m_syn_LOP_R, y_LOP, 'Color', col_orange, 'LineWidth', 2);
plot(m_up_LOP_R,  y_LOP, 'Color', col_blue,   'LineWidth', 2);
plot(m_up_LOP_L,  y_LOP, 'Color', col_green,  'LineWidth', 2);

set(gca,'YDir','reverse');  % top -> bottom
xlabel('Normalized value');
ylabel('Innvervation depth');
title('LOP innervation depth total');
set(gca,'TickDir','out','Box','off');
set(gca,'YTick',0.5:10:100+0.5,'YTickLabel',edges(1:10:end));
ylim([60.5 100.5]);
xlim([0 1]);

%% ===== Local function: draw mean±SD shaded band =====
function h = shaded_band(x_mean, x_sd, yvec, rgb)
    % x_mean, x_sd : (N x 1) vectors along depth
    % yvec         : (N x 1) depth index
    % rgb          : 1x3 color (same as line color)
    xl = x_mean - x_sd;
    xu = x_mean + x_sd;

    ok = ~(isnan(xl) | isnan(xu) | isnan(yvec));
    xl = xl(ok); xu = xu(ok); y = yvec(ok);

    h = fill([xl; flipud(xu)], [y; flipud(y)], rgb, ...
        'FaceAlpha', 0.1, 'EdgeColor', 'none'); % adjust transparency here
end
