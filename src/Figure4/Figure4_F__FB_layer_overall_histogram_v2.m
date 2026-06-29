clear all; close all; clc;

%% ===== Load =====
load('FB_Synapse_Layer.mat');
load('FB_Neuropil_Type.mat');

%% ===== Depth tick labels =====
edges = -7e-5:1e-6:5e-5;

%% ===== NaN 처리 (너의 의도 유지: NaN=시각입력 없음 -> 0으로 간주) =====
FB_Synapse_ME_R           = fillmissing(FB_Synapse_ME_R,'constant',0);
FB_Synapse_ME_L           = fillmissing(FB_Synapse_ME_L,'constant',0);
FB_Upstream_Synapse_ME_R  = fillmissing(FB_Upstream_Synapse_ME_R,'constant',0);
FB_Upstream_Synapse_ME_L  = fillmissing(FB_Upstream_Synapse_ME_L,'constant',0);

FB_Synapse_LO_R           = fillmissing(FB_Synapse_LO_R,'constant',0);
FB_Synapse_LO_L           = fillmissing(FB_Synapse_LO_L,'constant',0);
FB_Upstream_Synapse_LO_R  = fillmissing(FB_Upstream_Synapse_LO_R,'constant',0);
FB_Upstream_Synapse_LO_L  = fillmissing(FB_Upstream_Synapse_LO_L,'constant',0);

FB_Synapse_LOP_R          = fillmissing(FB_Synapse_LOP_R,'constant',0);
FB_Synapse_LOP_L          = fillmissing(FB_Synapse_LOP_L,'constant',0);
FB_Upstream_Synapse_LOP_R = fillmissing(FB_Upstream_Synapse_LOP_R,'constant',0);
FB_Upstream_Synapse_LOP_L = fillmissing(FB_Upstream_Synapse_LOP_L,'constant',0);

%% ===== Colors (keep your palette) =====
col_orange = [0.8500 0.3250 0.0980];   % Synapse_R
col_blue   = [0      0.4470 0.7410];   % Upstream_R
col_green  = [0.4660 0.6740 0.1880];   % Upstream_L

%% ======================= Figure 1: ME =======================
figure(1); clf; set(gcf,'Color','w'); hold on;

y_ME = (1:size(FB_Synapse_ME_R,1))';

% IQR bands + median lines
shaded_iqr(FB_Synapse_ME_R,          y_ME, col_orange); % 주황
shaded_iqr(FB_Upstream_Synapse_ME_R, y_ME, col_blue);   % 파랑
shaded_iqr(FB_Upstream_Synapse_ME_L, y_ME, col_green);  % 초록

xlabel('Normalized value');
ylabel('Innvervation depth');
title('ME innervation depth (median & IQR)');
set(gca,'TickDir','out','Box','off');
ylim([0.5 90.5]); xlim([0 0.8]);
set(gca,'YTick',0.5:10:90+0.5,'YTickLabel',edges(1:10:end));

% print(gcf,'-depsc2','-vector','Me_Layer_IQR.eps');

%% ======================= Figure 2: LO =======================
figure(2); clf; set(gcf,'Color','w'); hold on;

y_LO = (1:size(FB_Synapse_LO_R,1))';

shaded_iqr(FB_Synapse_LO_R,          y_LO, col_orange);
shaded_iqr(FB_Upstream_Synapse_LO_R, y_LO, col_blue);
shaded_iqr(FB_Upstream_Synapse_LO_L, y_LO, col_green);

xlabel('Normalized value');
ylabel('Innvervation depth');
title('LO innervation depth (median & IQR)');
set(gca,'TickDir','out','Box','off');
set(gca,'YTick',0.5:10:100+0.5,'YTickLabel',edges(1:10:end));
ylim([25.5 90.5]); xlim([0 0.8]);

% print(gcf,'-depsc2','-vector','Lo_Layer_IQR.eps');

%% ======================= Figure 3: LOP =======================
figure(3); clf; set(gcf,'Color','w'); hold on;

y_LOP = (1:size(FB_Synapse_LOP_R,1))';

shaded_iqr(FB_Synapse_LOP_R,          y_LOP, col_orange);
shaded_iqr(FB_Upstream_Synapse_LOP_R, y_LOP, col_blue);
shaded_iqr(FB_Upstream_Synapse_LOP_L, y_LOP, col_green);

set(gca,'YDir','reverse');  % LOP은 위->아래
xlabel('Normalized value');
ylabel('Innvervation depth');
title('LOP innervation depth (median & IQR)');
set(gca,'TickDir','out','Box','off');
set(gca,'YTick',0.5:10:100+0.5,'YTickLabel',edges(1:10:end));
ylim([60.5 100.5]); xlim([0 0.8]);

% print(gcf,'-depsc2','-vector','LoP_Layer_IQR.eps');

%% ===== Local function: draw median + IQR (25–75%) shaded band =====
function [hBand,hMed] = shaded_iqr(X, yvec, rgb)
    % X    : (depth x samples) 행=깊이, 열=개체
    % yvec : (depth x 1)
    % rgb  : 1x3 color
    % NOTE: 위 본문에서 이미 NaN->0으로 대체했으니, 여기선 그대로 사용.
    %      (만약 NaN을 그대로 유지했다면, prctile이 NaN에 민감하므로 행별로 rmmissing이 필요)
    q25 = prctile(X, 25, 2);
    q75 = prctile(X, 75, 2);
    med = median(X, 2);  % 중앙값 선

    ok = ~(isnan(q25) | isnan(q75) | isnan(yvec));
    q25 = q25(ok); q75 = q75(ok); y = yvec(ok);
    med = med(ok);

    % IQR shaded band
    hBand = fill([q25; flipud(q75)], [y; flipud(y)], rgb, ...
                 'FaceAlpha', 0.18, 'EdgeColor', 'none'); hold on
    % median line
    hMed  = plot(med, y, 'Color', rgb, 'LineWidth', 2);
end
