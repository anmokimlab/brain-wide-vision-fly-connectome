clear all; close all; clc;
load('FB_Synapse_Layer.mat');
load('FB_Neuropil_Type.mat');
%% Nan 값처리
%%%%%% nan 인 것은 아예 visual 잆력이 하나도 없다는 뜻 0000 은 좌우에 하나는 있다는뜻
FB_Synapse_ME_R = fillmissing(FB_Synapse_ME_R,'constant',0);
FB_Synapse_ME_L = fillmissing(FB_Synapse_ME_L ,'constant',0);
FB_Upstream_Synapse_ME_R = fillmissing(FB_Upstream_Synapse_ME_R,'constant',0);
FB_Upstream_Synapse_ME_L = fillmissing(FB_Upstream_Synapse_ME_L,'constant',0);


FB_Synapse_LO_R = fillmissing(FB_Synapse_LO_R,'constant',0);
FB_Synapse_LO_L = fillmissing(FB_Synapse_LO_L ,'constant',0);
FB_Upstream_Synapse_LO_R = fillmissing(FB_Upstream_Synapse_LO_R,'constant',0);
FB_Upstream_Synapse_LO_L = fillmissing(FB_Upstream_Synapse_LO_L,'constant',0);

FB_Synapse_LOP_R = fillmissing(FB_Synapse_LOP_R,'constant',0);
FB_Synapse_LOP_L = fillmissing(FB_Synapse_LOP_L ,'constant',0);
FB_Upstream_Synapse_LOP_R = fillmissing(FB_Upstream_Synapse_LOP_R,'constant',0);
FB_Upstream_Synapse_LOP_L = fillmissing(FB_Upstream_Synapse_LOP_L,'constant',0);



%% 빈 갯수 반으로
% FB_Synapse_LOP_R_reduced = zeros(60, size(FB_Synapse_LOP_R,2));
% FB_Synapse_LO_R_reduced = zeros(60, size(FB_Synapse_LO_R,2));
% FB_Synapse_ME_R_reduced = zeros(60, size(FB_Synapse_ME_R,2));
% FB_Upstream_Synapse_LOP_R_reduced = zeros(60, size(FB_Upstream_Synapse_LOP_R,2));
% FB_Upstream_Synapse_LOP_L_reduced = zeros(60, size(FB_Upstream_Synapse_LOP_L,2));
% FB_Upstream_Synapse_LO_R_reduced = zeros(60, size(FB_Upstream_Synapse_LO_R,2));
% FB_Upstream_Synapse_LO_L_reduced = zeros(60, size(FB_Upstream_Synapse_LO_L,2));
% FB_Upstream_Synapse_ME_R_reduced = zeros(60, size(FB_Upstream_Synapse_ME_R,2));
% FB_Upstream_Synapse_ME_L_reduced = zeros(60, size(FB_Upstream_Synapse_ME_L,2));
% 
% for i = 1:60
%     FB_Synapse_LOP_R_reduced(i, :) = sum(FB_Synapse_LOP_R(2*i-1:2*i, :), 1);  % 행 2개씩 평균
%     FB_Synapse_LO_R_reduced(i, :) = sum(FB_Synapse_LO_R(2*i-1:2*i, :), 1);  % 행 2개씩 평균
%     FB_Synapse_ME_R_reduced(i, :) = sum(FB_Synapse_ME_R(2*i-1:2*i, :), 1);  % 행 2개씩 평균
%     FB_Upstream_Synapse_LOP_R_reduced(i, :) = sum(FB_Upstream_Synapse_LOP_R(2*i-1:2*i, :), 1);  % 행 2개씩 평균
%     FB_Upstream_Synapse_LOP_L_reduced(i, :) = sum(FB_Upstream_Synapse_LOP_L(2*i-1:2*i, :), 1);  % 행 2개씩 평균
%     FB_Upstream_Synapse_LO_R_reduced(i, :) = sum(FB_Upstream_Synapse_LO_R(2*i-1:2*i, :), 1);  % 행 2개씩 평균
%     FB_Upstream_Synapse_LO_L_reduced(i, :) = sum(FB_Upstream_Synapse_LO_L(2*i-1:2*i, :), 1);  % 행 2개씩 평균
%     FB_Upstream_Synapse_ME_R_reduced(i, :) = sum(FB_Upstream_Synapse_ME_R(2*i-1:2*i, :), 1);  % 행 2개씩 평균
%     FB_Upstream_Synapse_ME_L_reduced(i, :) = sum(FB_Upstream_Synapse_ME_L(2*i-1:2*i, :), 1);  % 행 2개씩 평균
% end
% FB_Synapse_LOP_R=FB_Synapse_LOP_R_reduced;
% FB_Synapse_LO_R=FB_Synapse_LO_R_reduced;
% FB_Synapse_ME_R=FB_Synapse_ME_R_reduced;
% FB_Upstream_Synapse_LOP_R=FB_Upstream_Synapse_LOP_R_reduced;
% FB_Upstream_Synapse_LOP_L=FB_Upstream_Synapse_LOP_L_reduced;
% FB_Upstream_Synapse_LO_R=FB_Upstream_Synapse_LO_R_reduced;
% FB_Upstream_Synapse_LO_L=FB_Upstream_Synapse_LO_L_reduced;
% FB_Upstream_Synapse_ME_R=FB_Upstream_Synapse_ME_R_reduced;
% FB_Upstream_Synapse_ME_L=FB_Upstream_Synapse_ME_L_reduced;
%%
%% Gaussian kernel 설정
sigma = 2;
x = -round(3*sigma):round(3*sigma);
gaussKernel = exp(-(x.^2)/(2*sigma^2));
gaussKernel = gaussKernel / sum(gaussKernel);  % 정규화
% sigma = 2, → kernel length = 13
%% Gaussian smoothing + overlap 계산

% LOP
for i = 1:size(FB_Synapse_LOP_R,2)
    X = FB_Synapse_LOP_R(:,i);
    X = conv(X, gaussKernel, 'same');
    X = X ./ sum(X);
    X=fillmissing(X,'constant',0);

    Y_L =  FB_Upstream_Synapse_LOP_L(:,i);
    Y_L = conv(Y_L, gaussKernel, 'same');
    Y_L = Y_L ./ sum(Y_L)*(sum(FB_Upstream_Synapse_LOP_L(:,i)))/(sum(FB_Upstream_Synapse_LOP_L(:,i))+sum(FB_Upstream_Synapse_LOP_R(:,i)));
    Y_L=fillmissing(Y_L,'constant',0);
    Y_R =  FB_Upstream_Synapse_LOP_R(:,i);
    Y_R = conv(Y_R, gaussKernel, 'same');
    Y_R = Y_R ./ sum(Y_R)*(sum(FB_Upstream_Synapse_LOP_R(:,i)))/(sum(FB_Upstream_Synapse_LOP_L(:,i))+sum(FB_Upstream_Synapse_LOP_R(:,i)));
    Y_R=fillmissing(Y_R,'constant',0);


    overlap_idx_LOP(i) = sum(min(X, (Y_R+Y_L)));
    overlap_idx_LOP_L(i) = sum(min(X, (Y_L)));
    overlap_idx_LOP_R(i) = sum(min(X, (Y_R)));


end

% LO
for i = 1:size(FB_Synapse_LO_R,2)
    X = FB_Synapse_LO_R(:,i);
    X = conv(X, gaussKernel, 'same');
    X = X ./ sum(X);
    X=fillmissing(X,'constant',0);

    Y_L =  FB_Upstream_Synapse_LO_L(:,i);
    Y_L = conv(Y_L, gaussKernel, 'same');
    Y_L = Y_L ./ sum(Y_L)*(sum(FB_Upstream_Synapse_LO_L(:,i)))/(sum(FB_Upstream_Synapse_LO_L(:,i))+sum(FB_Upstream_Synapse_LO_R(:,i)));
    Y_L=fillmissing(Y_L,'constant',0);

    Y_R =  FB_Upstream_Synapse_LO_R(:,i);
    Y_R = conv(Y_R, gaussKernel, 'same');
    Y_R = Y_R ./ sum(Y_R)*(sum(FB_Upstream_Synapse_LO_R(:,i)))/(sum(FB_Upstream_Synapse_LO_L(:,i))+sum(FB_Upstream_Synapse_LO_R(:,i)));
    Y_R=fillmissing(Y_R,'constant',0);
    % 
    % figure(1);
    % subplot(311);plot(X)
    % subplot(312);plot(Y_L)
    % subplot(313);plot(Y_R)

    overlap_idx_LO(i) = sum(min(X, (Y_R+Y_L)));
    overlap_idx_LO_L(i) = sum(min(X, (Y_L)));
    overlap_idx_LO_R(i) = sum(min(X, (Y_R)));


end

% ME
for i = 1:size(FB_Synapse_ME_R,2)
    X = FB_Synapse_ME_R(:,i);
    X = conv(X, gaussKernel, 'same');
    X = X ./ sum(X);
    X=fillmissing(X,'constant',0);

    Y_L =  FB_Upstream_Synapse_ME_L(:,i);
    Y_L = conv(Y_L, gaussKernel, 'same');
    Y_L = Y_L ./ sum(Y_L)*(sum(FB_Upstream_Synapse_ME_L(:,i)))/(sum(FB_Upstream_Synapse_ME_L(:,i))+sum(FB_Upstream_Synapse_ME_R(:,i)));
    Y_L=fillmissing(Y_L,'constant',0);

    Y_R =  FB_Upstream_Synapse_ME_R(:,i);
    Y_R = conv(Y_R, gaussKernel, 'same');
    Y_R = Y_R ./ sum(Y_R)*(sum(FB_Upstream_Synapse_ME_R(:,i)))/(sum(FB_Upstream_Synapse_ME_L(:,i))+sum(FB_Upstream_Synapse_ME_R(:,i)));
    Y_R=fillmissing(Y_R,'constant',0);


    overlap_idx_ME(i) = sum(min(X, (Y_R+Y_L)));
    overlap_idx_ME_L(i) = sum(min(X, (Y_L)));
    overlap_idx_ME_R(i) = sum(min(X, (Y_R)));


end

%% 평균 계산
mean_overlap_ME = mean(overlap_idx_ME, 'omitnan');
mean_overlap_ME_L = mean(overlap_idx_ME_L, 'omitnan');
mean_overlap_ME_R = mean(overlap_idx_ME_R, 'omitnan');

mean_overlap_LO = mean(overlap_idx_LO, 'omitnan');
mean_overlap_LO_L = mean(overlap_idx_LO_L, 'omitnan');
mean_overlap_LO_R = mean(overlap_idx_LO_R, 'omitnan');

mean_overlap_LOP = mean(overlap_idx_LOP, 'omitnan');
mean_overlap_LOP_L = mean(overlap_idx_LOP_L, 'omitnan');
mean_overlap_LOP_R = mean(overlap_idx_LOP_R, 'omitnan');


%% 시각화
figure(1); set(gcf,'Color','w'); hold on;
y = [overlap_idx_ME_L' overlap_idx_ME_R'];

b = bar(FB_Me.type, y);
b(1).FaceColor = [0.4660 0.6740 0.1880];  % 첫째 열 색 (RGB 0~1)
b(2).FaceColor = [0 0.4470 0.7410];  % 둘째 열 색

set(gca, 'TickDir', 'out', 'Box', 'off');
ylim([0 1]);
ylabel('Overlap coefficient (Gaussian smoothed)');
title('Overlap of input and output depth distributions (with Gaussian smoothing)');
print(gcf, '-depsc2', '-vector', 'FB_Layer_overlap_EachNeuron_ME_bar.eps');

%%
figure(2); set(gcf,'Color','w'); hold on;
y = [overlap_idx_LO_L' overlap_idx_LO_R'];

b = bar(FB_Lo.type, y);
b(1).FaceColor = [0.4660 0.6740 0.1880];  % 첫째 열 색 (RGB 0~1)
b(2).FaceColor = [0 0.4470 0.7410];  % 둘째 열 색

set(gca, 'TickDir', 'out', 'Box', 'off');
ylim([0 1]);
ylabel('Overlap coefficient (Gaussian smoothed)');
title('Overlap of input and output depth distributions (with Gaussian smoothing)');
print(gcf, '-depsc2', '-vector', 'FB_Layer_overlap_EachNeuron_LO_bar.eps');

%%

figure(3); set(gcf,'Color','w'); hold on;
y = [overlap_idx_LOP_L' overlap_idx_LOP_R'];

b = bar(FB_Lop.type, y);
b(1).FaceColor = [0.4660 0.6740 0.1880];  % 첫째 열 색 (RGB 0~1)
b(2).FaceColor = [0 0.4470 0.7410];  % 둘째 열 색

set(gca, 'TickDir', 'out', 'Box', 'off');
ylim([0 1]);
ylabel('Overlap coefficient (Gaussian smoothed)');
title('Overlap of input and output depth distributions (with Gaussian smoothing)');
print(gcf, '-depsc2', '-vector', 'FB_Layer_overlap_EachNeuron_LOP_bar.eps');

%% 전체 뉴런 기준 평균 및 표준편차 계산

% 모든 overlap 값을 하나의 배열로 결합
all_overlap = [overlap_idx_ME, overlap_idx_LO, overlap_idx_LOP];
all_overlap_L = [overlap_idx_ME_L, overlap_idx_LO_L, overlap_idx_LOP_L];
all_overlap_R = [overlap_idx_ME_R, overlap_idx_LO_R, overlap_idx_LOP_R];

% 전체 평균
mean_overlap_all = mean(all_overlap, 'omitnan');
mean_overlap_all_L = mean(all_overlap_L, 'omitnan');
mean_overlap_all_R = mean(all_overlap_R, 'omitnan');

% 전체 표준편차
std_overlap_all = std(all_overlap, 0, 'omitnan');
std_overlap_all_L = std(all_overlap_L, 0, 'omitnan');
std_overlap_all_R = std(all_overlap_R, 0, 'omitnan');

%% 결과 출력
fprintf('전체 뉴런 기준 overlap 평균 및 표준편차:\n');
fprintf('- 전체: mean = %.4f, std = %.4f\n', mean_overlap_all, std_overlap_all);
fprintf('- Left : mean = %.4f, std = %.4f\n', mean_overlap_all_L, std_overlap_all_L);
fprintf('- Right: mean = %.4f, std = %.4f\n', mean_overlap_all_R, std_overlap_all_R);

%% 영역별 평균 및 표준편차 정리
mean_vals = [
    mean(overlap_idx_ME_L, 'omitnan'), mean(overlap_idx_ME_R, 'omitnan');
    mean(overlap_idx_LO_L, 'omitnan'), mean(overlap_idx_LO_R, 'omitnan');
    mean(overlap_idx_LOP_L, 'omitnan'), mean(overlap_idx_LOP_R, 'omitnan')
];

std_vals = [
    std(overlap_idx_ME_L, 0, 'omitnan'), std(overlap_idx_ME_R, 0, 'omitnan');
    std(overlap_idx_LO_L, 0, 'omitnan'), std(overlap_idx_LO_R, 0, 'omitnan');
    std(overlap_idx_LOP_L, 0, 'omitnan'), std(overlap_idx_LOP_R, 0, 'omitnan')
];

%% 바 그래프 그리기
figure(4); clf; set(gcf,'Color','w');

bar_handle = bar(mean_vals, 'grouped');
hold on;

% 색 설정
bar_handle(1).FaceColor = [0.4660 0.6740 0.1880]; % Left
bar_handle(2).FaceColor = [0 0.4470 0.7410];      % Right

% 에러바 추가
ngroups = size(mean_vals, 1);
nbars = size(mean_vals, 2);
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, mean_vals(:,i), std_vals(:,i), 'k', 'linestyle', 'none', 'LineWidth', 1);
end

% 라벨 및 제목
set(gca, 'XTickLabel', {'ME', 'LO', 'LOP'});
ylabel('Mean overlap coefficient ± std');
ylim([0 0.8]);
legend({'Left', 'Right'}, 'Location', 'northeast');
title('Overlap between input and output depths (mean ± std)');
set(gca, 'TickDir', 'out', 'Box', 'off');

print(gcf, '-depsc2', '-vector', 'FB_Layer_overlap_grouped_bar.eps');

