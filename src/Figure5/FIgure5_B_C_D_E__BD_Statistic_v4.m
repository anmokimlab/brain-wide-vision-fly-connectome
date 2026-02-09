%% Clean
clear all; close all; clc

%% Paths & Loads
opt = detectImportOptions('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\consolidated_cell_types.csv');
opt = setvartype(opt,'root_id','int64');
FAFBConsolidated_type = readtable('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\consolidated_cell_types.csv',opt);

load('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\Figure1_Overall\Right_Neurons_Thr0.mat')
load('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\Figure1_Overall\FAFB_NPI_Thr0.mat')

opt = detectImportOptions('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\connections_no_threshold.csv');
opt = setvartype(opt,'pre_root_id','int64');
opt = setvartype(opt,'post_root_id','int64');
FAFBConnections = readtable('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\connections_no_threshold.csv',opt);

FAFBNeuropils = unique(FAFBConnections.neuropil);
FAFBNeuropil_OpticLobeRight = {'LA_R';'AME_R';'ME_R';'LO_R';'LOP_R'};
FAFBNeuropil_OpticLobeLeft  = {'LA_L';'AME_L';'ME_L';'LO_L';'LOP_L'};
FAFBNeuropil_Central = FAFBNeuropils(~ismember(FAFBNeuropils,FAFBNeuropil_OpticLobeRight));
FAFBNeuropil_Central = FAFBNeuropil_Central(~ismember(FAFBNeuropil_Central,FAFBNeuropil_OpticLobeLeft));
FAFBNeuropil_Central = FAFBNeuropil_Central(~ismember(FAFBNeuropil_Central,'UNASGD'));

%% 여러 타입 입력
% 예: WantToSee = {'LC9','LC10b'};
WantToSee = {'LC9','LT43','LT52','aMe1','cLM01'};

% 각 타입별 right side root_id 수집
Want = struct('type',[],'root_ids',[],'tbl',[]);
all_want_root_ids = [];

% 안전한 소스 테이블 리스트 구성 (존재하는 것만 사용)
sourceTbls = {};
if exist('FAFBNPIs','var') && istable(FAFBNPIs), sourceTbls{end+1} = FAFBNPIs; end
if exist('RightFF_NPIs','var') && istable(RightFF_NPIs), sourceTbls{end+1} = RightFF_NPIs; end
if exist('RightFB_NPIs','var') && istable(RightFB_NPIs), sourceTbls{end+1} = RightFB_NPIs; end

for k = 1:numel(WantToSee)
    Want(k).type = WantToSee{k};
    rid = [];

    % 1차: NPI 기반 테이블에서 타입 매칭
    for si = 1:numel(sourceTbls)
        T = sourceTbls{si};
        hasInField = ismember('In_Synapse_Optic_R', T.Properties.VariableNames);
        if hasInField
            rid = [rid; T.root_id(strcmp(T.type, Want(k).type) & T.In_Synapse_Optic_R>0)];
        else
            rid = [rid; T.root_id(strcmp(T.type, Want(k).type))];
        end
    end

    % 2차: 그래도 없으면 consolidated_cell_types에서 primary_type으로 매칭
    if isempty(rid)
        if ismember('primary_type', FAFBConsolidated_type.Properties.VariableNames)
            rid = [rid; FAFBConsolidated_type.root_id(strcmp(FAFBConsolidated_type.primary_type, Want(k).type))];
        end
    end

    rid = unique(rid(:));
    Want(k).root_ids = rid;
    all_want_root_ids = [all_want_root_ids; rid]; %#ok<AGROW>
end
all_want_root_ids = unique(all_want_root_ids);

%% FF/FB 비교군에서 Wanted root_id 제거
FF_neuropils = RightFF_NPIs;
if ~isempty(all_want_root_ids)
    FF_neuropils(ismember(FF_neuropils.root_id, all_want_root_ids),:) = [];
end

FB_neuropils = RightFB_NPIs;
if ~isempty(all_want_root_ids)
    FB_neuropils(ismember(FB_neuropils.root_id, all_want_root_ids),:) = [];
end

%% 타입별 FF 요약 (in/out 비율)
[FF_Type,~,ic] = unique(FF_neuropils.type);
Type_FF_In_Out_ratio = table(FF_Type,'VariableNames',{'type'});
for i=1:size(Type_FF_In_Out_ratio,1)
    idx = (ic==i);
    Type_FF_In_Out_ratio.root_id{i} = FF_neuropils.root_id(idx);
    in_temp  = [FF_neuropils.In_Synapse_Optic_R(idx) FF_neuropils.In_Synapse_Central(idx)];
    out_temp = [FF_neuropils.Out_Synapse_Optic_R(idx) FF_neuropils.Out_Synapse_Central(idx)];
    Type_FF_In_Out_ratio.in_ratio{i}  = mean(bsxfun(@rdivide,in_temp,  sum(in_temp,2)),1,'omitnan');
    Type_FF_In_Out_ratio.out_ratio{i} = mean(bsxfun(@rdivide,out_temp, sum(out_temp,2)),1,'omitnan');
end

%% 타입별 FB 요약 (in/out 비율)
[FB_Type,~,ic] = unique(FB_neuropils.type);
Type_FB_In_Out_ratio = table(FB_Type,'VariableNames',{'type'});
for i=1:size(Type_FB_In_Out_ratio,1)
    idx = (ic==i);
    Type_FB_In_Out_ratio.root_id{i} = FB_neuropils.root_id(idx);
    in_temp  = [FB_neuropils.In_Synapse_Optic_R(idx) FB_neuropils.In_Synapse_Central(idx)];
    out_temp = [FB_neuropils.Out_Synapse_Optic_R(idx) FB_neuropils.Out_Synapse_Central(idx)];
    Type_FB_In_Out_ratio.in_ratio{i}  = mean(bsxfun(@rdivide,in_temp,  sum(in_temp,2)),1,'omitnan');
    Type_FB_In_Out_ratio.out_ratio{i} = mean(bsxfun(@rdivide,out_temp, sum(out_temp,2)),1,'omitnan');
end

%%
data1 = cell2mat(Type_FF_In_Out_ratio.in_ratio);
data2 = cell2mat(Type_FF_In_Out_ratio.out_ratio);
G = cell(1, 2 + numel(WantToSee)); % 미리 크기 지정 (작은 최적화)
G{1} = [data1(:,1) data2(:,1)];

data1 = cell2mat(Type_FB_In_Out_ratio.in_ratio);
data2 = cell2mat(Type_FB_In_Out_ratio.out_ratio);
G{2} = [data1(:,1) data2(:,1)];

%% 각 Want 타입별: root별 합계 테이블 만들고, In/Out 누적막대 그리기
fig_counter = 3;  % 5부터 이어서 그림
for k = 1:numel(Want)
    rid = Want(k).root_ids;
    Want_tbl = table(rid,'VariableNames',{'root_id'});
    Want_tbl.In_Synapse_Optic_R  = zeros(height(Want_tbl),1);
    Want_tbl.Out_Synapse_Optic_R = zeros(height(Want_tbl),1);
    Want_tbl.In_Synapse_Central  = zeros(height(Want_tbl),1);
    Want_tbl.Out_Synapse_Central = zeros(height(Want_tbl),1);
    Want_tbl.In_Synapse_Optic_L  = zeros(height(Want_tbl),1);
    Want_tbl.Out_Synapse_Optic_L = zeros(height(Want_tbl),1);

    for i=1:height(Want_tbl)
        current_root_id = Want_tbl.root_id(i);
        InConnections  = FAFBConnections(FAFBConnections.post_root_id==current_root_id, :);
        OutConnections = FAFBConnections(FAFBConnections.pre_root_id==current_root_id,  :);

        Want_tbl.In_Synapse_Optic_R(i)  = sum(InConnections.syn_count( ismember(InConnections.neuropil,FAFBNeuropil_OpticLobeRight) ));
        Want_tbl.Out_Synapse_Optic_R(i) = sum(OutConnections.syn_count(ismember(OutConnections.neuropil,FAFBNeuropil_OpticLobeRight)));
        Want_tbl.In_Synapse_Central(i)  = sum(InConnections.syn_count( ismember(InConnections.neuropil,FAFBNeuropil_Central) ));
        Want_tbl.Out_Synapse_Central(i) = sum(OutConnections.syn_count(ismember(OutConnections.neuropil,FAFBNeuropil_Central)));
        Want_tbl.In_Synapse_Optic_L(i)  = sum(InConnections.syn_count( ismember(InConnections.neuropil,FAFBNeuropil_OpticLobeLeft) ));
        Want_tbl.Out_Synapse_Optic_L(i) = sum(OutConnections.syn_count(ismember(OutConnections.neuropil,FAFBNeuropil_OpticLobeLeft)));
    end

    Want(k).tbl = Want_tbl;

    data = [Want_tbl.In_Synapse_Optic_R./(Want_tbl.In_Synapse_Optic_R+Want_tbl.In_Synapse_Central) ...
        Want_tbl.Out_Synapse_Optic_R./(Want_tbl.Out_Synapse_Optic_R+Want_tbl.Out_Synapse_Central)];
    G{fig_counter} = data;
    fig_counter=fig_counter+1;
end

%%
%% =====================================================================
%  GMM을 이용한 소속 확률 계산 (기존 코드의 G와 labels 변수 사용)
%  =====================================================================

%% 1. FF와 FB 데이터 준비
labels = {'FF','FB'};  % 그룹 라벨
labels=[labels WantToSee];
FF_data = G{1};
FB_data = G{2};
WantToSee_labels = labels(3:end);

% NaN 데이터 제거
FF_data(any(isnan(FF_data), 2), :) = [];
FB_data(any(isnan(FB_data), 2), :) = [];

% 안전 장치: 데이터가 충분한지 확인
if size(FF_data,1) < 2 || size(FB_data,1) < 2
    error('FF/FB 데이터가 너무 적습니다. FF=%d, FB=%d', size(FF_data,1), size(FB_data,1));
end

%% 2. 각 그룹에 대한 최적의 GMM 모델 학습 (BIC 기준)
rng(1); % 재현 가능성

% k_range를 표본 수에 맞춰 안전하게 제한
maxK_ff = max(1, min(10, size(FF_data,1)-1));
maxK_fb = max(1, min(10, size(FB_data,1)-1));
k_range_ff = 1:maxK_ff;
k_range_fb = 1:maxK_fb;

% --- FF 모델 학습 ---
bic_scores_ff = inf(size(k_range_ff));
gm_models_ff = cell(size(k_range_ff));
for i = 1:numel(k_range_ff)
    try
        gm_models_ff{i} = fitgmdist(FF_data, k_range_ff(i), ...
            'Options', statset('MaxIter', 2000), ...
            'SharedCovariance', false, ...
            'CovarianceType', 'full', ...
            'RegularizationValue', 1e-6, ...
            'Replicates', 10);
        bic_scores_ff(i) = gm_models_ff{i}.BIC;
    catch
        % 실패 시 해당 K는 건너뜀 (BIC는 inf 유지)
    end
end
[~, best_k_idx] = min(bic_scores_ff);
best_gm_ff = gm_models_ff{best_k_idx};
if isempty(best_gm_ff)
    error('FF GMM 적합 실패: 다른 파라미터를 시도하세요.');
end
fprintf('FF 그룹의 최적 클러스터 개수 (BIC 기준): %d\n', best_gm_ff.NumComponents);

% --- FB 모델 학습 ---
bic_scores_fb = inf(size(k_range_fb));
gm_models_fb = cell(size(k_range_fb));
for i = 1:numel(k_range_fb)
    try
        gm_models_fb{i} = fitgmdist(FB_data, k_range_fb(i), ...
            'Options', statset('MaxIter', 2000), ...
            'SharedCovariance', false, ...
            'CovarianceType', 'full', ...
            'RegularizationValue', 1e-6, ...
            'Replicates', 10);
        bic_scores_fb(i) = gm_models_fb{i}.BIC;
    catch
    end
end
[~, best_k_idx] = min(bic_scores_fb);
best_gm_fb = gm_models_fb{best_k_idx};
if isempty(best_gm_fb)
    error('FB GMM 적합 실패: 다른 파라미터를 시도하세요.');
end
fprintf('FB 그룹의 최적 클러스터 개수 (BIC 기준): %d\n', best_gm_fb.NumComponents);

%% 3. WantToSee 그룹들의 소속 확률 계산 및 결과 출력
% for i = 1:numel(WantToSee_labels)
%     group_label = WantToSee_labels{i};
%     group_data = G{i+2};
%     group_data(any(isnan(group_data), 2), :) = [];
% 
%     if isempty(group_data)
%         fprintf('\n--- [%s] 그룹은 분석할 데이터가 없습니다 ---\n', group_label);
%         continue;
%     end
% 
%     % --- 각 개별 뉴런의 모델 로그가능도 계산 (posterior 아님!) ---
%     % posterior의 합은 1이므로 모델 비교가 되지 않음 -> logpdf 사용
%     logp_ff = log(pdf(best_gm_ff, group_data));
%     logp_fb = log(pdf(best_gm_fb, group_data));
% 
%     % FF에 더 가깝다고 판단된 뉴런의 수
%     is_ff = logp_ff > logp_fb;
%     ff_closer_count = sum(is_ff);
%     fb_closer_count = size(group_data, 1) - ff_closer_count;
% 
%     % --- 그룹 평균 지점의 로그가능도 비교 ---
%     group_mean = mean(group_data, 1);
%     loglik_ff_mean = log(pdf(best_gm_ff, group_mean));
%     loglik_fb_mean = log(pdf(best_gm_fb, group_mean));
% 
%     % 결과 출력
%     fprintf('\n--- [%s] 그룹 분석 결과 ---\n', group_label);
%     fprintf('총 뉴런 수: %d\n', size(group_data, 1));
%     fprintf('개별 뉴런 분석:\n');
%     fprintf('  - FF 모델 로그가능도 우세: %d개 (%.1f%%)\n', ff_closer_count, (ff_closer_count/size(group_data, 1))*100);
%     fprintf('  - FB 모델 로그가능도 우세: %d개 (%.1f%%)\n', fb_closer_count, (fb_closer_count/size(group_data, 1))*100);
% 
%     fprintf('그룹 평균 지점 분석 (Log-Likelihood):\n');
%     fprintf('  - FF 모델 적합도: %.4f\n', loglik_ff_mean);
%     fprintf('  - FB 모델 적합도: %.4f\n', loglik_fb_mean);
% 
%     if loglik_ff_mean > loglik_fb_mean
%         fprintf('  >> 결론: 그룹 평균은 FF 모델에 더 가깝습니다.\n');
%     else
%         fprintf('  >> 결론: 그룹 평균은 FB 모델에 더 가깝습니다.\n');
%     end
% end
% --- [교체] 각 그룹의 사후 로그오즈(LLR*)/사후확률 계산 (사전 포함) ---
% for i = 1:numel(WantToSee_labels)
%     group_label = WantToSee_labels{i};
%     Xi = G{i+2};
%     Xi(any(isnan(Xi),2),:) = [];
%     if isempty(Xi)
%         fprintf('\n--- [%s] 그룹은 분석할 데이터가 없습니다 ---\n', group_label);
%         continue;
%     end
% 
%     % 사전 포함 점수
%     score_ff = log(pdf(best_gm_ff, Xi)) + log(max(pi_FF, realmin));
%     score_fb = log(pdf(best_gm_fb, Xi)) + log(max(pi_FB, realmin));
%     LLR_star_i = score_ff - score_fb;
%     Pff_i = 1 ./ (1 + exp(-LLR_star_i));
% 
%     % 개체 수준 요약
%     ff_closer_count = sum(LLR_star_i > 0);      % posterior > 0.5와 동일
%     fb_closer_count = numel(LLR_star_i) - ff_closer_count;
% 
%     % 그룹 수준 집계(권장): 평균 LLR*, 평균 posterior, 총 evidence(로그합)
%     meanLLR = mean(LLR_star_i);
%     meanP   = mean(Pff_i);
%     logEvFF = sum(score_ff);   % 독립 가정 하 집단 evidence
%     logEvFB = sum(score_fb);
% 
%     fprintf('\n--- [%s] 그룹 분석 결과 (사전 포함) ---\n', group_label);
%     fprintf('총 뉴런 수: %d\n', numel(Pff_i));
%     fprintf('개별 뉴런:\n');
%     fprintf('  - FF 우세(LLR*>0): %d개 (%.1f%%)\n', ff_closer_count, 100*ff_closer_count/numel(Pff_i));
%     fprintf('  - FB 우세(LLR*<=0): %d개 (%.1f%%)\n', fb_closer_count, 100*fb_closer_count/numel(Pff_i));
%     fprintf('그룹 요약:\n');
%     fprintf('  - 평균 LLR*: %+0.4f  (양수=FF 우세)\n', meanLLR);
%     fprintf('  - 평균 P(FF|x): %0.3f\n', meanP);
%     fprintf('  - 로그 evidence 합: FF=%0.4f, FB=%0.4f  => %s 우세\n', ...
%             logEvFF, logEvFB, ternary(logEvFF>logEvFB,'FF','FB'));
% end
% 
% function out = ternary(cond, a, b)
%     if cond, out = a; else, out = b; end
% end


%% =====================================================================
%  GMM 결과 시각화
%  =====================================================================

%% =====================================================================
%  2D 결정영역(FF/FB 근접) 시각화: P(FF|x) 히트맵 + 경계선 + 샘플/타원
%  (best_gm_ff, best_gm_fb, FF_data, FB_data, WantToSee_labels, G 필요)
%  =====================================================================

% 1) 사전확률(표본 크기 기반) 계산
N_FF = size(FF_data,1);
N_FB = size(FB_data,1);
pi_FF = N_FF / (N_FF + N_FB);
pi_FB = 1 - pi_FF;
log_prior_ratio = log(pi_FF) - log(pi_FB);   % log(pi_FF/pi_FB)

% 2) 격자 생성
n = 350;                             % 해상도(필요시 200~500 사이 조절)
xgrid = linspace(0,1,n);
[X,Y] = meshgrid(xgrid, xgrid);
P = [X(:) Y(:)];

% 3) 두 GMM의 로그가능도 및 LLR* 계산
logp_ff = log(pdf(best_gm_ff, P));
logp_fb = log(pdf(best_gm_fb, P));
LLR_star = logp_ff - logp_fb + log_prior_ratio;   % 최종 로그우도비
Pff = 1./(1+exp(-LLR_star));                      % P(FF|x) = sigma(LLR*)

Pff_map = reshape(Pff, n, n);     % (행:y, 열:x) 매핑
LLR_map = reshape(LLR_star, n, n);

% 4) 배경 히트맵 (확률) + 경계 (LLR*=0 <=> P=0.5)
figure('Color','w','Position',[100,100,900,780]); 
imagesc(xgrid, xgrid, Pff_map); 
set(gca,'YDir','normal'); 
axis([0 1 0 1]); axis square; grid on; box on;
c = colorbar; c.Label.String = 'P(FF | x)';
% 사용자 정의 colormap (gradient)
FF_color = [0.00, 0.45, 0.74];   % #0072BD (파랑)
FB_color = [0.47, 0.67, 0.19];   % #77AC30 (녹색)

ncol = 256;  % colormap 해상도
custom_cmap = [linspace(FB_color(1), FF_color(1), ncol)', ...
               linspace(FB_color(2), FF_color(2), ncol)', ...
               linspace(FB_color(3), FF_color(3), ncol)'];

colormap(custom_cmap);
% colormap("winter");

caxis([0 1]);   % 확률값 0→FB(녹색), 1→FF(파랑)
colorbar;

hold on;
contour(xgrid, xgrid, Pff_map, [0.5 0.5], 'k-', 'LineWidth', 1); % 결정경계


p_ff = scatter(FF_data(:,1), FF_data(:,2), 18, FF_color*0.5, 'filled', ...
    'MarkerFaceAlpha', 0.35, 'MarkerEdgeAlpha', 0.35);
p_fb = scatter(FB_data(:,1), FB_data(:,2), 18, FB_color*0.5, 'filled', ...
    'MarkerFaceAlpha', 0.35, 'MarkerEdgeAlpha', 0.35);
LC9_data= G{3};
LT43_data= G{4};
LT52_data= G{5};
aMe01_data= G{6};
cLM01_data = G{7};

p_cLM01 = scatter(cLM01_data(:,1), cLM01_data(:,2), 30, 'filled', ...
    'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1);
p_LC9 = scatter(LC9_data(:,1), LC9_data(:,2), 30, 'filled', ...
    'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1);
p_LT43 = scatter(LT43_data(:,1), LT43_data(:,2), 30, 'filled', ...
    'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1);
p_LT52 = scatter(LT52_data(:,1), LT52_data(:,2), 30, 'filled', ...
    'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1);
p_aMe01 = scatter(aMe01_data(:,1), aMe01_data(:,2), 30, 'filled', ...
    'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1);
legend({'Middle','FF','FB','cLM01','LC9','LT43','LT52','aMe01'},'Location','bestoutside')


%% ======================================================
%  WantToSee 전용 지터 플롯 (FF/FB 제외)
%  각 뉴런별 P(FF|x) + 그룹 평균 ± SD
% ======================================================

% 라벨: WantToSee_labels 있으면 사용, 없으면 WantToSee 사용
if exist('WantToSee_labels','var') && ~isempty(WantToSee_labels)
    lbls = string(WantToSee_labels(:));
else
    lbls = string(WantToSee(:));
end
lbls_cell = cellstr(lbls);

% 색 팔레트 (FF/FB 제외하고 3번째 색부터 사용)
C0 = [
    183  42  49;   % #B72A31  (LC9)
    152 122  63;   % #987A3F  (LT43)
    108  41 119;   % #6C2977  (LT52)
    43  170 184;    % #2BAAB8  (aMe1)
    228 172  41;   % #E4AC29  (cLM01)

] / 255;

% 라벨 수에 맞게 확장/절단
if size(C0,1) < numel(lbls)
    reps = ceil(numel(lbls)/size(C0,1));
    C = repmat(C0, reps, 1);
    C = C(1:numel(lbls),:);
else
    C = C0(1:numel(lbls),:);
end

figure('Color','w','Position',[100,100,1100,520]); 
ax = gca; hold(ax,'on'); box(ax,'off'); grid(ax,'on');
set(gca,'TickDir','out')
% 기준선
yline(0.5,'k--','LineWidth',1);

rng(42);                 % 지터 재현성
jitter_sigma = 0.06;     % x 지터 표준편차

for i = 1:numel(lbls)
    % --- 해당 그룹 데이터 (WantToSee 순서대로 G{3},G{4},...) ---
    % WantToSee_labels와 lbls가 동일 순서라면 바로 i+2 인덱스 사용
    Xi = G{i+2};                   
    Xi(any(isnan(Xi),2),:) = [];
    if isempty(Xi), continue; end

    % --- posterior P(FF|x) 계산 ---
    s_ff  = log(pdf(best_gm_ff, Xi)) + log(max(pi_FF, realmin));
    s_fb  = log(pdf(best_gm_fb, Xi)) + log(max(pi_FB, realmin));
    LLR   = s_ff - s_fb;
    Pff_i = 1 ./ (1 + exp(-LLR));   % 각 뉴런의 P(FF|x)

    % --- 원자료 지터 산점 ---
    xj = i + jitter_sigma*randn(numel(Pff_i),1);
    scatter(xj, Pff_i, 22, C(i,:), 'filled', ...
        'MarkerFaceAlpha', 0.38, 'MarkerEdgeColor','none');

    % --- 평균 + SD 에러바 ---
    m  = mean(Pff_i);
    if numel(Pff_i) >= 2
        sd = std(Pff_i);
        errorbar(i, m, sd, 'Color', C(i,:), 'CapSize', 6, 'LineWidth', 1.3);
    end
    scatter(i, m, 40, C(i,:), 'filled', 'MarkerEdgeColor','k');  % 평균점

    % 샘플 수 표기
    text(i, 1.02, sprintf('n=%d', numel(Pff_i)), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',9);
end

% 축/레이아웃
xlim([0.5, numel(lbls)+0.5]); 
ylim([0, 1.05]);
xticks(1:numel(lbls));
xticklabels(lbls_cell); xtickangle(15);
ylabel('Posterior  P(FF | x)');
title('Per-type posterior (WantToSee only): jitter + mean ± SD');
hold off;

%%
FF_postPI=G{1}(:,1);
FB_postPI=G{2}(:,1);
LC9_postPI=G{3}(:,1);
LT43_postPI=G{4}(:,1);
LT52_postPI=G{5}(:,1);
aMe1_postPI=G{6}(:,1);
cLM01_postPI=G{7}(:,1);

FF_prePI=G{1}(:,2);
FB_prePI=G{2}(:,2);
LC9_prePI=G{3}(:,2);
LT43_prePI=G{4}(:,2);
LT52_prePI=G{5}(:,2);
aMe1_prePI=G{6}(:,2);
cLM01_prePI=G{7}(:,2);
%%
C = [
    20  97 154;   % #14619A
    111 151  51;  % #6F9733
    183  42  49;   % #B72A31  (LC9)
    152 122  63;   % #987A3F  (LT43)
    108  41 119;   % #6C2977  (LT52)
    43  170 184    % #2BAAB8  (aMe1)
    228 172  41;   % #E4AC29  (cLM01)
] / 255;
figure('Color','w'); clf; ax = gca; hold(ax,'on');
% 각 그룹: 원자료 지터 + 평균점 + SD 에러바
lbls=["FF";"FB"; "LC9";"LT43";"LT52";"aMe1";"cLM01"];
lbls_str = string(lbls);
for i = 1:numel(lbls)
    switch lbls_str(i)
        case "FF",   yi = FF_postPI(:);
        case "FB",   yi = FB_postPI(:);
        case "cLM01", yi = cLM01_postPI(:);
        case "LC9",  yi = LC9_postPI(:);
        case "LT43", yi = LT43_postPI(:);
        case "LT52", yi = LT52_postPI(:);
        case "aMe1", yi = aMe1_postPI(:);
    end
    yi = yi(~isnan(yi));

    % 원자료 지터
    xj = i + 0.06*randn(numel(yi),1);
    scatter(xj, yi, 18, C(i,:), 'filled', 'MarkerFaceAlpha', 0.35, 'MarkerEdgeColor','none');

    % 평균 + SD 에러바
    n = numel(yi);
    if n>=2
        m  = mean(yi);
        sd = std(yi);
        errorbar(i, m, sd, 'Color', C(i,:), 'CapSize', 6, 'LineWidth', 1.2); % SD
        scatter(i, m, 36, C(i,:), 'filled', 'MarkerEdgeColor','k');          % mean 점
    elseif n==1
        scatter(i, yi, 36, C(i,:), 'filled', 'MarkerEdgeColor','k');
    end
end
ylim([0 1])
xlim([0.5 7.5])
title('postPI')
%% =========================
%  Crawford–Howell single-case t-tests (MAIN TEST ONLY)
%  - control: FF types (optionally FB too)
%  - two-sided, alpha = 0.05
%  =========================
alpha = 0.05;

% --- 1) 기준군 (FF) 준비
ctrl_FF = FF_postPI(:);
ctrl_FF = ctrl_FF(~isnan(ctrl_FF));
nFF = numel(ctrl_FF);
muFF = mean(ctrl_FF, 'omitnan');
sdFF = std(ctrl_FF, 'omitnan');

% (선택) FB를 기준군으로도 보고 싶다면 주석 해제
ctrl_FB = FB_postPI(:);
ctrl_FB = ctrl_FB(~isnan(ctrl_FB));
nFB = numel(ctrl_FB);
muFB = mean(ctrl_FB, 'omitnan');
sdFB = std(ctrl_FB, 'omitnan');

% --- 2) 케이스(관심 타입들)는 "타입 평균값"으로 통일
cases = struct( ...
    'name', {'LC9','LT43','LT52','aMe1','cLM01'}, ...
    'vals', {LC9_postPI(:), LT43_postPI(:), LT52_postPI(:), aMe1_postPI(:), cLM01_postPI(:)} );

case_names = strings(numel(cases),1);
x_case     = nan(numel(cases),1);

for i = 1:numel(cases)
    xi = cases(i).vals;
    xi = xi(~isnan(xi));
    case_names(i) = string(cases(i).name);
    x_case(i)     = mean(xi, 'omitnan');   % 타입 평균값 = 단일-사례 값
end

% --- 3) CH t-test 함수 (인라인)
ch_test = @(xcase, mu, sd, n) deal( ...
    (xcase - mu) ./ (sd .* sqrt((n+1)/n)), ...    % t
    n-1, ...                                       % df
    2*tcdf(-abs((xcase - mu) ./ (sd .* sqrt((n+1)/n))), n-1) ); % two-sided p

% --- 4) FF 기준군에 대해 계산/정리
if nFF < 3 || sdFF==0
    warning('FF 기준군 표본수/분산이 부족하여 CH-t 계산이 불안정할 수 있습니다. (n=%d, sd=%.3g)', nFF, sdFF);
end
[tFF, dfFF, pFF] = ch_test(x_case, muFF, sdFF, nFF);

T_FF = table(case_names, x_case, ...
             repmat(muFF,numel(cases),1), repmat(sdFF,numel(cases),1), ...
             repmat(nFF ,numel(cases),1), tFF, repmat(dfFF,numel(cases),1), pFF, ...
             'VariableNames', {'Case','CaseMean','CtrlMean','CtrlSD','CtrlN','t','df','p'});
T_FF.Control = repmat("FF types", height(T_FF), 1);
T_FF = movevars(T_FF, 'Control', 'Before', 'Case');

disp('=== Crawford–Howell single-case t-tests (two-sided, alpha=0.05) ===');
disp(T_FF);

% --- (선택) FB 기준군으로도 동시에 계산하려면 주석 해제
if nFB < 3 || sdFB==0
    warning('FB 기준군 표본수/분산이 부족할 수 있습니다. (n=%d, sd=%.3g)', nFB, sdFB);
end
[tFB, dfFB, pFB] = ch_test(x_case, muFB, sdFB, nFB);
T_FB = table(case_names, x_case, ...
             repmat(muFB,numel(cases),1), repmat(sdFB,numel(cases),1), ...
             repmat(nFB ,numel(cases),1), tFB, repmat(dfFB,numel(cases),1), pFB, ...
             'VariableNames', {'Case','CaseMean','CtrlMean','CtrlSD','CtrlN','t','df','p'});
T_FB.Control = repmat("FB types", height(T_FB), 1);
T_FB = movevars(T_FB, 'Control', 'Before', 'Case');
disp(T_FB);

%%
C = [
    20  97 154;   % #14619A
    111 151  51;  % #6F9733
    183  42  49;   % #B72A31  (LC9)
    152 122  63;   % #987A3F  (LT43)
    108  41 119;   % #6C2977  (LT52)
    43  170 184    % #2BAAB8  (aMe1)
    228 172  41;   % #E4AC29  (cLM01)
] / 255;
figure('Color','w'); clf; ax = gca; hold(ax,'on');
% 각 그룹: 원자료 지터 + 평균점 + SD 에러바
lbls=["FF";"FB"; "LC9";"LT43";"LT52";"aMe1";"cLM01"];
lbls_str = string(lbls);
for i = 1:numel(lbls)
    switch lbls_str(i)
        case "FF",   yi = FF_prePI(:);
        case "FB",   yi = FB_prePI(:);
        case "cLM01", yi = cLM01_prePI(:);
        case "LC9",  yi = LC9_prePI(:);
        case "LT43", yi = LT43_prePI(:);
        case "LT52", yi = LT52_prePI(:);
        case "aMe1", yi = aMe1_prePI(:);
    end
    yi = yi(~isnan(yi));

    % 원자료 지터
    xj = i + 0.06*randn(numel(yi),1);
    scatter(xj, yi, 18, C(i,:), 'filled', 'MarkerFaceAlpha', 0.35, 'MarkerEdgeColor','none');

    % 평균 + SD 에러바
    n = numel(yi);
    if n>=2
        m  = mean(yi);
        sd = std(yi);
        errorbar(i, m, sd, 'Color', C(i,:), 'CapSize', 6, 'LineWidth', 1.2); % SD
        scatter(i, m, 36, C(i,:), 'filled', 'MarkerEdgeColor','k');          % mean 점
    elseif n==1
        scatter(i, yi, 36, C(i,:), 'filled', 'MarkerEdgeColor','k');
    end
end
ylim([0 1])
xlim([0.5 7.5])
title('prePI')

%% =========================
%  Crawford–Howell single-case t-tests (MAIN TEST ONLY)
%  - control: FF types (optionally FB too)
%  - two-sided, alpha = 0.05
%  =========================
alpha = 0.05;

% --- 1) 기준군 (FF) 준비
ctrl_FF = FF_prePI(:);
ctrl_FF = ctrl_FF(~isnan(ctrl_FF));
nFF = numel(ctrl_FF);
muFF = mean(ctrl_FF, 'omitnan');
sdFF = std(ctrl_FF, 'omitnan');

% (선택) FB를 기준군으로도 보고 싶다면 주석 해제
ctrl_FB = FB_prePI(:);
ctrl_FB = ctrl_FB(~isnan(ctrl_FB));
nFB = numel(ctrl_FB);
muFB = mean(ctrl_FB, 'omitnan');
sdFB = std(ctrl_FB, 'omitnan');

% --- 2) 케이스(관심 타입들)는 "타입 평균값"으로 통일
cases = struct( ...
    'name', {'LC9','LT43','LT52','aMe1','cLM01'}, ...
    'vals', {LC9_prePI(:), LT43_prePI(:), LT52_prePI(:), aMe1_prePI(:), cLM01_prePI(:)} );

case_names = strings(numel(cases),1);
x_case     = nan(numel(cases),1);

for i = 1:numel(cases)
    xi = cases(i).vals;
    xi = xi(~isnan(xi));
    case_names(i) = string(cases(i).name);
    x_case(i)     = mean(xi, 'omitnan');   % 타입 평균값 = 단일-사례 값
end

% --- 3) CH t-test 함수 (인라인)
ch_test = @(xcase, mu, sd, n) deal( ...
    (xcase - mu) ./ (sd .* sqrt((n+1)/n)), ...    % t
    n-1, ...                                       % df
    2*tcdf(-abs((xcase - mu) ./ (sd .* sqrt((n+1)/n))), n-1) ); % two-sided p

% --- 4) FF 기준군에 대해 계산/정리
if nFF < 3 || sdFF==0
    warning('FF 기준군 표본수/분산이 부족하여 CH-t 계산이 불안정할 수 있습니다. (n=%d, sd=%.3g)', nFF, sdFF);
end
[tFF, dfFF, pFF] = ch_test(x_case, muFF, sdFF, nFF);

T_FF = table(case_names, x_case, ...
             repmat(muFF,numel(cases),1), repmat(sdFF,numel(cases),1), ...
             repmat(nFF ,numel(cases),1), tFF, repmat(dfFF,numel(cases),1), pFF, ...
             'VariableNames', {'Case','CaseMean','CtrlMean','CtrlSD','CtrlN','t','df','p'});
T_FF.Control = repmat("FF types", height(T_FF), 1);
T_FF = movevars(T_FF, 'Control', 'Before', 'Case');

disp('=== Crawford–Howell single-case t-tests (two-sided, alpha=0.05) ===');
disp(T_FF);

% --- (선택) FB 기준군으로도 동시에 계산하려면 주석 해제
if nFB < 3 || sdFB==0
    warning('FB 기준군 표본수/분산이 부족할 수 있습니다. (n=%d, sd=%.3g)', nFB, sdFB);
end
[tFB, dfFB, pFB] = ch_test(x_case, muFB, sdFB, nFB);
T_FB = table(case_names, x_case, ...
             repmat(muFB,numel(cases),1), repmat(sdFB,numel(cases),1), ...
             repmat(nFB ,numel(cases),1), tFB, repmat(dfFB,numel(cases),1), pFB, ...
             'VariableNames', {'Case','CaseMean','CtrlMean','CtrlSD','CtrlN','t','df','p'});
T_FB.Control = repmat("FB types", height(T_FB), 1);
T_FB = movevars(T_FB, 'Control', 'Before', 'Case');
disp(T_FB);


%% =====================================================================
%% =====================================================================
%  백분위수(Percentile) 분석 (단순 버전)
% ======================================================================

fprintf('\n\n--- [백분위수 분석 결과 (단순 버전)] ---\n');

%--------------------------------------------------------------------------
% 1. FF 그룹을 기준으로 분석
%--------------------------------------------------------------------------
fprintf('\n==================================================\n');
fprintf('   분석 기준: [FF 그룹] 분포\n');
fprintf('==================================================\n\n');

% --- [FF 기준] postPI (Input Ratio) 분석 ---
fprintf('--- postPI (Input Ratio) 분석 ---\n');
ctrl_FF = FF_postPI(~isnan(FF_postPI));

mean_val = mean(LC9_postPI, 'omitnan');
rank = sum(ctrl_FF <= mean_val) / numel(ctrl_FF) * 100;
fprintf('LC9 평균 postPI   -> FF 그룹의 %.2f percentile\n', rank);

mean_val = mean(LT43_postPI, 'omitnan');
rank = sum(ctrl_FF <= mean_val) / numel(ctrl_FF) * 100;
fprintf('LT43 평균 postPI  -> FF 그룹의 %.2f percentile\n', rank);

mean_val = mean(LT52_postPI, 'omitnan');
rank = sum(ctrl_FF <= mean_val) / numel(ctrl_FF) * 100;
fprintf('LT52 평균 postPI  -> FF 그룹의 %.2f percentile\n', rank);

mean_val = mean(aMe1_postPI, 'omitnan');
rank = sum(ctrl_FF <= mean_val) / numel(ctrl_FF) * 100;
fprintf('aMe1 평균 postPI  -> FF 그룹의 %.2f percentile\n', rank);

mean_val = mean(cLM01_postPI, 'omitnan');
rank = sum(ctrl_FF <= mean_val) / numel(ctrl_FF) * 100;
fprintf('cLM01 평균 postPI -> FF 그룹의 %.2f percentile\n', rank);
fprintf('\n');

% --- [FF 기준] prePI (Output Ratio) 분석 ---
fprintf('--- prePI (Output Ratio) 분석 ---\n');
ctrl_FF = FF_prePI(~isnan(FF_prePI));

mean_val = mean(LC9_prePI, 'omitnan');
rank = sum(ctrl_FF <= mean_val) / numel(ctrl_FF) * 100;
fprintf('LC9 평균 prePI   -> FF 그룹의 %.2f percentile\n', rank);

mean_val = mean(LT43_prePI, 'omitnan');
rank = sum(ctrl_FF <= mean_val) / numel(ctrl_FF) * 100;
fprintf('LT43 평균 prePI  -> FF 그룹의 %.2f percentile\n', rank);

mean_val = mean(LT52_prePI, 'omitnan');
rank = sum(ctrl_FF <= mean_val) / numel(ctrl_FF) * 100;
fprintf('LT52 평균 prePI  -> FF 그룹의 %.2f percentile\n', rank);

mean_val = mean(aMe1_prePI, 'omitnan');
rank = sum(ctrl_FF <= mean_val) / numel(ctrl_FF) * 100;
fprintf('aMe1 평균 prePI  -> FF 그룹의 %.2f percentile\n', rank);

mean_val = mean(cLM01_prePI, 'omitnan');
rank = sum(ctrl_FF <= mean_val) / numel(ctrl_FF) * 100;
fprintf('cLM01 평균 prePI -> FF 그룹의 %.2f percentile\n', rank);


%--------------------------------------------------------------------------
% 2. FB 그룹을 기준으로 분석
%--------------------------------------------------------------------------
fprintf('\n\n==================================================\n');
fprintf('   분석 기준: [FB 그룹] 분포\n');
fprintf('==================================================\n\n');

% --- [FB 기준] postPI (Input Ratio) 분석 ---
fprintf('--- postPI (Input Ratio) 분석 ---\n');
ctrl_FB = FB_postPI(~isnan(FB_postPI));

mean_val = mean(LC9_postPI, 'omitnan');
rank = sum(ctrl_FB <= mean_val) / numel(ctrl_FB) * 100;
fprintf('LC9 평균 postPI   -> FB 그룹의 %.2f percentile\n', rank);

mean_val = mean(LT43_postPI, 'omitnan');
rank = sum(ctrl_FB <= mean_val) / numel(ctrl_FB) * 100;
fprintf('LT43 평균 postPI  -> FB 그룹의 %.2f percentile\n', rank);

mean_val = mean(LT52_postPI, 'omitnan');
rank = sum(ctrl_FB <= mean_val) / numel(ctrl_FB) * 100;
fprintf('LT52 평균 postPI  -> FB 그룹의 %.2f percentile\n', rank);

mean_val = mean(aMe1_postPI, 'omitnan');
rank = sum(ctrl_FB <= mean_val) / numel(ctrl_FB) * 100;
fprintf('aMe1 평균 postPI  -> FB 그룹의 %.2f percentile\n', rank);

mean_val = mean(cLM01_postPI, 'omitnan');
rank = sum(ctrl_FB <= mean_val) / numel(ctrl_FB) * 100;
fprintf('cLM01 평균 postPI -> FB 그룹의 %.2f percentile\n', rank);
fprintf('\n');

% --- [FB 기준] prePI (Output Ratio) 분석 ---
fprintf('--- prePI (Output Ratio) 분석 ---\n');
ctrl_FB = FB_prePI(~isnan(FB_prePI));

mean_val = mean(LC9_prePI, 'omitnan');
rank = sum(ctrl_FB <= mean_val) / numel(ctrl_FB) * 100;
fprintf('LC9 평균 prePI   -> FB 그룹의 %.2f percentile\n', rank);

mean_val = mean(LT43_prePI, 'omitnan');
rank = sum(ctrl_FB <= mean_val) / numel(ctrl_FB) * 100;
fprintf('LT43 평균 prePI  -> FB 그룹의 %.2f percentile\n', rank);

mean_val = mean(LT52_prePI, 'omitnan');
rank = sum(ctrl_FB <= mean_val) / numel(ctrl_FB) * 100;
fprintf('LT52 평균 prePI  -> FB 그룹의 %.2f percentile\n', rank);

mean_val = mean(aMe1_prePI, 'omitnan');
rank = sum(ctrl_FB <= mean_val) / numel(ctrl_FB) * 100;
fprintf('aMe1 평균 prePI  -> FB 그룹의 %.2f percentile\n', rank);

mean_val = mean(cLM01_prePI, 'omitnan');
rank = sum(ctrl_FB <= mean_val) / numel(ctrl_FB) * 100;
fprintf('cLM01 평균 prePI -> FB 그룹의 %.2f percentile\n', rank);
fprintf('\n');
%%%%%

%% =====================================================================
%  Z-Score (1D Mahalanobis Distance) 분석
%  - 각 타입의 평균값이 FF/FB 그룹의 분포에서 얼마나 떨어져 있는지 계산
% ======================================================================
fprintf('\n\n--- [Z-Score / 1D Mahalanobis Distance 분석 결과] ---\n');

% --- 데이터 준비 (이미 Workspace에 변수가 있다면 이 부분은 생략 가능) ---
% FF/FB 그룹 데이터
ctrl_FF_post = G{1}(:,1);
ctrl_FF_pre  = G{1}(:,2);
ctrl_FB_post = G{2}(:,1);
ctrl_FB_pre  = G{2}(:,2);
% 분석할 케이스 그룹 데이터
cases = struct( ...
    'name', {'LC9','LT43','LT52','aMe1','cLM01'}, ...
    'postPI', {G{3}(:,1), G{4}(:,1), G{5}(:,1), G{6}(:,1), G{7}(:,1)}, ...
    'prePI',  {G{3}(:,2), G{4}(:,2), G{5}(:,2), G{6}(:,2), G{7}(:,2)} );

% --- Z-Score 계산 함수 정의 ---
% Z-점수의 절댓값이 곧 단변량 마할라노비스 거리 D 입니다.
z_score_func = @(x, mu, sigma) (x - mu) / sigma;

%--------------------------------------------------------------------------
% 1. FF 그룹을 기준으로 분석
%--------------------------------------------------------------------------
fprintf('\n==================================================\n');
fprintf('   분석 기준: [FF 그룹] 분포\n');
fprintf('==================================================\n\n');

% FF 그룹의 통계량 계산
mu_FF_post = mean(ctrl_FF_post, 'omitnan');
sd_FF_post = std(ctrl_FF_post, 'omitnan');
mu_FF_pre  = mean(ctrl_FF_pre, 'omitnan');
sd_FF_pre  = std(ctrl_FF_pre, 'omitnan');

% --- [FF 기준] postPI (Input Ratio) 분석 ---
fprintf('--- postPI (Input Ratio) Z-Scores ---\n');
for i = 1:numel(cases)
    case_mean = mean(cases(i).postPI, 'omitnan');
    z = z_score_func(case_mean, mu_FF_post, sd_FF_post);
    fprintf('%s 평균 postPI   -> FF 그룹 평균에서 Z = %+.4f (%.2f 표준편차 거리)\n', cases(i).name, z, abs(z));
end
fprintf('\n');

% --- [FF 기준] prePI (Output Ratio) 분석 ---
fprintf('--- prePI (Output Ratio) Z-Scores ---\n');
for i = 1:numel(cases)
    case_mean = mean(cases(i).prePI, 'omitnan');
    z = z_score_func(case_mean, mu_FF_pre, sd_FF_pre);
    fprintf('%s 평균 prePI    -> FF 그룹 평균에서 Z = %+.4f (%.2f 표준편차 거리)\n', cases(i).name, z, abs(z));
end

%--------------------------------------------------------------------------
% 2. FB 그룹을 기준으로 분석
%--------------------------------------------------------------------------
fprintf('\n\n==================================================\n');
fprintf('   분석 기준: [FB 그룹] 분포\n');
fprintf('==================================================\n\n');

% FB 그룹의 통계량 계산
mu_FB_post = mean(ctrl_FB_post, 'omitnan');
sd_FB_post = std(ctrl_FB_post, 'omitnan');
mu_FB_pre  = mean(ctrl_FB_pre, 'omitnan');
sd_FB_pre  = std(ctrl_FB_pre, 'omitnan');

% --- [FB 기준] postPI (Input Ratio) 분석 ---
fprintf('--- postPI (Input Ratio) Z-Scores ---\n');
for i = 1:numel(cases)
    case_mean = mean(cases(i).postPI, 'omitnan');
    z = z_score_func(case_mean, mu_FB_post, sd_FB_post);
    fprintf('%s 평균 postPI   -> FB 그룹 평균에서 Z = %+.4f (%.2f 표준편차 거리)\n', cases(i).name, z, abs(z));
end
fprintf('\n');

% --- [FB 기준] prePI (Output Ratio) 분석 ---
fprintf('--- prePI (Output Ratio) Z-Scores ---\n');
for i = 1:numel(cases)
    case_mean = mean(cases(i).prePI, 'omitnan');
    z = z_score_func(case_mean, mu_FB_pre, sd_FB_pre);
    fprintf('%s 평균 prePI    -> FB 그룹 평균에서 Z = %+.4f (%.2f 표준편차 거리)\n', cases(i).name, z, abs(z));
end

fprintf('\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
