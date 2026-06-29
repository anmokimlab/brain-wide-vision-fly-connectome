close all; clear vars; clc

load('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\Figure1_Overall\FAFB_NPI_Thr0.mat')


opt = detectImportOptions('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\Figure4_BD\Segregation\RightFF_NPIs_segregation_index.csv');
opt = setvartype(opt,'root_id','int64');
RightFF_NPIs = readtable('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\Figure4_BD\Segregation\RightFF_NPIs_segregation_index.csv',opt);
RightFF_NPIs = removevars(RightFF_NPIs, "Var1");


opt = detectImportOptions('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\Figure4_BD\Segregation\RightFB_NPIs_segregation_index.csv');
opt = setvartype(opt,'root_id','int64');
RightFB_NPIs = readtable('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\Figure4_BD\Segregation\RightFB_NPIs_segregation_index.csv',opt);
RightFB_NPIs = removevars(RightFB_NPIs, "Var1");


opt = detectImportOptions('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\Figure4_BD\Segregation\RightBD_NPIs_segregation_index.csv');
opt = setvartype(opt,'root_id','int64');
RightBD_NPIs = readtable('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\Figure4_BD\Segregation\RightBD_NPIs_segregation_index.csv',opt);
RightBD_NPIs = removevars(RightBD_NPIs, "Var1");

%% 여러 타입 입력
% 예: WantToSee = {'LC9','LC10b'};
WantToSee = {'LC9','LT43','LT52','aMe1','cLM01'};

% 각 타입별 right side root_id 수집
Want = struct('type',[],'root_ids',[],'tbl',[]);
all_want_root_ids = [];

for k = 1:numel(WantToSee)
    Want(k).type = WantToSee{k};
    rid = FAFBNPIs.root_id( strcmp(FAFBNPIs.type, Want(k).type) & (FAFBNPIs.In_Synapse_Optic_R>0) );
    Want(k).root_ids = rid(:);
    all_want_root_ids = [all_want_root_ids; rid(:)]; %#ok<AGROW>
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
Type_FF = table(FF_Type,'VariableNames',{'type'});
for i=1:size(Type_FF,1)
    idx = (ic==i);
    Type_FF.root_id{i} = FF_neuropils.root_id(idx);
    Type_FF.mean_Right_PostPI(i) = mean(FF_neuropils.Right_PostPI(idx),'omitnan');
    Type_FF.mean_Right_PrePI(i) = mean(FF_neuropils.Right_PrePI(idx),'omitnan');

    Type_FF.segregation_index{i}  = FF_neuropils.segregation_index(idx);
    Type_FF.mean_segregation_index(i)  = mean(FF_neuropils.segregation_index(idx),'omitnan');
end

%% 타입별 FB 요약 (in/out 비율)
[FB_Type,~,ic] = unique(FB_neuropils.type);
Type_FB = table(FB_Type,'VariableNames',{'type'});
for i=1:size(Type_FB,1)
    idx = (ic==i);
    Type_FB.root_id{i} = FB_neuropils.root_id(idx);
    Type_FB.mean_Right_PostPI(i) = mean(FB_neuropils.Right_PostPI(idx),'omitnan');
    Type_FB.mean_Right_PrePI(i) = mean(FB_neuropils.Right_PrePI(idx),'omitnan');
    Type_FB.segregation_index{i}  = FB_neuropils.segregation_index(idx);
    Type_FB.mean_segregation_index(i)  = mean(FB_neuropils.segregation_index(idx),'omitnan');
end
%%
BD_neuropils = RightBD_NPIs;
[BD_Type,~,ic] = unique(BD_neuropils.type);
Type_BD = table(BD_Type,'VariableNames',{'type'});
for i=1:size(Type_BD,1)
    idx = (ic==i);
    Type_BD.root_id{i} = BD_neuropils.root_id(idx);
    Type_BD.mean_Right_PostPI(i) = mean(BD_neuropils.Right_PostPI(idx),'omitnan');
    Type_BD.mean_Right_PrePI(i) = mean(BD_neuropils.Right_PrePI(idx),'omitnan');

    Type_BD.segregation_index{i}  = BD_neuropils.segregation_index(idx);
    Type_BD.mean_segregation_index(i)  = mean(BD_neuropils.segregation_index(idx),'omitnan');
end

%% Want 에 대해서 
for i=1:1:size(Want,2)
    root_ids=Want(i).root_ids;
    segregation_index=[];
    for j=1:1:size(root_ids,1)
        current_root_id=root_ids(j);
        idx=RightBD_NPIs.root_id==current_root_id;
        segregation_index(j,1)=RightBD_NPIs.segregation_index(idx); 

    end
    Want(i).segregation_index=segregation_index;
end

%%
labels = {'FF','FB'};  % 그룹 라벨
labels=[labels WantToSee];

FF_SegIdx=Type_FF.mean_segregation_index;
FB_SegIdx=Type_FB.mean_segregation_index;
LC9_SegIdx=Want(1).segregation_index;
LT43_SegIdx=Want(2).segregation_index;
LT52_SegIdx=Want(3).segregation_index;
aMe1_SegIdx=Want(4).segregation_index;
cLM01_SegIdx=Want(5).segregation_index;

%%
if exist('labels','var') && ~isempty(labels) && numel(labels)==7
    if isstring(labels)
        lbls = labels(:);
    elseif iscell(labels)
        lbls = string(labels(:));
    else
        % numeric/categorical도 호환되게 문자열로 변환
        lbls = string(labels(:));
    end
else
    lbls = unique(g, 'stable');   % ["FF","FB","cLM01","LC9","LT43","LT52","aMe1"] 형태의 string
end

% (그림 함수용) cellstr로 필요 시 변환
lbls_cell = cellstr(lbls);
C = [
    20  97 154;   % #14619A
    111 151  51;  % #6F9733
    183  42  49;  % #B72A31
    152 122  63;  % #987A3F
    108  41 119;  % #6C2977
    43 170 184   % #2BAAB8
    228 172  41;  % #E4AC29  ← 수정

] / 255;
figure('Color','w'); clf; ax = gca; hold(ax,'on');
% 각 그룹: 원자료 지터 + 평균점 + SD 에러바
lbls_str = string(lbls);
for i = 1:numel(lbls)
    switch lbls_str(i)
        case "FF",   yi = FF_SegIdx(:);
        case "FB",   yi = FB_SegIdx(:);
        case "cLM01", yi = cLM01_SegIdx(:);
        case "LC9",  yi = LC9_SegIdx(:);
        case "LT43", yi = LT43_SegIdx(:);
        case "LT52", yi = LT52_SegIdx(:);
        case "aMe1", yi = aMe1_SegIdx(:);
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
%% =========================
%  Crawford–Howell single-case t-tests (MAIN TEST ONLY)
%  - control: FF types (optionally FB too)
%  - two-sided, alpha = 0.05
%  =========================
alpha = 0.05;

% --- 1) 기준군 (FF) 준비
ctrl_FF = FF_SegIdx(:);
ctrl_FF = ctrl_FF(~isnan(ctrl_FF));
nFF = numel(ctrl_FF);
muFF = mean(ctrl_FF, 'omitnan');
sdFF = std(ctrl_FF, 'omitnan');

% (선택) FB를 기준군으로도 보고 싶다면 주석 해제
ctrl_FB = FB_SegIdx(:);
ctrl_FB = ctrl_FB(~isnan(ctrl_FB));
nFB = numel(ctrl_FB);
muFB = mean(ctrl_FB, 'omitnan');
sdFB = std(ctrl_FB, 'omitnan');

% --- 2) 케이스(관심 타입들)는 "타입 평균값"으로 통일
cases = struct( ...
    'name', {'LC9','LT43','LT52','aMe1','cLM01'}, ...
    'vals', {LC9_SegIdx(:), LT43_SegIdx(:), LT52_SegIdx(:), aMe1_SegIdx(:),cLM01_SegIdx(:)} );

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure; clf; hold on;

% 0) 데이터 준비 (Type_*의 mean 사용)
RightFF_PostPIMinustPrePI = (Type_FF.mean_Right_PostPI) - (Type_FF.mean_Right_PrePI);
RightFF_SegregationIdx    =  Type_FF.mean_segregation_index;

RightFB_PostPIMinustPrePI = (Type_FB.mean_Right_PostPI) - (Type_FB.mean_Right_PrePI);
RightFB_SegregationIdx    =  Type_FB.mean_segregation_index;

RightBD_PostPIMinustPrePI = (Type_BD.mean_Right_PostPI) - (Type_BD.mean_Right_PrePI);
RightBD_SegregationIdx    =  Type_BD.mean_segregation_index;

% 1) 그룹별 산점도
scatter(abs(RightFF_PostPIMinustPrePI), RightFF_SegregationIdx, 'filled', ...
    'MarkerFaceAlpha',0.5,'MarkerFaceColor',[0 0.4470 0.7410], ...
    'DisplayName',sprintf('Feedforward (n=%d)', numel(RightFF_SegregationIdx)));
scatter(abs(RightFB_PostPIMinustPrePI), RightFB_SegregationIdx, 'filled', ...
    'MarkerFaceAlpha',0.5,'MarkerFaceColor',[0.4660 0.6740 0.1880], ...
    'DisplayName',sprintf('Feedback (n=%d)', numel(RightFB_SegregationIdx)));
scatter(abs(RightBD_PostPIMinustPrePI), RightBD_SegregationIdx, 'filled', ...
    'MarkerFaceAlpha',0.5,'MarkerFaceColor',[0.65 0.65 0.65], ...
    'DisplayName',sprintf('Bidirectional (n=%d)', numel(RightBD_SegregationIdx)));

% 2) 전체 회귀 데이터 병합 (+ 그룹 레이블)
xAll = [abs(RightFF_PostPIMinustPrePI);
        abs(RightFB_PostPIMinustPrePI);
        abs(RightBD_PostPIMinustPrePI)];
yAll = [RightFF_SegregationIdx;
        RightFB_SegregationIdx;
        RightBD_SegregationIdx];

grp0 = [repmat("FF", numel(RightFF_SegregationIdx),1);
        repmat("FB", numel(RightFB_SegregationIdx),1);
        repmat("BD", numel(RightBD_SegregationIdx),1)];

% 3) NaN/Inf 제거 (grp에도 동일 마스크 적용!)
mask = isfinite(xAll) & isfinite(yAll);
xAll = xAll(mask);  yAll = yAll(mask);
grp  = grp0(mask);

% 4) 선형회귀 + 예측구간(PI)
mdl  = fitlm(xAll, yAll);                    % y = b0 + b1*x
xfit = linspace(min(xAll), max(xAll), 200)';

% 평균의 신뢰구간(CI)와 새 관측치의 예측구간(PI)
[yfit_CI, yCI] = predict(mdl, xfit);                            % CI (mean)
[~,        yPI] = predict(mdl, xfit, 'Prediction','observation');% PI (obs)

% 회귀선 + PI 밴드 (보고 목적엔 PI가 더 적합)
plot(xfit, yfit_CI, 'k-', 'LineWidth', 2, 'DisplayName','Linear fit (all)');
patch([xfit; flipud(xfit)], [yPI(:,1); flipud(yPI(:,2))], [0 0 0], ...
      'FaceAlpha', 0.06, 'EdgeColor','none', 'DisplayName','95% PI (new obs)');
% 필요시 CI도 보이려면 아래 주석 해제
% patch([xfit; flipud(xfit)], [yCI(:,1); flipud(yCI(:,2))], [0 0 0], ...
%       'FaceAlpha', 0.08, 'EdgeColor','none', 'DisplayName','95% CI (mean)');

% 5) 레이블/범례/스타일
xlabel('|Right PostPI - Right PrePI|');
ylabel('Segregation index');
set(gca, 'TickDir','out'); box off;
legend('Location','bestoutside');

% 6) 주석(계수/상관/유의성)
b0 = mdl.Coefficients.Estimate(1);
b1 = mdl.Coefficients.Estimate(2);
R2 = mdl.Rsquared.Ordinary;
[rp, p_r] = corr(xAll, yAll, 'Type','Pearson', 'Rows','complete');
p_slope   = mdl.Coefficients.pValue(2);

txt = sprintf('y = %.3f + %.3f x   R^2=%.3f   r=%.3f (p_r=%.3g)   p_{slope}=%.3g', ...
              b0, b1, R2, rp, p_r, p_slope);
text(min(xAll) + 0.02*range(xAll), max(yAll) - 0.05*range(yAll), txt, ...
     'FontSize',10, 'BackgroundColor','w','Margin',2);

cLM01_idx=strcmp(Type_BD.type,'cLM01');
LC9_idx=strcmp(Type_BD.type,'LC9');
LT43_idx=strcmp(Type_BD.type,'LT43');
LT52_idx=strcmp(Type_BD.type,'LT52');
aMe1_idx=strcmp(Type_BD.type,'aMe1');


scatter(abs(Type_BD.mean_Right_PostPI(LC9_idx)-Type_BD.mean_Right_PrePI(LC9_idx)), Type_BD.mean_segregation_index(LC9_idx), 'filled', ...
    'MarkerFaceAlpha',1,'MarkerFaceColor',C(3,:));
scatter(abs(Type_BD.mean_Right_PostPI(LT43_idx)-Type_BD.mean_Right_PrePI(LT43_idx)), Type_BD.mean_segregation_index(LT43_idx), 'filled', ...
    'MarkerFaceAlpha',1,'MarkerFaceColor',C(4,:));
scatter(abs(Type_BD.mean_Right_PostPI(LT52_idx)-Type_BD.mean_Right_PrePI(LT52_idx)), Type_BD.mean_segregation_index(LT52_idx), 'filled', ...
    'MarkerFaceAlpha',1,'MarkerFaceColor',C(5,:));
scatter(abs(Type_BD.mean_Right_PostPI(aMe1_idx)-Type_BD.mean_Right_PrePI(aMe1_idx)), Type_BD.mean_segregation_index(aMe1_idx), 'filled', ...
    'MarkerFaceAlpha',1,'MarkerFaceColor',C(6,:));
scatter(abs(Type_BD.mean_Right_PostPI(cLM01_idx)-Type_BD.mean_Right_PrePI(cLM01_idx)), Type_BD.mean_segregation_index(cLM01_idx), 'filled', ...
    'MarkerFaceAlpha',1,'MarkerFaceColor',C(7,:));
hold off;
xlim([0 2])
ylim([0 1])
%%
% ====== 1) 중심/산포(강건) 추정: mad(...) 대신 직접 계산 ======
ctrl = FF_SegIdx(:);
ctrl = ctrl(isfinite(ctrl));       % NaN/Inf 제거

if isempty(ctrl)
    error('FF_SegIdx에서 유효한 값이 없습니다.');
end

% 중앙값
ctrl_med = median(ctrl);           % omitnan 불필요: 이미 NaN 제거

% 강건 MAD (median absolute deviation) → 표준정규 일치 스케일 1.4826 곱
mad_raw    = median(abs(ctrl - ctrl_med));
ctrl_MADn  = 1.4826 * mad_raw;

% MADn이 0이면(IQR 기반 대체) — 값이 거의 동일하거나 분포가 너무 뾰족한 경우
if ctrl_MADn == 0
    ctrl_MADn = iqr(ctrl) / 1.349;     % IQR → σ 강건근사
    warning('FF 기준분포의 MADn이 0입니다. IQR/1.349 대체치를 사용합니다.');
    if ctrl_MADn == 0
        ctrl_MADn = eps;               % 최후 안전장치
        warning('산포가 사실상 0입니다. RobustZ 해석에 주의하세요.');
    end
end

% ====== 2) 강건 z 점수 ======
RobustZ = (x_case - ctrl_med) ./ ctrl_MADn;

% ====== 3) 경험적 퍼센타일(중복 mid-rank) ======
sorted_ctrl = sort(ctrl);
n_ctrl = numel(sorted_ctrl);

percentile_midrank = @(x) ( sum(sorted_ctrl < x) + 0.5*sum(sorted_ctrl == x) ) / n_ctrl;

Percentile       = arrayfun(percentile_midrank, x_case) * 100;  % 0–100(%)
RightTailPct     = 100 - Percentile;
TwoSidedTailPct  = 2 * min(Percentile, RightTailPct);

% ====== 4) 테이블 정리 ======
T_MAD = table( ...
    case_names, x_case, ...
    repmat(ctrl_med, numel(x_case),1), repmat(ctrl_MADn, numel(x_case),1), ...
    RobustZ, Percentile, RightTailPct, TwoSidedTailPct, ...
    'VariableNames', {'Case','CaseMean','CtrlMedian','MADn','RobustZ','Percentile','RightTailPct','TwoSidedTailPct'});

T_MAD.PLabel = arrayfun(@(p) sprintf('%.1fp', p), T_MAD.Percentile, 'UniformOutput', false);
T_MAD.IsOutlierCandidate = (abs(T_MAD.RobustZ) >= 2.5) | (T_MAD.Percentile >= 99) | (T_MAD.Percentile <= 1);

disp('=== MAD 기반 퍼센타일 / Robust Z (FF 기준; NaN 제거, 직접계산) ===');
disp(T_MAD);
%%
% ====== 1) 중심/산포(강건) 추정: mad(...) 대신 직접 계산 ======
ctrl = FB_SegIdx(:);
ctrl = ctrl(isfinite(ctrl));       % NaN/Inf 제거

if isempty(ctrl)
    error('FF_SegIdx에서 유효한 값이 없습니다.');
end

% 중앙값
ctrl_med = median(ctrl);           % omitnan 불필요: 이미 NaN 제거

% 강건 MAD (median absolute deviation) → 표준정규 일치 스케일 1.4826 곱
mad_raw    = median(abs(ctrl - ctrl_med));
ctrl_MADn  = 1.4826 * mad_raw;

% MADn이 0이면(IQR 기반 대체) — 값이 거의 동일하거나 분포가 너무 뾰족한 경우
if ctrl_MADn == 0
    ctrl_MADn = iqr(ctrl) / 1.349;     % IQR → σ 강건근사
    warning('FF 기준분포의 MADn이 0입니다. IQR/1.349 대체치를 사용합니다.');
    if ctrl_MADn == 0
        ctrl_MADn = eps;               % 최후 안전장치
        warning('산포가 사실상 0입니다. RobustZ 해석에 주의하세요.');
    end
end

% ====== 2) 강건 z 점수 ======
RobustZ = (x_case - ctrl_med) ./ ctrl_MADn;

% ====== 3) 경험적 퍼센타일(중복 mid-rank) ======
sorted_ctrl = sort(ctrl);
n_ctrl = numel(sorted_ctrl);

percentile_midrank = @(x) ( sum(sorted_ctrl < x) + 0.5*sum(sorted_ctrl == x) ) / n_ctrl;

Percentile       = arrayfun(percentile_midrank, x_case) * 100;  % 0–100(%)
RightTailPct     = 100 - Percentile;
TwoSidedTailPct  = 2 * min(Percentile, RightTailPct);

% ====== 4) 테이블 정리 ======
T_MAD = table( ...
    case_names, x_case, ...
    repmat(ctrl_med, numel(x_case),1), repmat(ctrl_MADn, numel(x_case),1), ...
    RobustZ, Percentile, RightTailPct, TwoSidedTailPct, ...
    'VariableNames', {'Case','CaseMean','CtrlMedian','MADn','RobustZ','Percentile','RightTailPct','TwoSidedTailPct'});

T_MAD.PLabel = arrayfun(@(p) sprintf('%.1fp', p), T_MAD.Percentile, 'UniformOutput', false);
T_MAD.IsOutlierCandidate = (abs(T_MAD.RobustZ) >= 2.5) | (T_MAD.Percentile >= 99) | (T_MAD.Percentile <= 1);

disp('=== MAD 기반 퍼센타일 / Robust Z (FB 기준; NaN 제거, 직접계산) ===');
disp(T_MAD);

%% =======================================================
%  Z-score Analysis (vs FF control)
%  - Compares each case mean to the distribution of FF type means.
%  - Calculates standard z-score and two-sided p-value.
%  =======================================================

% --- 1) Prepare FF control group data
ctrl_FF = FF_SegIdx(:);
ctrl_FF = ctrl_FF(isfinite(ctrl_FF)); % Remove NaNs and Infs

if numel(ctrl_FF) < 2
    error('Not enough data in the FF control group to perform z-score analysis.');
end

% --- 2) Calculate control group statistics
mu_FF = mean(ctrl_FF);
sd_FF = std(ctrl_FF);

if sd_FF == 0
    warning('Standard deviation of the FF control group is zero. Z-scores will be Inf or NaN.');
end

% --- 3) Calculate Z-scores and p-values for each case
z_scores_FF = (x_case - mu_FF) ./ sd_FF;
p_values_FF = 2 * normcdf(-abs(z_scores_FF)); % Two-sided p-value

% --- 4) Display results in a table
T_Z_FF = table(case_names, x_case, ...
    repmat(mu_FF, numel(x_case), 1), repmat(sd_FF, numel(x_case), 1), ...
    z_scores_FF, p_values_FF, ...
    'VariableNames', {'Case', 'CaseMean', 'CtrlMean', 'CtrlSD', 'Z_score', 'P_value_2sided'});

disp('=== Z-score Analysis (Control: FF types) ===');
disp(T_Z_FF);

%% =======================================================
%  Z-score 분석 (FF 유형 그룹을 대조군으로)
% ========================================================

% --- 1) FF 대조군 데이터 준비
ctrl_FF = FF_SegIdx(:);
ctrl_FF = ctrl_FF(isfinite(ctrl_FF)); % NaN, Inf 등 유효하지 않은 값 제거

% --- 2) 대조군 통계량 계산 (평균, 표준편차)
mu_FF = mean(ctrl_FF);
sd_FF = std(ctrl_FF);

% 표준편차가 0일 경우 경고 메시지 출력
if sd_FF == 0
    warning('FF 대조군의 표준편차가 0입니다. Z-score가 Inf 또는 NaN이 될 수 있습니다.');
end

% --- 3) 각 케이스(뉴런 타입)의 Z-score 계산
z_scores_FF = (x_case - mu_FF) ./ sd_FF;

% --- 4) 결과를 테이블로 정리하여 출력
T_Z_FF = table(case_names, x_case, ...
    repmat(mu_FF, numel(x_case), 1), repmat(sd_FF, numel(x_case), 1), ...
    z_scores_FF, ...
    'VariableNames', {'케이스', '케이스_평균값', '대조군_평균', '대조군_표준편차', 'Z_score'});

disp('=== Z-score 분석 (대조군: FF 유형) ===');
disp(T_Z_FF);

%% =======================================================
%  Z-score 분석 (FB 유형 그룹을 대조군으로)
% ========================================================

% --- 1) FB 대조군 데이터 준비
ctrl_FB = FB_SegIdx(:);
ctrl_FB = ctrl_FB(isfinite(ctrl_FB)); % NaN, Inf 등 유효하지 않은 값 제거

% --- 2) 대조군 통계량 계산 (평균, 표준편차)
mu_FB = mean(ctrl_FB);
sd_FB = std(ctrl_FB);

% 표준편차가 0일 경우 경고 메시지 출력
if sd_FB == 0
    warning('FB 대조군의 표준편차가 0입니다. Z-score가 Inf 또는 NaN이 될 수 있습니다.');
end

% --- 3) 각 케이스(뉴런 타입)의 Z-score 계산
z_scores_FB = (x_case - mu_FB) ./ sd_FB;

% --- 4) 결과를 테이블로 정리하여 출력
T_Z_FB = table(case_names, x_case, ...
    repmat(mu_FB, numel(x_case), 1), repmat(sd_FB, numel(x_case), 1), ...
    z_scores_FB, ...
    'VariableNames', {'케이스', '케이스_평균값', '대조군_평균', '대조군_표준편차', 'Z_score'});

disp('=== Z-score 분석 (대조군: FB 유형) ===');
disp(T_Z_FB);