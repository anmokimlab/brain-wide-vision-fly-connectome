%% Clean

clear all; clc; close all
load('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\Figure1_Overall\Right_Neurons_Thr0.mat')
load('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\Figure1_Overall\FAFB_NPI_Thr0.mat')

opt = detectImportOptions('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\connections_no_threshold.csv');
opt = setvartype(opt,'pre_root_id','int64');
opt = setvartype(opt,'post_root_id','int64');
FAFBConnections = readtable('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\connections_no_threshold.csv',opt);

opt = detectImportOptions('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\consolidated_cell_types.csv');
opt = setvartype(opt,'root_id','int64');
FAFBConsolidated_type = readtable('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\consolidated_cell_types.csv',opt);
load('Reciprocal_Data.mat')
% 
% %% 여러 타입 입력
% % 예: WantToSee = {'LC9','LC10b'};
WantToSee = {'LC9','LT43','LT52','aMe1','cLM01'};
% 
% % % 각 타입별 right side root_id 수집
% Want = struct('type',[],'root_ids',[],'tbl',[]);
% 
% for k = 1:numel(WantToSee)
%     Want(k).type = WantToSee{k};
%     rid = FAFBNPIs.root_id( strcmp(FAFBNPIs.type, Want(k).type) & (FAFBNPIs.In_Synapse_Optic_R>0) );
%     Want(k).root_ids = rid(:);
% end

%% FF/FB 비교군에서 Wanted root_id 제거
% FF_neuropils = RightFF_NPIs;
% 
% 
% FB_neuropils = RightFB_NPIs;
% 
% 
% %%
% for i=1:1:size(FF_neuropils,1)
%     current_root_id=FF_neuropils.root_id(i);
% 
%     CurrentConnections_In=FAFBConnections(ismember(FAFBConnections.post_root_id,current_root_id),:);
% 
%     CurrentConnections_Out=FAFBConnections(ismember(FAFBConnections.pre_root_id,current_root_id),:);
% 
%     [InputNeurons,~,ic_in]=unique(CurrentConnections_In.pre_root_id);
% 
%     [Lia,Locb]=ismember(InputNeurons,FAFBConsolidated_type.root_id);
%     InputNeuronTypes=FAFBConsolidated_type.primary_type(Locb(Lia));
% 
% 
%     [OutputNeurons,~,ic_out]=unique(CurrentConnections_Out.post_root_id);
% 
%     [Lia,Locb]=ismember(OutputNeurons,FAFBConsolidated_type.root_id);
%     OutputNeuronTypes=FAFBConsolidated_type.primary_type(Locb(Lia));
% 
%     for j=1:1:size(InputNeurons,1)
%         idx=ic_in==j;
%         InputNeurons(j,2)=sum(CurrentConnections_In.syn_count(idx));
%     end
% 
%     for j=1:1:size(OutputNeurons,1)
%         idx=ic_out==j;
%         OutputNeurons(j,2)=sum(CurrentConnections_Out.syn_count(idx));
%     end
% 
%     [uniqueInputTypes,~,ic_in_types]=unique(InputNeuronTypes);
%     for j=1:1:size(uniqueInputTypes,1)
%         idx=ic_in_types==j;
%         uniqueInputTypes{j,2}=sum(InputNeurons(idx,2));
%     end
%     [uniqueOutputTypes,~,ic_out_types]=unique(OutputNeuronTypes);
%     for j=1:1:size(uniqueOutputTypes,1)
%         idx=ic_out_types==j;
%         uniqueOutputTypes{j,2}=sum(OutputNeurons(idx,2));
%     end
% 
%     [J, Jw] = jaccard_calculate(InputNeurons, OutputNeurons);
%     [J_type, Jw_type] = jaccard_calculate_type(uniqueInputTypes, uniqueOutputTypes);
%     FF_neuropils.jaccard(i)=J;
%     FF_neuropils.weighted_Jaccard(i)=Jw;
%     FF_neuropils.jaccard_type(i)=J_type;
%     FF_neuropils.weighted_Jaccard_type(i)=Jw_type;
% end
% 
% %%
% for i=1:1:size(FB_neuropils,1)
%     current_root_id=FB_neuropils.root_id(i);
% 
%     CurrentConnections_In=FAFBConnections(ismember(FAFBConnections.post_root_id,current_root_id),:);
% 
%     CurrentConnections_Out=FAFBConnections(ismember(FAFBConnections.pre_root_id,current_root_id),:);
% 
%     [InputNeurons,~,ic_in]=unique(CurrentConnections_In.pre_root_id);
% 
%     [Lia,Locb]=ismember(InputNeurons,FAFBConsolidated_type.root_id);
%     InputNeuronTypes=FAFBConsolidated_type.primary_type(Locb(Lia));
% 
% 
%     [OutputNeurons,~,ic_out]=unique(CurrentConnections_Out.post_root_id);
% 
%     [Lia,Locb]=ismember(OutputNeurons,FAFBConsolidated_type.root_id);
%     OutputNeuronTypes=FAFBConsolidated_type.primary_type(Locb(Lia));
% 
%     for j=1:1:size(InputNeurons,1)
%         idx=ic_in==j;
%         InputNeurons(j,2)=sum(CurrentConnections_In.syn_count(idx));
%     end
% 
%     for j=1:1:size(OutputNeurons,1)
%         idx=ic_out==j;
%         OutputNeurons(j,2)=sum(CurrentConnections_Out.syn_count(idx));
%     end
% 
%     [uniqueInputTypes,~,ic_in_types]=unique(InputNeuronTypes);
%     for j=1:1:size(uniqueInputTypes,1)
%         idx=ic_in_types==j;
%         uniqueInputTypes{j,2}=sum(InputNeurons(idx,2));
%     end
%     [uniqueOutputTypes,~,ic_out_types]=unique(OutputNeuronTypes);
%     for j=1:1:size(uniqueOutputTypes,1)
%         idx=ic_out_types==j;
%         uniqueOutputTypes{j,2}=sum(OutputNeurons(idx,2));
%     end
% 
%     [J, Jw] = jaccard_calculate(InputNeurons, OutputNeurons);
%     [J_type, Jw_type] = jaccard_calculate_type(uniqueInputTypes, uniqueOutputTypes);
%     FB_neuropils.jaccard(i)=J;
%     FB_neuropils.weighted_Jaccard(i)=Jw;
%     FB_neuropils.jaccard_type(i)=J_type;
%     FB_neuropils.weighted_Jaccard_type(i)=Jw_type;
% end
% 
% %%
% BD_neuropils = RightBD_NPIs;
% 
% 
% for i=1:1:size(BD_neuropils,1)
%     current_root_id=BD_neuropils.root_id(i);
% 
%     CurrentConnections_In=FAFBConnections(ismember(FAFBConnections.post_root_id,current_root_id),:);
% 
%     CurrentConnections_Out=FAFBConnections(ismember(FAFBConnections.pre_root_id,current_root_id),:);
% 
%     [InputNeurons,~,ic_in]=unique(CurrentConnections_In.pre_root_id);
% 
%     [Lia,Locb]=ismember(InputNeurons,FAFBConsolidated_type.root_id);
%     InputNeuronTypes=FAFBConsolidated_type.primary_type(Locb(Lia));
% 
% 
%     [OutputNeurons,~,ic_out]=unique(CurrentConnections_Out.post_root_id);
% 
%     [Lia,Locb]=ismember(OutputNeurons,FAFBConsolidated_type.root_id);
%     OutputNeuronTypes=FAFBConsolidated_type.primary_type(Locb(Lia));
% 
%     for j=1:1:size(InputNeurons,1)
%         idx=ic_in==j;
%         InputNeurons(j,2)=sum(CurrentConnections_In.syn_count(idx));
%     end
% 
%     for j=1:1:size(OutputNeurons,1)
%         idx=ic_out==j;
%         OutputNeurons(j,2)=sum(CurrentConnections_Out.syn_count(idx));
%     end
% 
%     [uniqueInputTypes,~,ic_in_types]=unique(InputNeuronTypes);
%     for j=1:1:size(uniqueInputTypes,1)
%         idx=ic_in_types==j;
%         uniqueInputTypes{j,2}=sum(InputNeurons(idx,2));
%     end
%     [uniqueOutputTypes,~,ic_out_types]=unique(OutputNeuronTypes);
%     for j=1:1:size(uniqueOutputTypes,1)
%         idx=ic_out_types==j;
%         uniqueOutputTypes{j,2}=sum(OutputNeurons(idx,2));
%     end
% 
%     [J, Jw] = jaccard_calculate(InputNeurons, OutputNeurons);
%     [J_type, Jw_type] = jaccard_calculate_type(uniqueInputTypes, uniqueOutputTypes);
%     BD_neuropils.jaccard(i)=J;
%     BD_neuropils.weighted_Jaccard(i)=Jw;
%     BD_neuropils.jaccard_type(i)=J_type;
%     BD_neuropils.weighted_Jaccard_type(i)=Jw_type;
% end
% 
% %% 타입별 FF 요약 (in/out 비율)
% [FF_Type,~,ic] = unique(FF_neuropils.type);
% Type_FF = table(FF_Type,'VariableNames',{'type'});
% for i=1:size(Type_FF,1)
%     idx = (ic==i);
%     Type_FF.root_id{i} = FF_neuropils.root_id(idx);
%     Type_FF.mean_Right_PostPI(i) = mean(FF_neuropils.Right_PostPI(idx),'omitnan');
%     Type_FF.mean_Right_PrePI(i) = mean(FF_neuropils.Right_PrePI(idx),'omitnan');
% 
%     Type_FF.jaccard(i)  = mean(FF_neuropils.jaccard(idx));
%     Type_FF.weighted_Jaccard(i)  = mean(FF_neuropils.weighted_Jaccard(idx));
% 
%     Type_FF.jaccard_type(i)  = mean(FF_neuropils.jaccard_type(idx));
%     Type_FF.weighted_Jaccard_type(i)  = mean(FF_neuropils.weighted_Jaccard_type(idx));
% end
% 
% %% 타입별 FB 요약 (in/out 비율)
% [FB_Type,~,ic] = unique(FB_neuropils.type);
% Type_FB = table(FB_Type,'VariableNames',{'type'});
% for i=1:size(Type_FB,1)
%     idx = (ic==i);
%     Type_FB.root_id{i} = FB_neuropils.root_id(idx);
%     Type_FB.mean_Right_PostPI(i) = mean(FB_neuropils.Right_PostPI(idx),'omitnan');
%     Type_FB.mean_Right_PrePI(i) = mean(FB_neuropils.Right_PrePI(idx),'omitnan');
% 
%     Type_FB.jaccard(i)  = mean(FB_neuropils.jaccard(idx));
%     Type_FB.weighted_Jaccard(i)  = mean(FB_neuropils.weighted_Jaccard(idx));
% 
%     Type_FB.jaccard_type(i) = mean(FB_neuropils.jaccard_type(idx));
%     Type_FB.weighted_Jaccard_type(i)  = mean(FB_neuropils.weighted_Jaccard_type(idx));
% end
% 
% %%
% [BD_Type,~,ic] = unique(BD_neuropils.type);
% Type_BD = table(BD_Type,'VariableNames',{'type'});
% for i=1:size(Type_BD,1)
%     idx = (ic==i);
%     Type_BD.root_id{i} = BD_neuropils.root_id(idx);
%     Type_BD.mean_Right_PostPI(i) = mean(BD_neuropils.Right_PostPI(idx),'omitnan');
%     Type_BD.mean_Right_PrePI(i) = mean(BD_neuropils.Right_PrePI(idx),'omitnan');
% 
%     Type_BD.jaccard(i)  = mean(BD_neuropils.jaccard(idx));
%     Type_BD.weighted_Jaccard(i)  = mean(BD_neuropils.weighted_Jaccard(idx));
% 
%     Type_BD.jaccard_type(i)  = mean(BD_neuropils.jaccard_type(idx));
%     Type_BD.weighted_Jaccard_type(i)  = mean(BD_neuropils.weighted_Jaccard_type(idx));
% end
% 
% % %%
% for i=1:1:size(Want,2)
%     root_ids=Want(i).root_ids;
%     jaccard=[];
%     weighted_Jaccard=[];
%     jaccard_type=[];
%     weighted_Jaccard_type=[];
%     for j=1:1:size(root_ids,1)
%         current_root_id=root_ids(j);
% 
% 
%         CurrentConnections_In=FAFBConnections(ismember(FAFBConnections.post_root_id,current_root_id),:);
% 
%         CurrentConnections_Out=FAFBConnections(ismember(FAFBConnections.pre_root_id,current_root_id),:);
% 
%         [InputNeurons,~,ic_in]=unique(CurrentConnections_In.pre_root_id);
% 
%         [Lia,Locb]=ismember(InputNeurons,FAFBConsolidated_type.root_id);
%         InputNeuronTypes=FAFBConsolidated_type.primary_type(Locb(Lia));
% 
% 
%         [OutputNeurons,~,ic_out]=unique(CurrentConnections_Out.post_root_id);
% 
%         [Lia,Locb]=ismember(OutputNeurons,FAFBConsolidated_type.root_id);
%         OutputNeuronTypes=FAFBConsolidated_type.primary_type(Locb(Lia));
% 
%         for k=1:1:size(InputNeurons,1)
%             idx=ic_in==k;
%             InputNeurons(k,2)=sum(CurrentConnections_In.syn_count(idx));
%         end
% 
%         for k=1:1:size(OutputNeurons,1)
%             idx=ic_out==k;
%             OutputNeurons(k,2)=sum(CurrentConnections_Out.syn_count(idx));
%         end
% 
%         [uniqueInputTypes,~,ic_in_types]=unique(InputNeuronTypes);
%         for k=1:1:size(uniqueInputTypes,1)
%             idx=ic_in_types==k;
%             uniqueInputTypes{k,2}=sum(InputNeurons(idx,2));
%         end
%         [uniqueOutputTypes,~,ic_out_types]=unique(OutputNeuronTypes);
%         for k=1:1:size(uniqueOutputTypes,1)
%             idx=ic_out_types==k;
%             uniqueOutputTypes{k,2}=sum(OutputNeurons(idx,2));
%         end
% 
%         [J, Jw] = jaccard_calculate(InputNeurons, OutputNeurons);
%         [J_type, Jw_type] = jaccard_calculate_type(uniqueInputTypes, uniqueOutputTypes);
% 
% 
%         jaccard(j,1)     = J;    % 양방향 평균(Chamfer)
%         weighted_Jaccard(j,1)=Jw;
%         jaccard_type(j,1)     = J_type;    % 양방향 평균(Chamfer)
%         weighted_Jaccard_type(j,1)=Jw_type;
%     end
%     Want(i).jaccard=jaccard;
%     Want(i).weighted_Jaccard=weighted_Jaccard;
%     Want(i).jaccard_type=jaccard_type;
%     Want(i).weighted_Jaccard_type=weighted_Jaccard_type;
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
labels = {'FF','FB'};  % 그룹 라벨
labels=[labels WantToSee];

FF_Reciprocal=Type_FF.weighted_Jaccard;
FB_Reciprocal=Type_FB.weighted_Jaccard;
LC9_Reciprocal=Want(1).weighted_Jaccard;
LT43_Reciprocal=Want(2).weighted_Jaccard;
LT52_Reciprocal=Want(3).weighted_Jaccard;
aMe1_Reciprocal=Want(4).weighted_Jaccard;
cLM01_Reciprocal=Want(5).weighted_Jaccard;

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
        case "FF",   yi = FF_Reciprocal(:);
        case "FB",   yi = FB_Reciprocal(:);
        case "cLM01", yi = cLM01_Reciprocal(:);
        case "LC9",  yi = LC9_Reciprocal(:);
        case "LT43", yi = LT43_Reciprocal(:);
        case "LT52", yi = LT52_Reciprocal(:);
        case "aMe1", yi = aMe1_Reciprocal(:);
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
set(gca,'TickDir','out')
ylim([0 0.3])
xlim([0.5 7.5])
%% =========================
%  Crawford–Howell single-case t-tests (MAIN TEST ONLY)
%  - control: FF types (optionally FB too)
%  - two-sided, alpha = 0.05
%  =========================
% --- 1) 기준군 (FF) 준비
ctrl_FF = FF_Reciprocal(:);
ctrl_FF = ctrl_FF(~isnan(ctrl_FF));
nFF = numel(ctrl_FF);
muFF = mean(ctrl_FF, 'omitnan');
sdFF = std(ctrl_FF, 'omitnan');

% (선택) FB를 기준군으로도 보고 싶다면 주석 해제
ctrl_FB = FB_Reciprocal(:);
ctrl_FB = ctrl_FB(~isnan(ctrl_FB));
nFB = numel(ctrl_FB);
muFB = mean(ctrl_FB, 'omitnan');
sdFB = std(ctrl_FB, 'omitnan');

% --- 2) 케이스(관심 타입들)는 "타입 평균값"으로 통일
cases = struct( ...
    'name', {'LC9','LT43','LT52','aMe1','cLM01'}, ...
    'vals', {LC9_Reciprocal(:), LT43_Reciprocal(:), LT52_Reciprocal(:), aMe1_Reciprocal(:),cLM01_Reciprocal(:)} );

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
%%
figure; clf; hold on;

% 0) 데이터 준비 (Type_*의 mean 사용)
RightFF_PostPIMinustPrePI = (Type_FF.mean_Right_PostPI) - (Type_FF.mean_Right_PrePI);
RightFF_Reciprocal    =  Type_FF.weighted_Jaccard;

RightFB_PostPIMinustPrePI = (Type_FB.mean_Right_PostPI) - (Type_FB.mean_Right_PrePI);
RightFB_Reciprocal    =  Type_FB.weighted_Jaccard;

RightBD_PostPIMinustPrePI = (Type_BD.mean_Right_PostPI) - (Type_BD.mean_Right_PrePI);
RightBD_Reciprocal    =  Type_BD.weighted_Jaccard;

% 1) 그룹별 산점도
scatter(abs(RightFF_PostPIMinustPrePI), RightFF_Reciprocal, 'filled', ...
    'MarkerFaceAlpha',0.5,'MarkerFaceColor',[0 0.4470 0.7410], ...
    'DisplayName',sprintf('Feedforward (n=%d)', numel(RightFF_Reciprocal)));
scatter(abs(RightFB_PostPIMinustPrePI), RightFB_Reciprocal, 'filled', ...
    'MarkerFaceAlpha',0.5,'MarkerFaceColor',[0.4660 0.6740 0.1880], ...
    'DisplayName',sprintf('Feedback (n=%d)', numel(RightFB_Reciprocal)));
scatter(abs(RightBD_PostPIMinustPrePI), RightBD_Reciprocal, 'filled', ...
    'MarkerFaceAlpha',0.5,'MarkerFaceColor',[0.65 0.65 0.65], ...
    'DisplayName',sprintf('Bidirectional (n=%d)', numel(RightBD_Reciprocal)));

% 2) 전체 회귀 데이터 병합 (+ 그룹 레이블)
xAll = [abs(RightFF_PostPIMinustPrePI);
        abs(RightFB_PostPIMinustPrePI);
        abs(RightBD_PostPIMinustPrePI)];
yAll = [RightFF_Reciprocal;
        RightFB_Reciprocal;
        RightBD_Reciprocal];

grp0 = [repmat("FF", numel(RightFF_Reciprocal),1);
        repmat("FB", numel(RightFB_Reciprocal),1);
        repmat("BD", numel(RightBD_Reciprocal),1)];

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

scatter(abs(Type_BD.mean_Right_PostPI(LC9_idx)-Type_BD.mean_Right_PrePI(LC9_idx)), Type_BD.weighted_Jaccard(LC9_idx), 'filled', ...
    'MarkerFaceAlpha',1,'MarkerFaceColor',C(3,:));
scatter(abs(Type_BD.mean_Right_PostPI(LT43_idx)-Type_BD.mean_Right_PrePI(LT43_idx)), Type_BD.weighted_Jaccard(LT43_idx), 'filled', ...
    'MarkerFaceAlpha',1,'MarkerFaceColor',C(4,:));
scatter(abs(Type_BD.mean_Right_PostPI(LT52_idx)-Type_BD.mean_Right_PrePI(LT52_idx)), Type_BD.weighted_Jaccard(LT52_idx), 'filled', ...
    'MarkerFaceAlpha',1,'MarkerFaceColor',C(5,:));
scatter(abs(Type_BD.mean_Right_PostPI(aMe1_idx)-Type_BD.mean_Right_PrePI(aMe1_idx)), Type_BD.weighted_Jaccard(aMe1_idx), 'filled', ...
    'MarkerFaceAlpha',1,'MarkerFaceColor',C(6,:));
scatter(abs(Type_BD.mean_Right_PostPI(cLM01_idx)-Type_BD.mean_Right_PrePI(cLM01_idx)), Type_BD.weighted_Jaccard(cLM01_idx), 'filled', ...
    'MarkerFaceAlpha',1,'MarkerFaceColor',C(7,:));
hold off;
xlim([0 2])
ylim([0 0.3])

%% =======================================================
%  Z-score 분석 (FF / FB 대조군)
% ========================================================

% --- FF 그룹을 대조군으로 Z-score 계산 ---
ctrl_FF = FF_Reciprocal(:);
ctrl_FF = ctrl_FF(isfinite(ctrl_FF)); % NaN 값 제거

mu_FF = mean(ctrl_FF);
sd_FF = std(ctrl_FF);

% Z-score 계산
z_scores_FF = (x_case - mu_FF) ./ sd_FF;

% 결과 테이블 생성 및 출력
T_Z_FF = table(case_names, x_case, ...
    repmat(mu_FF, numel(x_case), 1), repmat(sd_FF, numel(x_case), 1), ...
    z_scores_FF, ...
    'VariableNames', {'케이스', '케이스_평균값', '대조군_평균', '대조군_표준편차', 'Z_score'});

disp('=== Z-score 분석 (대조군: FF 유형) ===');
disp(T_Z_FF);

% --- FB 그룹을 대조군으로 Z-score 계산 ---
ctrl_FB = FB_Reciprocal(:);
ctrl_FB = ctrl_FB(isfinite(ctrl_FB)); % NaN 값 제거

mu_FB = mean(ctrl_FB);
sd_FB = std(ctrl_FB);

% Z-score 계산
z_scores_FB = (x_case - mu_FB) ./ sd_FB;

% 결과 테이블 생성 및 출력
T_Z_FB = table(case_names, x_case, ...
    repmat(mu_FB, numel(x_case), 1), repmat(sd_FB, numel(x_case), 1), ...
    z_scores_FB, ...
    'VariableNames', {'케이스', '케이스_평균값', '대조군_평균', '대조군_표준편차', 'Z_score'});

disp('=== Z-score 분석 (대조군: FB 유형) ===');
disp(T_Z_FB);
%% ---- Local helper (네가 기존에 ternary를 안만들었을 경우 대비) ----
function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end

function [J, Jw] = jaccard_calculate(InputNeurons, OutputNeurons)
% InputNeurons / OutputNeurons: n×2, [root_id, weight] (weight≥0 가정)

idA = (InputNeurons(:,1));  wA = InputNeurons(:,2);
idB = (OutputNeurons(:,1)); wB = OutputNeurons(:,2);

% --- (1) Binary Jaccard: 존재 유무만 비교 ---
setA = idA(wA > 0);
setB = idB(wB > 0);
uniAB = union(setA, setB);
if isempty(uniAB)
    J = 0;
else
    J = numel(intersect(setA, setB)) / numel(uniAB);
end

% --- (2) Weighted Jaccard: Σmin / Σmax ---
allIDs = union(idA, idB);
[lia, locA] = ismember(allIDs, idA);
[lib, locB] = ismember(allIDs, idB);

wAfull = zeros(size(allIDs));
wBfull = zeros(size(allIDs));
wAfull(lia) = wA(locA(lia));
wBfull(lib) = wB(locB(lib));

num = sum(min(wAfull, wBfull));          % Σ min
den = sum(max(wAfull, wBfull));          % Σ max (= sum(wA)+sum(wB)-num)
if den == 0
    Jw = 0;
else
    Jw = num / den;
end
end

function [J, Jw] = jaccard_calculate_type(InputNeurons, OutputNeurons)
% InputNeurons / OutputNeurons: n×2, [root_id, weight] (weight≥0 가정)

idA = (InputNeurons(:,1));  wA = cell2mat(InputNeurons(:,2));
idB = (OutputNeurons(:,1)); wB = cell2mat(OutputNeurons(:,2));

% --- (1) Binary Jaccard: 존재 유무만 비교 ---
setA = idA(wA > 0);
setB = idB(wB > 0);
uniAB = union(setA, setB);
if isempty(uniAB)
    J = 0;
else
    J = numel(intersect(setA, setB)) / numel(uniAB);
end

% --- (2) Weighted Jaccard: Σmin / Σmax ---
allIDs = union(idA, idB);
[lia, locA] = ismember(allIDs, idA);
[lib, locB] = ismember(allIDs, idB);

wAfull = zeros(size(allIDs));
wBfull = zeros(size(allIDs));
wAfull(lia) = wA(locA(lia));
wBfull(lib) = wB(locB(lib));

num = sum(min(wAfull, wBfull));          % Σ min
den = sum(max(wAfull, wBfull));          % Σ max (= sum(wA)+sum(wB)-num)
if den == 0
    Jw = 0;
else
    Jw = num / den;
end
end

