clear all; close all; clc
load('FAFB_NPI_Thr0.mat')
FAFBNPIs_All=FAFBNPIs;

load('Right_Neurons_All.mat')
load('Left_Neurons_All.mat')

%%
[UniqueTypesRight,~,ic] = unique(Right_Neurons_All.type);
Right_by_type=table(UniqueTypesRight,'VariableNames',{'type'});

for i=1:size(Right_by_type,1)
    idx=strcmp(Right_Neurons_All.type,Right_by_type.type{i});
    currentSuperclass=Right_Neurons_All.superclass(idx);
    [uniquecurrentSuperclass,~,ic]=unique(currentSuperclass);
    a_counts = accumarray(ic,1);
    [~,maxidx] = max(a_counts);
    Right_by_type.SuperClass{i}=uniquecurrentSuperclass{maxidx};
    Right_by_type.Mean_Right_PostPI{i}=mean(Right_Neurons_All.Right_PostPI(idx),'omitmissing');
    Right_by_type.Std_Right_PostPI{i}=std(Right_Neurons_All.Right_PostPI(idx),'omitmissing');
    Right_by_type.Mean_Right_PrePI{i}=mean(Right_Neurons_All.Right_PrePI(idx),'omitmissing');
    Right_by_type.Std_Right_PrePI{i}=std(Right_Neurons_All.Right_PrePI(idx),'omitmissing');
    Right_by_type.Number_of_neurons(i)=sum(strcmp(Right_Neurons_All.type,Right_by_type.type{i}));
    Right_by_type.Total_number_of_neurons(i)=sum(strcmp(FAFBNPIs_All.type,Right_by_type.type{i}));
end
%%
[UniqueTypesLeft,~,ic] = unique(Left_Neurons_All.type);
Left_by_type=table(UniqueTypesLeft,'VariableNames',{'type'});

for i=1:size(Left_by_type,1)
    idx=strcmp(Left_Neurons_All.type,Left_by_type.type{i});
    currentSuperclass=Left_Neurons_All.superclass(idx);
    [uniquecurrentSuperclass,~,ic]=unique(currentSuperclass);
    a_counts = accumarray(ic,1);
    [~,maxidx] = max(a_counts);
    Left_by_type.SuperClass{i}=uniquecurrentSuperclass{maxidx};
    Left_by_type.Mean_Left_PostPI{i}=mean(Left_Neurons_All.Left_PostPI(idx),'omitmissing');
    Left_by_type.Std_Left_PostPI{i}=std(Left_Neurons_All.Left_PostPI(idx),'omitmissing');
    Left_by_type.Mean_Left_PrePI{i}=mean(Left_Neurons_All.Left_PrePI(idx),'omitmissing');
    Left_by_type.Std_Left_PrePI{i}=std(Left_Neurons_All.Left_PrePI(idx),'omitmissing');
    Left_by_type.Number_of_neurons(i)=sum(strcmp(Left_Neurons_All.type,Left_by_type.type{i}));
    Left_by_type.Total_number_of_neurons(i)=sum(strcmp(FAFBNPIs_All.type,Left_by_type.type{i}));
end

%%
Distance_LR=[];

idx=ismember(Left_by_type.type,Right_by_type.type);
uniqueTypes(:,1)=Left_by_type.type(idx);
uniqueTypes=table(uniqueTypes,'VariableNames',{'type'});

uniqueTypes.SuperClass=Left_by_type.SuperClass(idx);
for i=1:1:size(uniqueTypes,1)
    idx_L=strcmp(Left_by_type.type,uniqueTypes.type{i});
    idx_R=strcmp(Right_by_type.type,uniqueTypes.type{i});
    uniqueTypes.Distance(i)=sqrt((Left_by_type.Mean_Left_PostPI{idx_L}-Right_by_type.Mean_Right_PostPI{idx_R})^2+...
        (Left_by_type.Mean_Left_PrePI{idx_L}-Right_by_type.Mean_Right_PrePI{idx_R})^2);
end
%%
idx_VPN=strcmp(uniqueTypes.SuperClass,'visual_projection');
VPN_Distance=uniqueTypes.Distance(idx_VPN);
VPN_Table=uniqueTypes(idx_VPN,:);

idx_VCN=strcmp(uniqueTypes.SuperClass,'visual_centrifugal');
VCN_Distance=uniqueTypes.Distance(idx_VCN);
VCN_Table=uniqueTypes(idx_VCN,:);

idx_Optic=strcmp(uniqueTypes.SuperClass,'optic');
Optic_Distance=uniqueTypes.Distance(idx_Optic);
Optic_Table=uniqueTypes(idx_Optic,:);

idx_Central=strcmp(uniqueTypes.SuperClass,'central');
Central_Distance=uniqueTypes.Distance(idx_Central);
Central_Table=uniqueTypes(idx_Central,:);

%%
% 색상 정의 (각 그룹에 대응)
colors = [ 0.0000, 0.4470, 0.7410;  % VPN (visual projection)
           0.4660, 0.6740, 0.1880;  % VCN (visual centrifugal)
           0.9290, 0.6940, 0.1250;  % Optic
           0.4940, 0.1840, 0.5560]; % Central
% 점 색상을 어둡게 하기 위한 스케일
darken_ratio = 0.8;  % 0~1 사이 값 (작을수록 더 어두움)
scatter_colors = colors * darken_ratio;

figure(1); clf; set(gcf,'Color','w'); hold on;

% 평균과 표준편차 계산
data = [mean(VPN_Table.Distance), mean(VCN_Table.Distance), mean(Optic_Table.Distance), mean(Central_Table.Distance)];
err  = [std(VPN_Table.Distance),  std(VCN_Table.Distance),  std(Optic_Table.Distance),  std(Central_Table.Distance)];

x = ["VPN", "VCN", "Optic", "Central"];
bar_handle = bar(x, data, 'FaceColor', 'flat','EdgeColor','flat');  % FaceColor를 'flat'으로 설정해야 CData가 적용됨

% 색상 적용
for i = 1:numel(data)
    bar_handle.CData(i, :) = colors(i, :);
end

% x 위치 추출
xtick_positions = bar_handle.XEndPoints;

% 에러바 추가
errorbar(xtick_positions-0.05, data, err, 'k.', 'LineWidth', 1.5);

% 개별 데이터 점 찍기 (scatter)
scatter(ones(size(VPN_Table.Distance))   * xtick_positions(1)+0.05, VPN_Table.Distance,   10, scatter_colors(1,:), 'filled', 'MarkerFaceAlpha', 0.9);
scatter(ones(size(VCN_Table.Distance))   * xtick_positions(2)+0.05, VCN_Table.Distance,   10, scatter_colors(2,:), 'filled', 'MarkerFaceAlpha', 0.9);
scatter(ones(size(Optic_Table.Distance)) * xtick_positions(3)+0.05, Optic_Table.Distance, 10, scatter_colors(3,:), 'filled', 'MarkerFaceAlpha', 0.9);
scatter(ones(size(Central_Table.Distance))* xtick_positions(4)+0.05, Central_Table.Distance, 10, scatter_colors(4,:), 'filled', 'MarkerFaceAlpha', 0.9);

% 시각적 설정
set(gca, 'Box', 'off', 'TickDir', 'out');
ylim([0 2]);
ylabel('Distance');
title('Mean Distance with Standard Deviation');
grid on
print(gcf,'-depsc2','-vector','figureS1_LR_comparison.eps')

% %%
% % 두 점 집합: 각각 [Pre PI, Post PI] 좌표들
% X = [Right_Neurons_All.Right_PrePI, Right_Neurons_All.Right_PostPI];
% Y = [Left_Neurons_All.Left_PrePI, Left_Neurons_All.Left_PostPI];
% 
% % Procrustes 분석
% [d, Z, transform] = procrustes(X, Y);
% 
% % d: 유사도 지표 (0~1)
% % Z: Y를 변형시켜 X에 맞춘 점들
% % transform: 회전, 스케일, 이동 정보 포함
