clear all; close all; clc
load('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\Figure1_Overall\FAFB_NPI_Thr0.mat')
%%
% FAFBNPIs.In_Synapse_Total=FAFBNPIs.In_Synapse_Optic_R+FAFBNPIs.In_Synapse_Central;
% FAFBNPIs.Out_Synapse_Total=FAFBNPIs.Out_Synapse_Optic_R+FAFBNPIs.Out_Synapse_Central;
FAFBNPIs.In_Synapse_Total=FAFBNPIs.In_Synapse_Optic_R+FAFBNPIs.In_Synapse_Central+FAFBNPIs.In_Synapse_Optic_L;
FAFBNPIs.Out_Synapse_Total=FAFBNPIs.Out_Synapse_Optic_R+FAFBNPIs.Out_Synapse_Central+FAFBNPIs.Out_Synapse_Optic_L;

Bilateral_PostPI=FAFBNPIs.In_Synapse_Optic_R./(FAFBNPIs.In_Synapse_Optic_R+FAFBNPIs.In_Synapse_Optic_L)-FAFBNPIs.In_Synapse_Optic_L./(FAFBNPIs.In_Synapse_Optic_R+FAFBNPIs.In_Synapse_Optic_L);
Bilateral_PrePI=FAFBNPIs.Out_Synapse_Optic_R./(FAFBNPIs.Out_Synapse_Optic_R+FAFBNPIs.Out_Synapse_Optic_L)-FAFBNPIs.Out_Synapse_Optic_L./(FAFBNPIs.Out_Synapse_Optic_R+FAFBNPIs.Out_Synapse_Optic_L);

FAFBNPIs.Bilateral_PostPI=Bilateral_PostPI;
FAFBNPIs.Bilateral_PrePI=Bilateral_PrePI;
FAFBNPIs_All=FAFBNPIs;

%%%% condition 1
% idx= (FAFBNPIs.In_Synapse_Optic_R>=10 & FAFBNPIs.Out_Synapse_Central>=10)|(FAFBNPIs.Out_Synapse_Optic_R>=10&FAFBNPIs.In_Synapse_Central>=10);
idx= (FAFBNPIs.In_Synapse_Optic_R>=5 & FAFBNPIs.Out_Synapse_Optic_L>=5)|(FAFBNPIs.Out_Synapse_Optic_R>=5&FAFBNPIs.In_Synapse_Optic_L>=5);

FAFBNPIs=FAFBNPIs(idx,:); 

%% 가장 큰 입력이 Optic L 이면 거르기 condition 2

% idx= (FAFBNPIs.In_Synapse_Optic_L>FAFBNPIs.In_Synapse_Optic_R)&(FAFBNPIs.In_Synapse_Optic_L>FAFBNPIs.In_Synapse_Central);
% FAFBNPIs(idx,:)=[];

Thr_L=0.7;
idx= (FAFBNPIs.In_Synapse_Central>(FAFBNPIs.In_Synapse_Optic_R+FAFBNPIs.In_Synapse_Optic_L+FAFBNPIs.In_Synapse_Central)*Thr_L)...
    |(FAFBNPIs.Out_Synapse_Central>(FAFBNPIs.Out_Synapse_Optic_R+FAFBNPIs.Out_Synapse_Optic_L+FAFBNPIs.Out_Synapse_Central)*Thr_L);
temp=FAFBNPIs(idx,:);
FAFBNPIs(idx,:)=[];

for i=1:1:size(FAFBNPIs,1)
    FAFBNPIs.type{i}=[FAFBNPIs.type{i} '_' FAFBNPIs.side{i}];
end
%% 타입별 계산해보기
[UniqueTypes,~,ic] = unique(FAFBNPIs.type);
FAFBNPIs_by_type=table(UniqueTypes,'VariableNames',{'type'});

for i=1:size(FAFBNPIs_by_type,1)
    idx=strcmp(FAFBNPIs.type,FAFBNPIs_by_type.type{i});
    currentSuperclass=FAFBNPIs.superclass(idx);
    [uniquecurrentSuperclass,~,ic]=unique(currentSuperclass);
    a_counts = accumarray(ic,1);
    [~,maxidx] = max(a_counts);
    FAFBNPIs_by_type.SuperClass{i}=currentSuperclass{maxidx};
    FAFBNPIs_by_type.Mean_Bilateral_PostPI(i)=mean(FAFBNPIs.Bilateral_PostPI(idx));
    FAFBNPIs_by_type.Std_Bilateral_PostPI(i)=std(FAFBNPIs.Bilateral_PostPI(idx));
    FAFBNPIs_by_type.Mean_Bilateral_PrePI(i)=mean(FAFBNPIs.Bilateral_PrePI(idx));
    FAFBNPIs_by_type.Std_Bilateral_PrePI(i)=std(FAFBNPIs.Bilateral_PrePI(idx));
    FAFBNPIs_by_type.Number_of_neurons(i)=sum(strcmp(FAFBNPIs.type,FAFBNPIs_by_type.type{i}));
    FAFBNPIs_by_type.Total_number_of_neurons(i)=sum(strcmp(FAFBNPIs_All.type,FAFBNPIs_by_type.type{i}));
end
%% 10퍼센트 미만 거르기 condition 3
idx=(FAFBNPIs_by_type.Number_of_neurons./FAFBNPIs_by_type.Total_number_of_neurons*100)<20;
FAFBNPIs_by_type(idx,:)=[];
FAFBNPIs(~ismember(FAFBNPIs.type,FAFBNPIs_by_type.type),:)=[];

%% superclass 별 그림..
%%
% 'visual_projection'  0.0000, 0.4470, 0.7410;  % 파랑
% 'visual_centrifugal' 0.4660, 0.6740, 0.1880;  % 초록
% 'sensory' 0.6350, 0.5090, 0.2540;  % 갈색
% 'optic' 노랑 0.9290, 0.6940, 0.1250;  % 노랑
% % 'motor' 0.2780, 0.6000, 0.8000;  % 하늘색
% % 'endocrine'  0.3010, 0.7450, 0.9330;  % 청록
% 'descending' 0.8500, 0.3250, 0.0980;  % 주황
% 'central'  0.4940, 0.1840, 0.5560;  % 보라
% 'ascending' 0.6350, 0.0780, 0.1840;  % 빨강
uniqueSuperclass={'endocrine',...
    'motor',...
    'sensory',...
    'visual_projection',...
    'visual_centrifugal',...
    'optic',...
    'central',...
    'descending',...
    'ascending'};
% uniqueSuperclass={'descending', 'ascending'};
colors = [ 0.3010, 0.7450, 0.9330; %endocrine
    0.2780, 0.6000, 0.8000; %motor
    0.8500, 0.3250, 0.0980; %sensory
    0.0000, 0.4470, 0.7410; %visual_projection
    0.4660, 0.6740, 0.1880; %visual_centrifugal
    0.9290, 0.6940, 0.1250; % 'optic',...
    0.4940, 0.1840, 0.5560;  % 'central',...
    0.8500, 0.1500, 0.2000; % 'descending',...
    0.6350, 0.5090, 0.2540 % 'ascending'
    ];

figure(1);set(gcf,'Color','w'); hold on;
plot(linspace(-1,0.8,10),linspace(-0.8,1,10),'Color','#EAEBEB','LineWidth',1);
plot(linspace(-0.8,1,10),linspace(-1,0.8,10),'Color','#EAEBEB','LineWidth',1);

for i=1:1:length(uniqueSuperclass)
    WantToSee=uniqueSuperclass{i};
    idx=strcmp(FAFBNPIs.superclass,WantToSee);
    scatter(FAFBNPIs.Bilateral_PostPI(idx),FAFBNPIs.Bilateral_PrePI(idx),27,colors(i,:),"filled",'MarkerFaceAlpha',0.65);
    hold on;

end
hold off;
axis square
% grid on
xlabel('Input CenterOfMass')
ylabel('Output CenterOfMass')
set(gca,'TickDir','Out','Box','off','XTick',-1:0.2:1,'YTick',-1:0.2:1)
xlim([-1.05 1.05])
ylim([-1.05 1.05])

figure(2);set(gcf,'Color','w')
plot(linspace(-1,0.8,10),linspace(-0.8,1,10),'Color','#EAEBEB','LineWidth',1); hold on;
plot(linspace(-0.8,1,10),linspace(-1,0.8,10),'Color','#EAEBEB','LineWidth',1);

uniqueSuperclass={'endocrine',...
    'motor',...
    'sensory',...
    'visual_projection',...
    'visual_centrifugal',...
    'optic',...
    'central',...
    'descending',...
    'ascending'};
for i=1:1:length(uniqueSuperclass)
    WantToSee=uniqueSuperclass{i};
    idx=strcmp(FAFBNPIs_by_type.SuperClass,WantToSee);
    scatter((FAFBNPIs_by_type.Mean_Bilateral_PostPI(idx)),(FAFBNPIs_by_type.Mean_Bilateral_PrePI(idx)),27,colors(i,:),"filled",'MarkerFaceAlpha',0.65);
    hold on;

    
end
axis square
grid on
xlabel('Input CenterOfMass')
ylabel('Output CenterOfMass')
set(gca,'TickDir','Out','Box','off','XTick',-1:0.2:1,'YTick',-1:0.2:1)
xlim([-1.05 1.05])
ylim([-1.05 1.05])
print(gcf,'-depsc2','-vector','figureBL_Dots_Type.eps')

%%
Bilateral_R_by_type=FAFBNPIs_by_type((FAFBNPIs_by_type.Mean_Bilateral_PostPI-FAFBNPIs_by_type.Mean_Bilateral_PrePI)>=0.2,:);
Bilateral_L_by_type=FAFBNPIs_by_type((FAFBNPIs_by_type.Mean_Bilateral_PostPI-FAFBNPIs_by_type.Mean_Bilateral_PrePI)<=-0.2,:);
Bilateral_BD_type=FAFBNPIs_by_type(((FAFBNPIs_by_type.Mean_Bilateral_PostPI-FAFBNPIs_by_type.Mean_Bilateral_PrePI)<0.2)&((FAFBNPIs_by_type.Mean_Bilateral_PostPI-FAFBNPIs_by_type.Mean_Bilateral_PrePI)>-0.2),:);


Bilateral_R_NPIs=FAFBNPIs(ismember(FAFBNPIs.type,Bilateral_R_by_type.type),:);
Bilateral_L_NPIs=FAFBNPIs(ismember(FAFBNPIs.type,Bilateral_L_by_type.type),:);
Bilateral_BD_NPIs=FAFBNPIs(ismember(FAFBNPIs.type,Bilateral_BD_type.type),:);


%%
% 'visual_projection'  0.0000, 0.4470, 0.7410;  % 파랑
% 'visual_centrifugal' 0.4660, 0.6740, 0.1880;  % 초록
% 'sensory' 0.6350, 0.5090, 0.2540;  % 갈색
% 'optic' 노랑 0.9290, 0.6940, 0.1250;  % 노랑
% % 'motor' 0.2780, 0.6000, 0.8000;  % 하늘색
% % 'endocrine'  0.3010, 0.7450, 0.9330;  % 청록
% 'descending' 0.8500, 0.3250, 0.0980;  % 주황
% 'central'  0.4940, 0.1840, 0.5560;  % 보라
% 'ascending' 0.6350, 0.0780, 0.1840;  % 빨강
% uniqueSuperclass={'endocrine',...
%     'motor',...
%     'sensory',...
%     'visual_projection',...
%     'visual_centrifugal',...
%     'optic',...
%     'central',...
%     'descending',...
%     'ascending'};
uniqueSuperclass={
    'visual_projection',...
    'visual_centrifugal',...
    'optic',...
    'central'};
colors = [ 0.3010, 0.7450, 0.9330; %endocrine
    0.2780, 0.6000, 0.8000; %motor
    0.8500, 0.3250, 0.0980; %sensory
    0.0000, 0.4470, 0.7410; %visual_projection
    0.4660, 0.6740, 0.1880; %visual_centrifugal
    0.9290, 0.6940, 0.1250; % 'optic',...
    0.4940, 0.1840, 0.5560;  % 'central',...
    0.8500, 0.1500, 0.2000; % 'descending',...
    0.6350, 0.5090, 0.2540 % 'ascending'
    ];

figure(3);set(gcf,'Color','w')

for i=1:1:length(uniqueSuperclass)
    WantToSee=uniqueSuperclass{i}
    idx=strcmp(FAFBNPIs.superclass,WantToSee);
    hold on;
    if any(idx)
        h=histogram(FAFBNPIs.Bilateral_PostPI(idx)+FAFBNPIs.Bilateral_PrePI(idx),-2:0.4:2,'Normalization','probability');
        h.FaceColor=colors(i+3,:);
        h.EdgeColor=colors(i+3,:);
        h.FaceAlpha = 0.8; % 반투명 효과 (0: 완전 투명, 1: 불투명)
        
        % 히스토그램 계산
        % edges = -2:0.2:2; % 구간 정의
        % [N, edges] = histcounts(FAFBNPIs.Bilateral_PostPI(idx)+FAFBNPIs.Bilateral_PrePI(idx), edges, 'Normalization', 'probability');
        
        % 구간 중심 계산
        % binCenters = edges(1:end-1) + diff(edges)/2;
        % plot(binCenters, N,'Color', colors(i+3,:), 'LineWidth', 2);

    end

end
hold off;
axis square
xlabel('PostPI+PrePI')
set(gca,'TickDir','Out','Box','off','XTick',-2:0.4:2,'YTick',0:0.2:1)

xlim([-2.05 2.05])
ylim([0 1])
print(gcf,'-depsc2','-vector','figureBL_Plus.eps')

%%
uniqueSuperclass={
    'visual_projection',...
    'visual_centrifugal',...
    'optic',...
    'central'};
colors = [ 0.3010, 0.7450, 0.9330; %endocrine
    0.2780, 0.6000, 0.8000; %motor
    0.8500, 0.3250, 0.0980; %sensory
    0.0000, 0.4470, 0.7410; %visual_projection
    0.4660, 0.6740, 0.1880; %visual_centrifugal
    0.9290, 0.6940, 0.1250; % 'optic',...
    0.4940, 0.1840, 0.5560;  % 'central',...
    0.8500, 0.1500, 0.2000; % 'descending',...
    0.6350, 0.5090, 0.2540 % 'ascending'
    ];

figure(4);set(gcf,'Color','w')

for i=1:1:length(uniqueSuperclass)
    WantToSee=uniqueSuperclass{i}
    idx=strcmp(FAFBNPIs.superclass,WantToSee);
    hold on;
    if any(idx)
        h=histogram(FAFBNPIs.Bilateral_PostPI(idx)-FAFBNPIs.Bilateral_PrePI(idx),-2:0.4:2,'Normalization','probability');
        h.FaceColor=colors(i+3,:);
        h.EdgeColor=colors(i+3,:);
        h.FaceAlpha = 0.8; % 반투명 효과 (0: 완전 투명, 1: 불투명)
        
        % % 히스토그램 계산
        % edges = -2:0.2:2; % 구간 정의
        % [N, edges] = histcounts(FAFBNPIs.Bilateral_PostPI(idx)-FAFBNPIs.Bilateral_PrePI(idx), edges, 'Normalization', 'probability');
        
        % 구간 중심 계산
        % binCenters = edges(1:end-1) + diff(edges)/2;
        % plot(binCenters, N,'Color', colors(i+3,:), 'LineWidth', 2);

    end

end
hold off;
axis square
xlabel('PostPI-PrePI')
set(gca,'TickDir','Out','Box','off','XTick',-2:0.4:2,'YTick',0:0.2:1)
xlim([-2.05 2.05])
ylim([0 1])
print(gcf,'-depsc2','-vector','figureBL_Minus.eps')

%%
uniqueSuperclass={
    'visual_projection',...
    'visual_centrifugal',...
    'optic',...
    'central'};
colors = [ 0.3010, 0.7450, 0.9330; %endocrine
    0.2780, 0.6000, 0.8000; %motor
    0.8500, 0.3250, 0.0980; %sensory
    0.0000, 0.4470, 0.7410; %visual_projection
    0.4660, 0.6740, 0.1880; %visual_centrifugal
    0.9290, 0.6940, 0.1250; % 'optic',...
    0.4940, 0.1840, 0.5560;  % 'central',...
    0.8500, 0.1500, 0.2000; % 'descending',...
    0.6350, 0.5090, 0.2540 % 'ascending'
    ];

figure(5);set(gcf,'Color','w')

for i=1:1:length(uniqueSuperclass)
    WantToSee=uniqueSuperclass{i}
    idx=strcmp(FAFBNPIs.superclass,WantToSee);
    hold on;
    if any(idx)
        % h=histogram(FAFBNPIs.Bilateral_PostPI(idx),-2:0.1:2,'Normalization','probability');
        % h.FaceColor=colors(i+3,:);
        % h.EdgeColor=colors(i+3,:);
        % h.FaceAlpha = 0.5; % 반투명 효과 (0: 완전 투명, 1: 불투명)
        
        % 히스토그램 계산
        edges = -1:0.2:1; % 구간 정의
        [N, edges] = histcounts(FAFBNPIs.Bilateral_PostPI(idx), edges, 'Normalization', 'probability');
        
        % 구간 중심 계산
        binCenters = edges(1:end-1) + diff(edges)/2;
        plot(binCenters, N,'Color', colors(i+3,:), 'LineWidth', 2);

    end

end
hold off;
% axis square
xlabel('PostPI')
set(gca,'TickDir','Out','Box','off','XTick',-2:0.2:2,'YTick',0:0.2:1)

% xlim([-1 1])
ylim([0 1])
print(gcf,'-depsc2','-vector','figureBL_PostPI.eps')

%%
uniqueSuperclass={
    'visual_projection',...
    'visual_centrifugal',...
    'optic',...
    'central'};
colors = [ 0.3010, 0.7450, 0.9330; %endocrine
    0.2780, 0.6000, 0.8000; %motor
    0.8500, 0.3250, 0.0980; %sensory
    0.0000, 0.4470, 0.7410; %visual_projection
    0.4660, 0.6740, 0.1880; %visual_centrifugal
    0.9290, 0.6940, 0.1250; % 'optic',...
    0.4940, 0.1840, 0.5560;  % 'central',...
    0.8500, 0.1500, 0.2000; % 'descending',...
    0.6350, 0.5090, 0.2540 % 'ascending'
    ];

figure(6);set(gcf,'Color','w')

for i=1:1:length(uniqueSuperclass)
    WantToSee=uniqueSuperclass{i}
    idx=strcmp(FAFBNPIs.superclass,WantToSee);
    hold on;
    if any(idx)
        % h=histogram(FAFBNPIs.Bilateral_PostPI(idx),-2:0.2:2,'Normalization','probability');
        % h.FaceColor=colors(i+3,:);
        % h.EdgeColor=colors(i+3,:);
        % h.FaceAlpha = 0.5; % 반투명 효과 (0: 완전 투명, 1: 불투명)
        
        % 히스토그램 계산
        edges = -1:0.2:1; % 구간 정의
        [N, edges] = histcounts(FAFBNPIs.Bilateral_PrePI(idx), edges, 'Normalization', 'probability');
        
        % 구간 중심 계산
        binCenters = edges(1:end-1) + diff(edges)/2;
        plot(binCenters, N,'Color', colors(i+3,:), 'LineWidth', 2);

    end

end
hold off;
% axis square
xlabel('PrePI')
set(gca,'TickDir','Out','Box','off','XTick',-2:0.2:2,'YTick',0:0.2:1)

% xlim([-1 1])
ylim([0 1])
print(gcf,'-depsc2','-vector','figureBL_PrePI.eps')

%%
FAFBNPIs_by_type_LR=FAFBNPIs_by_type;
FAFBNPIs_by_type_LR(57,:) = [];

uniqueSuperclass={
    'visual_projection',...
    'visual_centrifugal',...
    'optic',...
    'ascending'};
Distance_LR=[];
for i=1:2:size(FAFBNPIs_by_type_LR,1)
    left_pos=[FAFBNPIs_by_type_LR.Mean_Bilateral_PostPI(i) FAFBNPIs_by_type_LR.Mean_Bilateral_PrePI(i)];
    right_pos=[FAFBNPIs_by_type_LR.Mean_Bilateral_PostPI(i+1) FAFBNPIs_by_type_LR.Mean_Bilateral_PrePI(i+1)];
    Distance_LR=[Distance_LR norm(left_pos+right_pos)];

end
%%
FAFBNPIs_by_type_LR_SuperClass=FAFBNPIs_by_type_LR.SuperClass(1:2:end);
%%
idx_VPN=strcmp(FAFBNPIs_by_type_LR_SuperClass,'visual_projection');
VPN_Distance=Distance_LR(idx_VPN);

idx_VCN=strcmp(FAFBNPIs_by_type_LR_SuperClass,'visual_centrifugal');
VCN_Distance=Distance_LR(idx_VCN);

idx_Optic=strcmp(FAFBNPIs_by_type_LR_SuperClass,'optic');
Optic_Distance=Distance_LR(idx_Optic);

idx_Ascending=strcmp(FAFBNPIs_by_type_LR_SuperClass,'ascending');
Ascending_Distance=Distance_LR(idx_Ascending);

%%
% 색상 정의 (각 그룹에 대응)
colors = [
    0.0000, 0.4470, 0.7410; %visual_projection
    0.4660, 0.6740, 0.1880; %visual_centrifugal
    0.9290, 0.6940, 0.1250; % 'optic',...
    ];
% 점 색상을 어둡게 하기 위한 스케일
darken_ratio = 0.8;  % 0~1 사이 값 (작을수록 더 어두움)
scatter_colors = colors * darken_ratio;

figure(1); clf; set(gcf,'Color','w'); hold on;

% 평균과 표준편차 계산
data = [mean(VPN_Distance), mean(VCN_Distance), mean(Optic_Distance)];
err  = [std(VPN_Distance),  std(VCN_Distance),  std(Optic_Distance)];

x = ["VPN", "VCN", "Optic"];
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
scatter(ones(size(VPN_Distance))   * xtick_positions(1)+0.05, VPN_Distance,   10, scatter_colors(1,:), 'filled', 'MarkerFaceAlpha', 0.9);
scatter(ones(size(VCN_Distance))   * xtick_positions(2)+0.05, VCN_Distance,   10, scatter_colors(2,:), 'filled', 'MarkerFaceAlpha', 0.9);
scatter(ones(size(Optic_Distance)) * xtick_positions(3)+0.05, Optic_Distance, 10, scatter_colors(3,:), 'filled', 'MarkerFaceAlpha', 0.9);
% scatter(ones(size(Ascending_Distance))* xtick_positions(4)+0.05, Ascending_Distance, 10, scatter_colors(4,:), 'filled', 'MarkerFaceAlpha', 0.9);
% 시각적 설정
set(gca, 'Box', 'off', 'TickDir', 'out');
ylim([0 2]);
ylabel('Distance');
title('Mean Distance with Standard Deviation');
grid off
%%
print(gcf,'-depsc2','-vector','figureBL_LRcomparison.eps')
