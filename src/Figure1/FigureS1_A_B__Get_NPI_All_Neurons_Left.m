clear all; close all; clc
load('FAFB_NPI_Thr0.mat')
%%
% FAFBNPIs.In_Synapse_Total=FAFBNPIs.In_Synapse_Optic_R+FAFBNPIs.In_Synapse_Central;
% FAFBNPIs.Out_Synapse_Total=FAFBNPIs.Out_Synapse_Optic_R+FAFBNPIs.Out_Synapse_Central;
FAFBNPIs.In_Synapse_Total=FAFBNPIs.In_Synapse_Optic_R+FAFBNPIs.In_Synapse_Central+FAFBNPIs.In_Synapse_Optic_L;
FAFBNPIs.Out_Synapse_Total=FAFBNPIs.Out_Synapse_Optic_R+FAFBNPIs.Out_Synapse_Central+FAFBNPIs.Out_Synapse_Optic_L;

Right_PostPI=FAFBNPIs.In_Synapse_Optic_R./(FAFBNPIs.In_Synapse_Optic_R+FAFBNPIs.In_Synapse_Central)-FAFBNPIs.In_Synapse_Central./(FAFBNPIs.In_Synapse_Optic_R+FAFBNPIs.In_Synapse_Central);
Right_PrePI=FAFBNPIs.Out_Synapse_Optic_R./(FAFBNPIs.Out_Synapse_Optic_R+FAFBNPIs.Out_Synapse_Central)-FAFBNPIs.Out_Synapse_Central./(FAFBNPIs.Out_Synapse_Optic_R+FAFBNPIs.Out_Synapse_Central);

Left_PostPI=FAFBNPIs.In_Synapse_Optic_L./(FAFBNPIs.In_Synapse_Optic_L+FAFBNPIs.In_Synapse_Central)-FAFBNPIs.In_Synapse_Central./(FAFBNPIs.In_Synapse_Optic_L+FAFBNPIs.In_Synapse_Central);
Left_PrePI=FAFBNPIs.Out_Synapse_Optic_L./(FAFBNPIs.Out_Synapse_Optic_L+FAFBNPIs.Out_Synapse_Central)-FAFBNPIs.Out_Synapse_Central./(FAFBNPIs.Out_Synapse_Optic_L+FAFBNPIs.Out_Synapse_Central);


FAFBNPIs.Right_PostPI=Right_PostPI;
FAFBNPIs.Right_PrePI=Right_PrePI;
FAFBNPIs.Left_PostPI=Left_PostPI;
FAFBNPIs.Left_PrePI=Left_PrePI;

FAFBNPIs_All=FAFBNPIs;

%%%% condition 1
% idx= (FAFBNPIs.In_Synapse_Optic_R>=10 & FAFBNPIs.Out_Synapse_Central>=10)|(FAFBNPIs.Out_Synapse_Optic_R>=10&FAFBNPIs.In_Synapse_Central>=10);
idx= (FAFBNPIs.In_Synapse_Optic_L>=5 & FAFBNPIs.Out_Synapse_Central>=5)|(FAFBNPIs.Out_Synapse_Optic_L>=5&FAFBNPIs.In_Synapse_Central>=5);

FAFBNPIs=FAFBNPIs(idx,:); 

%% 가장 큰 입력이 Optic R 이면 거르기 condition 2

% idx= (FAFBNPIs.In_Synapse_Optic_L>FAFBNPIs.In_Synapse_Optic_R)&(FAFBNPIs.In_Synapse_Optic_L>FAFBNPIs.In_Synapse_Central);
% FAFBNPIs(idx,:)=[];

Thr_R=0.7;
idx= (FAFBNPIs.In_Synapse_Optic_R>(FAFBNPIs.In_Synapse_Optic_R+FAFBNPIs.In_Synapse_Optic_L+FAFBNPIs.In_Synapse_Central)*Thr_R)...
    |(FAFBNPIs.Out_Synapse_Optic_R>(FAFBNPIs.Out_Synapse_Optic_R+FAFBNPIs.Out_Synapse_Optic_L+FAFBNPIs.Out_Synapse_Central)*Thr_R);
    
FAFBNPIs(idx,:)=[];

%% 타입별 계산해보기
[UniqueTypes,~,ic] = unique(FAFBNPIs.type);
FAFBNPIs_by_type=table(UniqueTypes,'VariableNames',{'type'});

for i=1:size(FAFBNPIs_by_type,1)
    idx=strcmp(FAFBNPIs.type,FAFBNPIs_by_type.type{i});
    currentSuperclass=FAFBNPIs.superclass(idx);
    [uniquecurrentSuperclass,~,ic]=unique(currentSuperclass);
    a_counts = accumarray(ic,1);
    [~,maxidx] = max(a_counts);
    FAFBNPIs_by_type.SuperClass{i}=uniquecurrentSuperclass{maxidx};
    FAFBNPIs_by_type.Mean_Left_PostPI(i)=mean(FAFBNPIs.Left_PostPI(idx));
    FAFBNPIs_by_type.Std_Left_PostPI(i)=std(FAFBNPIs.Left_PostPI(idx));
    FAFBNPIs_by_type.Mean_Left_PrePI(i)=mean(FAFBNPIs.Left_PrePI(idx));
    FAFBNPIs_by_type.Std_Left_PrePI(i)=std(FAFBNPIs.Left_PrePI(idx));
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
    scatter(FAFBNPIs.Left_PostPI(idx),FAFBNPIs.Left_PrePI(idx),27,colors(i,:),"filled",'MarkerFaceAlpha',0.65);
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
    scatter((FAFBNPIs_by_type.Mean_Left_PostPI(idx)),(FAFBNPIs_by_type.Mean_Left_PrePI(idx)),27,colors(i,:),"filled",'MarkerFaceAlpha',0.65);
    hold on;

    
end
axis square
grid on
xlabel('Input CenterOfMass')
ylabel('Output CenterOfMass')
set(gca,'TickDir','Out','Box','off','XTick',-1:0.2:1,'YTick',-1:0.2:1)
xlim([-1.05 1.05])
ylim([-1.05 1.05])
print(gcf,'-depsc2','-vector','figureS1_Dots.eps')

%%
LeftFF_by_type=FAFBNPIs_by_type((FAFBNPIs_by_type.Mean_Left_PostPI-FAFBNPIs_by_type.Mean_Left_PrePI)>=0.2,:);
LeftFB_by_type=FAFBNPIs_by_type((FAFBNPIs_by_type.Mean_Left_PostPI-FAFBNPIs_by_type.Mean_Left_PrePI)<=-0.2,:);
LeftBD_by_type=FAFBNPIs_by_type(((FAFBNPIs_by_type.Mean_Left_PostPI-FAFBNPIs_by_type.Mean_Left_PrePI)<0.2)&((FAFBNPIs_by_type.Mean_Left_PostPI-FAFBNPIs_by_type.Mean_Left_PrePI)>-0.2),:);


LeftFF_NPIs=FAFBNPIs(ismember(FAFBNPIs.type,LeftFF_by_type.type),:);
LeftFB_NPIs=FAFBNPIs(ismember(FAFBNPIs.type,LeftFB_by_type.type),:);
LeftBD_NPIs=FAFBNPIs(ismember(FAFBNPIs.type,LeftBD_by_type.type),:);


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
        h=histogram(FAFBNPIs.Left_PostPI(idx)+FAFBNPIs.Left_PrePI(idx),-2:0.4:2,'Normalization','probability');
        h.FaceColor=colors(i+3,:);
        h.EdgeColor=colors(i+3,:);
        h.FaceAlpha = 0.8; % 반투명 효과 (0: 완전 투명, 1: 불투명)
        
        % 히스토그램 계산
        % edges = -2:0.2:2; % 구간 정의
        % [N, edges] = histcounts(FAFBNPIs.Left_PostPI(idx)+FAFBNPIs.Left_PrePI(idx), edges, 'Normalization', 'probability');
        
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
print(gcf,'-depsc2','-vector','figureS1_Plus.eps')

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
        h=histogram(FAFBNPIs.Left_PostPI(idx)-FAFBNPIs.Left_PrePI(idx),-2:0.4:2,'Normalization','probability');
        h.FaceColor=colors(i+3,:);
        h.EdgeColor=colors(i+3,:);
        h.FaceAlpha = 0.8; % 반투명 효과 (0: 완전 투명, 1: 불투명)
        
        % % 히스토그램 계산
        % edges = -2:0.2:2; % 구간 정의
        % [N, edges] = histcounts(FAFBNPIs.Left_PostPI(idx)-FAFBNPIs.Left_PrePI(idx), edges, 'Normalization', 'probability');
        
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
print(gcf,'-depsc2','-vector','figureS1_Minus.eps')

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
        % h=histogram(FAFBNPIs.Left_PostPI(idx),-2:0.1:2,'Normalization','probability');
        % h.FaceColor=colors(i+3,:);
        % h.EdgeColor=colors(i+3,:);
        % h.FaceAlpha = 0.5; % 반투명 효과 (0: 완전 투명, 1: 불투명)
        
        % 히스토그램 계산
        edges = -1:0.2:1; % 구간 정의
        [N, edges] = histcounts(FAFBNPIs.Left_PostPI(idx), edges, 'Normalization', 'probability');
        
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
% ylim([0 1])
print(gcf,'-depsc2','-vector','figureS1_PostPI.eps')

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
        % h=histogram(FAFBNPIs.Left_PostPI(idx),-2:0.2:2,'Normalization','probability');
        % h.FaceColor=colors(i+3,:);
        % h.EdgeColor=colors(i+3,:);
        % h.FaceAlpha = 0.5; % 반투명 효과 (0: 완전 투명, 1: 불투명)
        
        % 히스토그램 계산
        edges = -1:0.2:1; % 구간 정의
        [N, edges] = histcounts(FAFBNPIs.Left_PrePI(idx), edges, 'Normalization', 'probability');
        
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
print(gcf,'-depsc2','-vector','figureS1_PrePI.eps')

%%
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

figure(9);set(gcf,'Color','w')
plot(linspace(-1,0.8,10),linspace(-0.8,1,10),'Color','#EAEBEB','LineWidth',1); hold on;
plot(linspace(-0.8,1,10),linspace(-1,0.8,10),'Color','#EAEBEB','LineWidth',1);

for i=1:1:length(uniqueSuperclass)
    WantToSee=uniqueSuperclass{i}
    idx=strcmp(FAFBNPIs.superclass,WantToSee);
    scatter(FAFBNPIs.Left_PostPI(idx)+FAFBNPIs.Left_PrePI(idx),FAFBNPIs.Left_PostPI(idx)-FAFBNPIs.Left_PrePI(idx),27,colors(i,:),"filled",'MarkerFaceAlpha',0.65);
    hold on;

end
hold off;
axis square
% grid on
xlabel('postPI+prePI')
ylabel('postPI-prePI')
set(gca,'TickDir','Out','Box','off','XTick',-2:0.2:2,'YTick',-2:0.2:2)
xlim([-2.05 2.05])
ylim([-2.05 2.05])
% print(gcf,'-depsc2','-vector','figure1_PM_Dia_v2.eps')


