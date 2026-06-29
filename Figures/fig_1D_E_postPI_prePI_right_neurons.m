clear all; close all; clc

baseDir = fileparts(fileparts(mfilename('fullpath')));  % resolved from this script's location (portable; no absolute paths)
load(fullfile(baseDir, 'Processed_Data', 'FAFB_NPI_thr0.mat'))
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

%% condition 1
idx= (FAFBNPIs.In_Synapse_Optic_R>=5 & FAFBNPIs.Out_Synapse_Central>=5)|(FAFBNPIs.Out_Synapse_Optic_R>=5&FAFBNPIs.In_Synapse_Central>=5);

FAFBNPIs=FAFBNPIs(idx,:);

%% condition 2: drop neurons dominated by the left optic lobe (Optic L)
Thr_L=0.7;
idx= (FAFBNPIs.In_Synapse_Optic_L>(FAFBNPIs.In_Synapse_Optic_R+FAFBNPIs.In_Synapse_Optic_L+FAFBNPIs.In_Synapse_Central)*Thr_L)...
    |(FAFBNPIs.Out_Synapse_Optic_L>(FAFBNPIs.Out_Synapse_Optic_R+FAFBNPIs.Out_Synapse_Optic_L+FAFBNPIs.Out_Synapse_Central)*Thr_L);
temp=FAFBNPIs(idx,:);
FAFBNPIs(idx,:)=[];

%% Compute per-type statistics
[UniqueTypes,~,ic] = unique(FAFBNPIs.type);
FAFBNPIs_by_type=table(UniqueTypes,'VariableNames',{'type'});

for i=1:size(FAFBNPIs_by_type,1)
    idx=strcmp(FAFBNPIs.type,FAFBNPIs_by_type.type{i});
    currentSuperclass=FAFBNPIs.superclass(idx);
    [uniquecurrentSuperclass,~,ic]=unique(currentSuperclass);
    a_counts = accumarray(ic,1);
    [~,maxidx] = max(a_counts);
    FAFBNPIs_by_type.SuperClass{i}=uniquecurrentSuperclass{maxidx};
    FAFBNPIs_by_type.Mean_Right_PostPI(i)=mean(FAFBNPIs.Right_PostPI(idx));
    FAFBNPIs_by_type.Std_Right_PostPI(i)=std(FAFBNPIs.Right_PostPI(idx));
    FAFBNPIs_by_type.Mean_Right_PrePI(i)=mean(FAFBNPIs.Right_PrePI(idx));
    FAFBNPIs_by_type.Std_Right_PrePI(i)=std(FAFBNPIs.Right_PrePI(idx));
    FAFBNPIs_by_type.Number_of_neurons(i)=sum(strcmp(FAFBNPIs.type,FAFBNPIs_by_type.type{i}));
    FAFBNPIs_by_type.Total_number_of_neurons(i)=sum(strcmp(FAFBNPIs_All.type,FAFBNPIs_by_type.type{i}));
end
%%%%%%%%%%%%%%%%%%%%%%%
%% condition 3: drop types where less than 20% of the type's neurons remain
idx=(FAFBNPIs_by_type.Number_of_neurons./FAFBNPIs_by_type.Total_number_of_neurons*100)<20;
FAFBNPIs_by_type(idx,:)=[];
FAFBNPIs(~ismember(FAFBNPIs.type,FAFBNPIs_by_type.type),:)=[];

%% Plots by superclass
%%
% 'visual_projection'  0.0000, 0.4470, 0.7410;  % blue
% 'visual_centrifugal' 0.4660, 0.6740, 0.1880;  % green
% 'sensory' 0.6350, 0.5090, 0.2540;  % brown
% 'optic' 0.9290, 0.6940, 0.1250;  % yellow
% % 'motor' 0.2780, 0.6000, 0.8000;  % sky blue
% % 'endocrine'  0.3010, 0.7450, 0.9330;  % cyan
% 'descending' 0.8500, 0.3250, 0.0980;  % orange
% 'central'  0.4940, 0.1840, 0.5560;  % purple
% 'ascending' 0.6350, 0.0780, 0.1840;  % red
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

figure(1);set(gcf,'Color','w')
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
    scatter((FAFBNPIs_by_type.Mean_Right_PostPI(idx)),(FAFBNPIs_by_type.Mean_Right_PrePI(idx)),27,colors(i,:),"filled",'MarkerFaceAlpha',0.65);
    hold on;
end
axis square
grid on
xlabel('postPI')
ylabel('prePI')
set(gca,'TickDir','Out','Box','off','XTick',-1:0.2:1,'YTick',-1:0.2:1)
xlim([-1.05 1.05])
ylim([-1.05 1.05])


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

figure(2);set(gcf,'Color','w')

for i=1:1:length(uniqueSuperclass)
    WantToSee=uniqueSuperclass{i}
    idx=strcmp(FAFBNPIs.superclass,WantToSee);
    hold on;
    if any(idx)
        h=histogram(FAFBNPIs.Right_PostPI(idx)+FAFBNPIs.Right_PrePI(idx),-2:0.4:2,'Normalization','probability');
        h.FaceColor=colors(i+3,:);
        h.EdgeColor=colors(i+3,:);
        h.FaceAlpha = 0.8; % transparency (0: fully transparent, 1: opaque)
    end

end
hold off;
axis square
xlabel('PostPI+PrePI')
set(gca,'TickDir','Out','Box','off','XTick',-2:0.4:2,'YTick',0:0.2:1)

xlim([-2.05 2.05])
ylim([0 1])

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

figure(3);set(gcf,'Color','w')

for i=1:1:length(uniqueSuperclass)
    WantToSee=uniqueSuperclass{i}
    idx=strcmp(FAFBNPIs.superclass,WantToSee);
    hold on;
    if any(idx)
        h=histogram(FAFBNPIs.Right_PostPI(idx)-FAFBNPIs.Right_PrePI(idx),-2:0.4:2,'Normalization','probability');
        h.FaceColor=colors(i+3,:);
        h.EdgeColor=colors(i+3,:);
        h.FaceAlpha = 0.8; % transparency (0: fully transparent, 1: opaque)
    end

end
hold off;
axis square
xlabel('PostPI-PrePI')
set(gca,'TickDir','Out','Box','off','XTick',-2:0.4:2,'YTick',0:0.2:1)
xlim([-2.05 2.05])
ylim([0 1])

%% Classify right-hemisphere types (FFP / FBP / BDP) by postPI - prePI
RightFFP_by_type=FAFBNPIs_by_type((FAFBNPIs_by_type.Mean_Right_PostPI-FAFBNPIs_by_type.Mean_Right_PrePI)>=0.2,:);
RightFBP_by_type=FAFBNPIs_by_type((FAFBNPIs_by_type.Mean_Right_PostPI-FAFBNPIs_by_type.Mean_Right_PrePI)<=-0.2,:);
RightBDP_by_type=FAFBNPIs_by_type(((FAFBNPIs_by_type.Mean_Right_PostPI-FAFBNPIs_by_type.Mean_Right_PrePI)<0.2)&((FAFBNPIs_by_type.Mean_Right_PostPI-FAFBNPIs_by_type.Mean_Right_PrePI)>-0.2),:);

RightFFP_NPIs=FAFBNPIs(ismember(FAFBNPIs.type,RightFFP_by_type.type),:);
RightFBP_NPIs=FAFBNPIs(ismember(FAFBNPIs.type,RightFBP_by_type.type),:);
RightBDP_NPIs=FAFBNPIs(ismember(FAFBNPIs.type,RightBDP_by_type.type),:);

%% Split the bidirectional (BDP) types into optic / central / real by postPI + prePI
RightBDP_optic_type=RightBDP_by_type((RightBDP_by_type.Mean_Right_PostPI+RightBDP_by_type.Mean_Right_PrePI)>=1.2,:);
RightBDP_central_type=RightBDP_by_type((RightBDP_by_type.Mean_Right_PostPI+RightBDP_by_type.Mean_Right_PrePI)<=-1.2,:);
RightBDP_real_type=RightBDP_by_type(((RightBDP_by_type.Mean_Right_PostPI+RightBDP_by_type.Mean_Right_PrePI)<1.2&(RightBDP_by_type.Mean_Right_PostPI+RightBDP_by_type.Mean_Right_PrePI)>-1.2),:);

RightBDP_optic_NPIs=RightBDP_NPIs(ismember(RightBDP_NPIs.type,RightBDP_optic_type.type),:);
RightBDP_central_NPIs=RightBDP_NPIs(ismember(RightBDP_NPIs.type,RightBDP_central_type.type),:);
RightBDP_real_NPIs=RightBDP_NPIs(ismember(RightBDP_NPIs.type,RightBDP_real_type.type),:);

%% Save right-hemisphere neuron classification
save(fullfile(baseDir, 'Processed_Data', 'right_neurons_thr0.mat'), ...
    'RightFFP_by_type','RightFFP_NPIs', ...
    'RightFBP_by_type','RightFBP_NPIs', ...
    'RightBDP_optic_type','RightBDP_optic_NPIs', ...
    'RightBDP_central_type','RightBDP_central_NPIs', ...
    'RightBDP_real_type','RightBDP_real_NPIs');

%% Save the root_ids of FFP / FBP / BDP / Others neurons as CSV (no header, root_id only)
% Others = other bidirectional neurons (optic + central BD, i.e. everything but real BD)
writematrix(RightFFP_NPIs.root_id, fullfile(baseDir, 'Processed_Data', 'right_FFP_root_ids.csv'));
writematrix(RightFBP_NPIs.root_id, fullfile(baseDir, 'Processed_Data', 'right_FBP_root_ids.csv'));
writematrix(RightBDP_real_NPIs.root_id, fullfile(baseDir, 'Processed_Data', 'right_BDP_root_ids.csv'));
writematrix([RightBDP_optic_NPIs.root_id; RightBDP_central_NPIs.root_id], ...
    fullfile(baseDir, 'Processed_Data', 'right_Others_root_ids.csv'));
