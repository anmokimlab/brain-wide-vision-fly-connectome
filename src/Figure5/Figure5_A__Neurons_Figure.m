
clear all; clc; close all

opt = detectImportOptions('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\fafb_v783_princeton_synapse_table.csv');
opt = setvartype(opt,'pre_root_id_720575940','int64');
opt = setvartype(opt,'post_root_id_720575940','int64');
FAFB_synapse_coordinates = readtable('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\fafb_v783_princeton_synapse_table.csv',opt);
FAFB_synapse_coordinates.Properties.VariableNames(11) = "pre_root_id";
FAFB_synapse_coordinates.Properties.VariableNames(12) = "post_root_id";
FAFB_synapse_coordinates.pre_root_id=FAFB_synapse_coordinates.pre_root_id+int64(720575940000000000);
FAFB_synapse_coordinates.post_root_id=FAFB_synapse_coordinates.post_root_id+int64(720575940000000000);

%%
close all 
clearvars -except FAFB_synapse_coordinates
clc

NeuronFaces = readmatrix('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\Figure4_BD\Neuron_Mor\720575940608367746_cLM01_faces.csv');
NeuronVertices=readmatrix('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\Figure4_BD\Neuron_Mor\720575940608367746_cLM01_vertices.csv');
NeuronFaces=NeuronFaces+1;
Neuron_root_id=int64(720575940608367746);
%{
LC9 720575940635224943
LC10b 720575940624847623
aMe8 720575940614973708
cL16 720575940633573836
LT52 720575940631938647
aMe1 720575940632780359
cLM01 720575940608367746

%}
%%%

PreSynapse_x=FAFB_synapse_coordinates.pre_x(FAFB_synapse_coordinates.pre_root_id==Neuron_root_id);
PreSynapse_y=FAFB_synapse_coordinates.pre_y(FAFB_synapse_coordinates.pre_root_id==Neuron_root_id);
PreSynapse_z=FAFB_synapse_coordinates.pre_z(FAFB_synapse_coordinates.pre_root_id==Neuron_root_id);

Output_location=[PreSynapse_x  PreSynapse_y PreSynapse_z];

%%%
PostSynapse_x=FAFB_synapse_coordinates.post_x(FAFB_synapse_coordinates.post_root_id==Neuron_root_id);
PostSynapse_y=FAFB_synapse_coordinates.post_y(FAFB_synapse_coordinates.post_root_id==Neuron_root_id);
PostSynapse_z=FAFB_synapse_coordinates.post_z(FAFB_synapse_coordinates.post_root_id==Neuron_root_id);

Input_location=[PostSynapse_x  PostSynapse_y PostSynapse_z];

%%%
% figure(1); set(gcf,'color','w')
% trisurf(NeuronFaces, NeuronVertices(:,1),NeuronVertices(:,2),NeuronVertices(:,3),'EdgeColor','none','FaceAlpha',1,'FaceColor',[0.15 0.15 0.15]); hold on;
% scatter3(Input_location(:,1),Input_location(:,2),Input_location(:,3),40,'filled','MarkerFaceAlpha',0.3)
% scatter3(Output_location(:,1),Output_location(:,2),Output_location(:,3),40,'filled','MarkerFaceAlpha',0.3)
% % set(gca,'XColor','none','YColor','none','ZColor','none')
% hold off;
% axis equal
% xlabel('x')
% ylabel('y')
% zlabel('z')
% % set(gca,'XColor','none','YColor','none','ZColor','none')
% view([0 0 1]); 
% %%%
% minSynapse=min(min(Input_location(:,1)),min(Output_location(:,1)));
% maxSynapse=max(max(Input_location(:,1)),max(Output_location(:,1)));
% 
% figure(2); set(gcf,'color','w')
% histogram(Input_location(:,1),linspace(minSynapse,maxSynapse,100),'EdgeColor','none','FaceAlpha',0.7); hold on;
% histogram(Output_location(:,1),linspace(minSynapse,maxSynapse,100),'EdgeColor','none','FaceAlpha',0.7)
% set(gca,'XColor','none','YColor','none','ZColor','none')
% % print(gcf,'-depsc2','-vector','figure1_720575940635224943_histogram.eps')

%%% PCA 해서 그려보자
Ref_x_mean=mean([Input_location(:,1);Output_location(:,1)]);
Ref_y_mean=mean([Input_location(:,2);Output_location(:,2)]);
Ref_z_mean=mean([Input_location(:,3);Output_location(:,3)]);
[coeffs, transformedData] = pca([Input_location; Output_location]);

NeuronVertices_PCA(:,1)=NeuronVertices(:,1)-Ref_x_mean;
NeuronVertices_PCA(:,2)=NeuronVertices(:,2)-Ref_y_mean;
NeuronVertices_PCA(:,3)=NeuronVertices(:,3)-Ref_z_mean;
NeuronVertices_PCA=NeuronVertices_PCA*coeffs;

Input_location_PCA(:,1)=Input_location(:,1)-Ref_x_mean;
Input_location_PCA(:,2)=Input_location(:,2)-Ref_y_mean;
Input_location_PCA(:,3)=Input_location(:,3)-Ref_z_mean;
%%% pca axis
Input_location_PCA=Input_location_PCA*coeffs;

Output_location_PCA(:,1)=Output_location(:,1)-Ref_x_mean;
Output_location_PCA(:,2)=Output_location(:,2)-Ref_y_mean;
Output_location_PCA(:,3)=Output_location(:,3)-Ref_z_mean;
%%% pca axis
Output_location_PCA=Output_location_PCA*coeffs;
%%%
figure(3); set(gcf,'color','w','Position',[444,50,1209,945])
trisurf(NeuronFaces, NeuronVertices_PCA(:,1),NeuronVertices_PCA(:,2),NeuronVertices_PCA(:,3),'EdgeColor','none','FaceAlpha',1,'FaceColor',[0.15 0.15 0.15]); hold on;
scatter3(Input_location_PCA(:,1),Input_location_PCA(:,2),Input_location_PCA(:,3),40,'filled','MarkerFaceAlpha',0.5)
scatter3(Output_location_PCA(:,1),Output_location_PCA(:,2),Output_location_PCA(:,3),40,'filled','MarkerFaceAlpha',0.5)
set(gca,'XColor','none','YColor','none','ZColor','none')
% set(gca,'TickDir','out')

hold off;
axis equal
% grid off
xlabel('x')
ylabel('y')
zlabel('z')
% view([3.769370935261506,-90]) %LC9
% view([-0.538206501722826,-63.563181048454595])% LT52
% view([-5.931305284241514,21.714545459478511]); %cL16
view([0 -90]) %cLM01
% view([-1.307161461171119,-14.932896704854537]); %aMe1
 % print(gcf,'-depsc2','-vector','figure1_720575940640851363_neuron.eps')

 % yz 평면에서 보기
%%%
minSynapse=min(min(Input_location_PCA(:,1)),min(Output_location_PCA(:,1)));
maxSynapse=max(max(Input_location_PCA(:,1)),max(Output_location_PCA(:,1)));

figure(4); set(gcf,'color','w')
histogram(Input_location_PCA(:,1),linspace(minSynapse,maxSynapse,100),'EdgeColor','none','FaceAlpha',0.7); hold on;
histogram(Output_location_PCA(:,1),linspace(minSynapse,maxSynapse,100),'EdgeColor','none','FaceAlpha',0.7)
set(gca,'XColor','none','YColor','none','ZColor','none')
% print(gcf,'-depsc2','-vector','LC10b_histogram.eps')

%%% PCA 적용된 xyz 좌표축 추가
origin = [0, 0, 0]; % 원점
axis_length = 1000; % 축 길이

% 기존 XYZ 축을 PCA 좌표계로 변환
x_axis_pca = [axis_length, 0, 0] * coeffs;
y_axis_pca = [0, axis_length, 0] * coeffs;
z_axis_pca = [0, 0, axis_length] * coeffs;

figure(5);hold on;
quiver3(origin(1), origin(2), origin(3), x_axis_pca(1), x_axis_pca(2), x_axis_pca(3), 'r', 'LineWidth', 2, 'MaxHeadSize', 1); % X축
quiver3(origin(1), origin(2), origin(3), y_axis_pca(1), y_axis_pca(2), y_axis_pca(3), 'g', 'LineWidth', 2, 'MaxHeadSize', 1); % Y축
quiver3(origin(1), origin(2), origin(3), z_axis_pca(1), z_axis_pca(2), z_axis_pca(3), 'b', 'LineWidth', 2, 'MaxHeadSize', 1); % Z축
legend({'X-axis (PCA)', 'Y-axis (PCA)', 'Z-axis (PCA)'},'Location','northeastoutside');

hold off;
% view([0 -1 0]); 
% view([3.769370935261506,-90]) %LC9
% view([-0.538206501722826,-63.563181048454595])% LT52
% view([-5.931305284241514,21.714545459478511]); %cL16
view([0 -90]) %cLM01
% view([-1.307161461171119,-14.932896704854537]); %aMe1
% xlim([-2000 2000])
% ylim([-2000 2000])
% zlim([-2000 2000])
axis equal
% print(gcf,'-depsc2','-vector','figure1_720575940635224943_arrow.eps')
