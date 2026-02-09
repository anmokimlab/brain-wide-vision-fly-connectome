%% 데이터 가져오기 / synapse / cL rootids
clear all; clc; close all


opt = detectImportOptions('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\fafb_v783_princeton_synapse_table.csv');
opt = setvartype(opt,'pre_root_id_720575940','int64');
opt = setvartype(opt,'post_root_id_720575940','int64');
FAFB_synapse_coordinates = readtable('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\fafb_v783_princeton_synapse_table.csv',opt);
FAFB_synapse_coordinates.Properties.VariableNames(11) = "pre_root_id";
FAFB_synapse_coordinates.Properties.VariableNames(12) = "post_root_id";
FAFB_synapse_coordinates.pre_root_id=FAFB_synapse_coordinates.pre_root_id+int64(720575940000000000);
FAFB_synapse_coordinates.post_root_id=FAFB_synapse_coordinates.post_root_id+int64(720575940000000000);
FAFB_synapse_coordinates.pre_x=FAFB_synapse_coordinates.pre_x*1e-9;
FAFB_synapse_coordinates.pre_y=FAFB_synapse_coordinates.pre_y*1e-9;
FAFB_synapse_coordinates.pre_z=FAFB_synapse_coordinates.pre_z*1e-9;

opt = detectImportOptions('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\consolidated_cell_types.csv');
opt = setvartype(opt,'root_id','int64');
FAFBConsolidated_type = readtable('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\consolidated_cell_types.csv',opt);


LOP_R_Faces = readmatrix('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\Neuropil\LoP_R_faces.csv');
LOP_R_Vertices=readmatrix('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\Neuropil\LoP_R_vertices.csv');
LOP_R_Faces=LOP_R_Faces+1;
LOP_R_Vertices=LOP_R_Vertices*1e-9;

LOP_L_Faces = readmatrix('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\Neuropil\LoP_L_faces.csv');
LOP_L_Vertices=readmatrix('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\Neuropil\LoP_L_vertices.csv');
LOP_L_Faces=LOP_L_Faces+1;
LOP_L_Vertices=LOP_L_Vertices*1e-9;

LO_R_Faces = readmatrix('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\Neuropil\Lo_R_faces.csv');
LO_R_Vertices=readmatrix('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\Neuropil\Lo_R_vertices.csv');
LO_R_Faces=LO_R_Faces+1;
LO_R_Vertices=LO_R_Vertices*1e-9;

LO_L_Faces = readmatrix('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\Neuropil\Lo_L_faces.csv');
LO_L_Vertices=readmatrix('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\Neuropil\Lo_L_vertices.csv');
LO_L_Faces=LO_L_Faces+1;
LO_L_Vertices=LO_L_Vertices*1e-9;

ME_R_Faces = readmatrix('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\Neuropil\Me_R_faces.csv');
ME_R_Vertices=readmatrix('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\Neuropil\Me_R_vertices.csv');
ME_R_Faces=ME_R_Faces+1;
ME_R_Vertices=ME_R_Vertices*1e-9;

ME_L_Faces = readmatrix('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\Neuropil\Me_L_faces.csv');
ME_L_Vertices=readmatrix('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\Neuropil\Me_L_vertices.csv');
ME_L_Faces=ME_L_Faces+1;
ME_L_Vertices=ME_L_Vertices*1e-9;
%% custom color
% figure(1);set(gcf,'Color','w')
% colormap의 단계 수 설정
n = 256;

% 시작 색 (검은색)과 종료 색 ([0.8500, 0.3250, 0.0980]) 정의
startColor = [0, 0, 0];
endColor = [0.8500, 0.3250, 0.0980];

% 색상 보간 (선형 보간)
customColormapOut = zeros(n, 3);
for i = 1:3
    customColormapOut(:, i) = linspace(startColor(i), endColor(i), n);
end

% colormap 적용
colormap(customColormapOut);

% 색상 막대 표시 (테스트용)
% colorbar;
%
% colormap의 단계 수 설정
n = 256;

% 시작 색 (검은색)과 종료 색 ([0.8500, 0.3250, 0.0980]) 정의
startColor = [0, 0, 0];
endColor = [0 0.4470 0.7410];

% 색상 보간 (선형 보간)
customColormapIn = zeros(n, 3);
for i = 1:3
    customColormapIn(:, i) = linspace(startColor(i), endColor(i), n);
end

% colormap 적용
colormap(customColormapIn);

% 색상 막대 표시 (테스트용)
% colorbar;

% colormap의 단계 수 설정
n = 256;

% 시작 색 (검은색)과 종료 색 ([0.8500, 0.3250, 0.0980]) 정의
startColor = [0, 0, 0];
endColor = [0.4660 0.6740 0.1880]	;

% 색상 보간 (선형 보간)
customColormapIn_contra = zeros(n, 3);
for i = 1:3
    customColormapIn_contra(:, i) = linspace(startColor(i), endColor(i), n);
end

% colormap 적용
colormap(customColormapIn_contra);
colorbar
%%
load('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\Figure1_Overall\Right_Neurons_Thr0.mat','RightFB_NPIs')
load('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\Figure1_Overall\FAFB_NPI_Thr0.mat')

[RightFB_type,~,ic]=unique(RightFB_NPIs.type);

RightFB_type=table(RightFB_type,'VariableNames',{'type'});

for i=1:1:size(RightFB_type,1)
    idx=ic==i;
    RightFB_type.root_id{i}=RightFB_NPIs.root_id(idx);
end

Me_FB_idx=[4 22 25 26 29 30 31 32 33 68 69 70 71 72 73 74 75 76 78 81 82 83 84 85 86 87 93];
Lo_FB_idx=[2 10 11 12 13 14 15 16 17 18 19 21 34 35 36 37 38 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 95];
Lop_FB_idx=[3 7 8 9 28 62 63 64 65 66];
Multi_FB_idx=[1 5 6 23 24 27 39 60 61 67 77 79 80 88 89 90 91 92 94 96];



FB_Me=RightFB_type(Me_FB_idx,:);
FB_Lo=RightFB_type(Lo_FB_idx,:);
FB_Lop=RightFB_type(Lop_FB_idx,:);
FB_Multi=RightFB_type(Multi_FB_idx,:);

Target_Neuron=FB_Lo.type;
Target_What='axon';
Target_Where='Lo';

%% Ref and Target root ids;
Ref_Neuron={'LT1a'}; %% Dm6 dendrite LT1a dendrite T5a axon
Ref_What='dendrite';
Ref_Where='Lo';

Ref_Table=table(Ref_Neuron,'VariableNames',{'cell_type'});

for i=1:1:size(Ref_Table,1)
    idx=strcmpi(FAFBConsolidated_type.primary_type,Ref_Table.cell_type{i});
    Wantrootids=FAFBConsolidated_type.root_id(idx);
    Ref_Table.root_ids{i}=Wantrootids;
    Ref_Table.N{i}=sum(idx);
end

Target_Table=table(Target_Neuron,'VariableNames',{'cell_type'});

for i=1:1:size(Target_Table,1)
    idx=strcmpi(RightFB_type.type,Target_Table.cell_type{i});
    Wantrootids=RightFB_type.root_id{idx};
    Target_Table.root_ids{i}=Wantrootids;
    Target_Table.N{i}=size(Wantrootids,1);
end


%% Ref synapse 위치
for i=1:1:size(Ref_Table,1)
    ref_syn_loc=[];
    if strcmpi(Ref_What,'axon')
        Ref_Syn_pre_root_idx=find(ismember(FAFB_synapse_coordinates.pre_root_id,Ref_Table.root_ids{i}));
        ref_syn_loc= FAFB_synapse_coordinates(Ref_Syn_pre_root_idx,1:3);
        
    elseif strcmpi(Ref_What,'dendrite')
        Ref_Syn_post_root_idx=find(ismember(FAFB_synapse_coordinates.post_root_id,Ref_Table.root_ids{i}));
        ref_syn_loc=FAFB_synapse_coordinates(Ref_Syn_post_root_idx,1:3);
    else
        msg = 'axon or dendrite';
        error(msg)
    end

    ref_syn_loc=ref_syn_loc{:,:};
    if strcmpi(Ref_Where,'Lop')
        Ref_syn_loc_in_R=ref_syn_loc(intriangulation(LOP_R_Vertices,LOP_R_Faces,ref_syn_loc),:);
        Ref_syn_loc_in_L=ref_syn_loc(intriangulation(LOP_L_Vertices,LOP_L_Faces,ref_syn_loc),:);
    elseif strcmpi(Ref_Where,'Lo')
        Ref_syn_loc_in_R=ref_syn_loc(intriangulation(LO_R_Vertices,LO_R_Faces,ref_syn_loc),:);
        Ref_syn_loc_in_L=ref_syn_loc(intriangulation(LO_L_Vertices,LO_L_Faces,ref_syn_loc),:);
    elseif strcmpi(Ref_Where,'Me')
        Ref_syn_loc_in_R=ref_syn_loc(intriangulation(ME_R_Vertices,ME_R_Faces,ref_syn_loc),:);
        Ref_syn_loc_in_L=ref_syn_loc(intriangulation(ME_L_Vertices,ME_L_Faces,ref_syn_loc),:);
    else
        msg = 'Lop or Lo or Me';
        error(msg)
    end
    %%%% outlier 제거
    Ref_syn_loc_in_L_pointCloud = pointCloud(Ref_syn_loc_in_L);
    Ref_syn_loc_in_L_pointCloud_Denoise = pcdenoise(Ref_syn_loc_in_L_pointCloud,"NumNeighbors",20,"Threshold",2.5);
    Ref_syn_loc_in_L_Denoise=Ref_syn_loc_in_L_pointCloud_Denoise.Location;

    Ref_syn_loc_in_R_pointCloud = pointCloud(Ref_syn_loc_in_R);
    Ref_syn_loc_in_R_pointCloud_Denoise = pcdenoise(Ref_syn_loc_in_R_pointCloud,"NumNeighbors",20,"Threshold",2.5);
    Ref_syn_loc_in_R_Denoise=Ref_syn_loc_in_R_pointCloud_Denoise.Location;

    Ref_Table.syn_loc{i}=ref_syn_loc;
    Ref_Table.syn_loc_L{i}=Ref_syn_loc_in_L;
    Ref_Table.syn_loc_R{i}=Ref_syn_loc_in_R;
    Ref_Table.syn_loc_L_Denoise{i}=Ref_syn_loc_in_L_Denoise;
    Ref_Table.syn_loc_R_Denoise{i}=Ref_syn_loc_in_R_Denoise;
end

%% PCA 해보자
Ref_x_mean_R=mean(Ref_Table.syn_loc_R_Denoise{1}(:,1));
Ref_y_mean_R=mean(Ref_Table.syn_loc_R_Denoise{1}(:,2));
Ref_z_mean_R=mean(Ref_Table.syn_loc_R_Denoise{1}(:,3));
[coeffs_R, transformedData_R] = pca(Ref_Table.syn_loc_R_Denoise{1});

Ref_x_mean_L=mean(Ref_Table.syn_loc_L_Denoise{1}(:,1));
Ref_y_mean_L=mean(Ref_Table.syn_loc_L_Denoise{1}(:,2));
Ref_z_mean_L=mean(Ref_Table.syn_loc_L_Denoise{1}(:,3));
[coeffs_L, transformedData_L] = pca(Ref_Table.syn_loc_L_Denoise{1});
reverse='false';
if strcmpi(reverse,'true')
    transformedData_L(:,3) = -transformedData_L(:,3);   % 필요할 쪽만(R 또는 L)
    coeffs_L(:,3) = -coeffs_L(:,3);
end
%% fitting parabola fitting

[xData_R, yData_R, zData_R] = prepareSurfaceData(transformedData_R(:,1), transformedData_R(:,2), transformedData_R(:,3));
[xData_L, yData_L, zData_L] = prepareSurfaceData(transformedData_L(:,1), transformedData_L(:,2), transformedData_L(:,3));

if strcmpi(Ref_Where,'Lop')
    ft = fittype( 'poly33' );

    % 데이터에 모델을 피팅하십시오.
    [fitresult_R, gof_R] = fit( [xData_R, yData_R], zData_R, ft, 'Normalize', 'off'  );

    % % Define function representing the paraboloid's surface
    paraboloid_surface_R = @(x, y) fitresult_R.p00+fitresult_R.p10*x+fitresult_R.p01*y...
        +fitresult_R.p20*(x.^2)+fitresult_R.p11*x.*y+fitresult_R.p02*(y.^2)...
        +fitresult_R.p30*(x.^3)+fitresult_R.p21*(x.^2.*y)+fitresult_R.p12*(x.*y.^2)+fitresult_R.p03*(y.^3);
    fprintf('R-squared_R: %f\n', gof_R.adjrsquare);

    [fitresult_L, gof_L] = fit( [xData_L, yData_L], zData_L, ft, 'Normalize', 'off'  );

    % % Define function representing the paraboloid's surface
    paraboloid_surface_L = @(x, y) fitresult_L.p00+fitresult_L.p10*x+fitresult_L.p01*y...
        +fitresult_L.p20*(x.^2)+fitresult_L.p11*x.*y+fitresult_L.p02*(y.^2)...
        +fitresult_L.p30*(x.^3)+fitresult_L.p21*(x.^2.*y)+fitresult_L.p12*(x.*y.^2)+fitresult_L.p03*(y.^3);
    fprintf('R-squared_L: %f\n', gof_L.adjrsquare);


elseif strcmpi(Ref_Where,'Lo')
    % fittype과 옵션을 설정하십시오.
    ft = fittype( 'poly22' );

    % 데이터에 모델을 피팅하십시오.
    [fitresult_R, gof_R] = fit( [xData_R, yData_R], zData_R, ft, 'Normalize', 'off'  );

    % % Define function representing the paraboloid's surface
    paraboloid_surface_R = @(x, y) fitresult_R.p00+fitresult_R.p10*x+fitresult_R.p01*y...
        +fitresult_R.p20*x.^2+fitresult_R.p11*x.*y+fitresult_R.p02*y.^2;
    fprintf('R-squared_R: %f\n', gof_R.adjrsquare);

    % 데이터에 모델을 피팅하십시오.
    [fitresult_L, gof_L] = fit( [xData_L, yData_L], zData_L, ft, 'Normalize', 'off'  );

    % % Define function representing the paraboloid's surface
    paraboloid_surface_L = @(x, y) fitresult_L.p00+fitresult_L.p10*x+fitresult_L.p01*y...
        +fitresult_L.p20*x.^2+fitresult_L.p11*x.*y+fitresult_L.p02*y.^2;
    fprintf('R-squared_L: %f\n', gof_L.adjrsquare);
elseif strcmpi(Ref_Where,'me')
    ft = fittype( 'poly22' );

    % 데이터에 모델을 피팅하십시오.
    [fitresult_R, gof_R] = fit( [xData_R, yData_R], zData_R, ft, 'Normalize', 'off'  );

    % % Define function representing the paraboloid's surface
    paraboloid_surface_R = @(x, y) fitresult_R.p00+fitresult_R.p10*x+fitresult_R.p01*y...
        +fitresult_R.p20*x.^2+fitresult_R.p11*x.*y+fitresult_R.p02*y.^2;
    fprintf('R-squared_R: %f\n', gof_R.adjrsquare);

    % 데이터에 모델을 피팅하십시오.
    [fitresult_L, gof_L] = fit( [xData_L, yData_L], zData_L, ft, 'Normalize', 'off'  );

    % % Define function representing the paraboloid's surface
    paraboloid_surface_L = @(x, y) fitresult_L.p00+fitresult_L.p10*x+fitresult_L.p01*y...
        +fitresult_L.p20*x.^2+fitresult_L.p11*x.*y+fitresult_L.p02*y.^2;
    fprintf('R-squared_L: %f\n', gof_L.adjrsquare);

else
    msg = 'Lop or Lo or Me';
    error(msg)
end
%%
close all
figure(1);set(gcf,'Color','w');
hold on;
[X,Y] = meshgrid(min(transformedData_R(:,1)):5e-7:max(transformedData_R(:,1)),min(transformedData_R(:,2)):5e-7:max(transformedData_R(:,2)));
Z = paraboloid_surface_R(X,Y);
surf(X,Y,Z,'EdgeColor','none')
scatter3(transformedData_R(:,1),transformedData_R(:,2),transformedData_R(:,3),20,[0.8500 0.3250 0.0980],'filled');
% scatter3(transformedData_R(:,1),transformedData_R(:,2),transformedData_R(:,3),20,[0 0.4470 0.7410]	,'filled'); hold on;

grid on;

FaceAlhpaValue=0.7;
if strcmpi(Ref_Where,'Me')
    ME_R_Vertices_PCA=ME_R_Vertices;
    ME_R_Vertices_PCA(:,1)=ME_R_Vertices_PCA(:,1)-Ref_x_mean_R;
    ME_R_Vertices_PCA(:,2)=ME_R_Vertices_PCA(:,2)-Ref_y_mean_R;
    ME_R_Vertices_PCA(:,3)=ME_R_Vertices_PCA(:,3)-Ref_z_mean_R;
    ME_R_Vertices_PCA=ME_R_Vertices_PCA*coeffs_R;
    trisurf(ME_R_Faces,ME_R_Vertices_PCA(:,1),ME_R_Vertices_PCA(:,2),ME_R_Vertices_PCA(:,3),'FaceAlpha',0.2,'EdgeColor','None');
elseif strcmpi(Ref_Where,'Lo')
    LO_R_Vertices_PCA=LO_R_Vertices;
    LO_R_Vertices_PCA(:,1)=LO_R_Vertices_PCA(:,1)-Ref_x_mean_R;
    LO_R_Vertices_PCA(:,2)=LO_R_Vertices_PCA(:,2)-Ref_y_mean_R;
    LO_R_Vertices_PCA(:,3)=LO_R_Vertices_PCA(:,3)-Ref_z_mean_R;
    LO_R_Vertices_PCA=LO_R_Vertices_PCA*coeffs_R;
    trisurf(LO_R_Faces,LO_R_Vertices_PCA(:,1),LO_R_Vertices_PCA(:,2),LO_R_Vertices_PCA(:,3),'FaceAlpha',0.2,'EdgeColor','None');
elseif strcmpi(Ref_Where,'Lop')
    LOP_R_Vertices_PCA=LOP_R_Vertices;
    LOP_R_Vertices_PCA(:,1)=LOP_R_Vertices_PCA(:,1)-Ref_x_mean_R;
    LOP_R_Vertices_PCA(:,2)=LOP_R_Vertices_PCA(:,2)-Ref_y_mean_R;
    LOP_R_Vertices_PCA(:,3)=LOP_R_Vertices_PCA(:,3)-Ref_z_mean_R;
    LOP_R_Vertices_PCA=LOP_R_Vertices_PCA*coeffs_R;
    trisurf(LOP_R_Faces,LOP_R_Vertices_PCA(:,1),LOP_R_Vertices_PCA(:,2),LOP_R_Vertices_PCA(:,3),'FaceAlpha',0.2,'EdgeColor','None');
else
    msg = 'Lop or Lo or Me';
    error(msg)
end

axis("equal")
xlabel( 'PC1', 'Interpreter', 'none' );
ylabel( 'PC2', 'Interpreter', 'none' );
zlabel( 'PC3', 'Interpreter', 'none' );
view([-31.459612469645986,5.86855640905495]) %%%%% me
% view([-59.6596132214191,9.018885905511171]) %%%%%%lop

ax1 = gca;  % 현재 figure(1)의 축 가져오기
x_limits = xlim(ax1);
y_limits = ylim(ax1);
z_limits = zlim(ax1);

figure(2);set(gcf,'Color','w');
hold on;
xyz_axes = eye(3); % 단위 행렬 (각 축을 나타냄)
arrow_length = 2e-5; % 2 * 10^(-5)

new_axes = arrow_length * (xyz_axes * coeffs_R); % 미리 스케일 적용
% 기존 XYZ 축을 새로운 PC 공간에 플로팅
quiver3(0, 0, 0, new_axes(1,1), new_axes(2,1), new_axes(3,1), 'r', 'LineWidth', 2); % X축
quiver3(0, 0, 0, new_axes(1,2), new_axes(2,2), new_axes(3,2), 'g', 'LineWidth', 2); % Y축
quiver3(0, 0, 0, new_axes(1,3), new_axes(2,3), new_axes(3,3), 'b', 'LineWidth', 2); % Z축

% 축 라벨링
% text(new_axes(1,1), new_axes(2,1), new_axes(3,1), 'X', 'Color', 'r', 'FontSize', 12);
% text(new_axes(1,2), new_axes(2,2), new_axes(3,2), 'Y', 'Color', 'g', 'FontSize', 12);
% text(new_axes(1,3), new_axes(2,3), new_axes(3,3), 'Z', 'Color', 'b', 'FontSize', 12);
ax2 = gca;  % 현재 figure(2)의 축 가져오기
xlim(ax2, x_limits);
ylim(ax2, y_limits);
zlim(ax2, z_limits);
% axis("equal")
xlabel( 'PC1', 'Interpreter', 'none' );
ylabel( 'PC2', 'Interpreter', 'none' );
zlabel( 'PC3', 'Interpreter', 'none' );
view([-31.459612469645986,5.86855640905495]) %%%%% me
% view([-59.6596132214191,9.018885905511171]) %%%%%%lop
%%
figure(3); set(gcf, 'Color', 'w');
hold on;
[X,Y] = meshgrid(min(transformedData_L(:,1)):5e-7:max(transformedData_L(:,1)), min(transformedData_L(:,2)):5e-7:max(transformedData_L(:,2)));
Z = paraboloid_surface_L(X,Y);
surf(X,Y,Z,'EdgeColor','none')
scatter3(transformedData_L(:,1), transformedData_L(:,2), transformedData_L(:,3), 20, [0.8500 0.3250 0.0980], 'filled');

grid on;

FaceAlphaValue = 0.7;
if strcmpi(Ref_Where,'Me')
    ME_L_Vertices_PCA = ME_L_Vertices;
    ME_L_Vertices_PCA(:,1) = ME_L_Vertices_PCA(:,1) - Ref_x_mean_L;
    ME_L_Vertices_PCA(:,2) = ME_L_Vertices_PCA(:,2) - Ref_y_mean_L;
    ME_L_Vertices_PCA(:,3) = ME_L_Vertices_PCA(:,3) - Ref_z_mean_L;
    ME_L_Vertices_PCA = ME_L_Vertices_PCA * coeffs_L;
    trisurf(ME_L_Faces, ME_L_Vertices_PCA(:,1), ME_L_Vertices_PCA(:,2), ME_L_Vertices_PCA(:,3), 'FaceAlpha', 0.2, 'EdgeColor', 'None');
elseif strcmpi(Ref_Where,'Lo')
    LO_L_Vertices_PCA = LO_L_Vertices;
    LO_L_Vertices_PCA(:,1) = LO_L_Vertices_PCA(:,1) - Ref_x_mean_L;
    LO_L_Vertices_PCA(:,2) = LO_L_Vertices_PCA(:,2) - Ref_y_mean_L;
    LO_L_Vertices_PCA(:,3) = LO_L_Vertices_PCA(:,3) - Ref_z_mean_L;
    LO_L_Vertices_PCA = LO_L_Vertices_PCA * coeffs_L;
    trisurf(LO_L_Faces, LO_L_Vertices_PCA(:,1), LO_L_Vertices_PCA(:,2), LO_L_Vertices_PCA(:,3), 'FaceAlpha', 0.2, 'EdgeColor', 'None');
elseif strcmpi(Ref_Where,'Lop')
    LOP_L_Vertices_PCA = LOP_L_Vertices;
    LOP_L_Vertices_PCA(:,1) = LOP_L_Vertices_PCA(:,1) - Ref_x_mean_L;
    LOP_L_Vertices_PCA(:,2) = LOP_L_Vertices_PCA(:,2) - Ref_y_mean_L;
    LOP_L_Vertices_PCA(:,3) = LOP_L_Vertices_PCA(:,3) - Ref_z_mean_L;
    LOP_L_Vertices_PCA = LOP_L_Vertices_PCA * coeffs_L;
    trisurf(LOP_L_Faces, LOP_L_Vertices_PCA(:,1), LOP_L_Vertices_PCA(:,2), LOP_L_Vertices_PCA(:,3), 'FaceAlpha', 0.2, 'EdgeColor', 'None');
else
    error('Lop or Lo or Me');
end

axis("equal")
xlabel('PC1', 'Interpreter', 'none');
ylabel('PC2', 'Interpreter', 'none');
zlabel('PC3', 'Interpreter', 'none');
view([-31.459612469645986,5.86855640905495]) %%%%% me

% view([-59.6596132214191,9.018885905511171]) %%%%%%lop

ax1 = gca;
x_limits = xlim(ax1);
y_limits = ylim(ax1);
z_limits = zlim(ax1);

figure(4); set(gcf,'Color','w');
hold on;
xyz_axes = eye(3);
arrow_length = 2e-5;

new_axes = arrow_length * (xyz_axes * coeffs_L);
quiver3(0, 0, 0, new_axes(1,1), new_axes(2,1), new_axes(3,1), 'r', 'LineWidth', 2);
quiver3(0, 0, 0, new_axes(1,2), new_axes(2,2), new_axes(3,2), 'g', 'LineWidth', 2);
quiver3(0, 0, 0, new_axes(1,3), new_axes(2,3), new_axes(3,3), 'b', 'LineWidth', 2);

ax2 = gca;
xlim(ax2, x_limits);
ylim(ax2, y_limits);
zlim(ax2, z_limits);
xlabel('PC1', 'Interpreter', 'none');
ylabel('PC2', 'Interpreter', 'none');
zlabel('PC3', 'Interpreter', 'none');
view([-31.459612469645986,5.86855640905495]) %%%%% me

% view([-59.6596132214191,9.018885905511171]) %%%%%%lop

%% Target Synapse위치

for i=1:1:size(Target_Table,1)
    target_syn_loc=[];
    if strcmpi(Target_What,'axon')

        Target_Syn_pre_root_idx=find(ismember(FAFB_synapse_coordinates.pre_root_id,Target_Table.root_ids{i}));
        target_syn_loc=FAFB_synapse_coordinates(Target_Syn_pre_root_idx,1:3);

    elseif strcmpi(Target_What,'dendrite')

        Target_Syn_post_root_idx=find(ismember(FAFB_synapse_coordinates.post_root_id,Target_Table.root_ids{i}));
        target_syn_loc=FAFB_synapse_coordinates(Target_Syn_post_root_idx,1:3);
    else
        msg = 'axon or dendrite';
        error(msg)
    end
    if isempty(target_syn_loc)
        continue;
    end
    target_syn_loc=target_syn_loc{:,:};
    if strcmpi(Target_Where,'Lop')
        Target_syn_loc_in_R=target_syn_loc(intriangulation(LOP_R_Vertices,LOP_R_Faces,target_syn_loc),:);
        Target_syn_loc_in_L=target_syn_loc(intriangulation(LOP_L_Vertices,LOP_L_Faces,target_syn_loc),:);
    elseif strcmpi(Target_Where,'Lo')
        Target_syn_loc_in_R=target_syn_loc(intriangulation(LO_R_Vertices,LO_R_Faces,target_syn_loc),:);
        Target_syn_loc_in_L=target_syn_loc(intriangulation(LO_L_Vertices,LO_L_Faces,target_syn_loc),:);
    elseif strcmpi(Target_Where,'Me')
        Target_syn_loc_in_R=target_syn_loc(intriangulation(ME_R_Vertices,ME_R_Faces,target_syn_loc),:);
        Target_syn_loc_in_L=target_syn_loc(intriangulation(ME_L_Vertices,ME_L_Faces,target_syn_loc),:);
    else
        msg = 'Lop or Lo or Me';
        error(msg)
    end
    Target_Table.syn_loc{i}=target_syn_loc;

    %%%% outlier 제거
    if ~isempty(Target_syn_loc_in_L)
        Target_syn_loc_in_L_pointCloud = pointCloud(Target_syn_loc_in_L);
        Target_syn_loc_in_L_pointCloud_Denoise = pcdenoise(Target_syn_loc_in_L_pointCloud,"NumNeighbors",20,"Threshold",2.5);
        Target_syn_loc_in_L_Denoise=Target_syn_loc_in_L_pointCloud_Denoise.Location;
        Target_Table.syn_loc_L{i}=Target_syn_loc_in_L;
        Target_Table.syn_loc_L_Denoise{i}=Target_syn_loc_in_L_Denoise;
    else
        Target_Table.syn_loc_L{i}=[];
        Target_Table.syn_loc_L_Denoise{i}=[];
    end
    if ~isempty(Target_syn_loc_in_R)
        Target_syn_loc_in_R_pointCloud = pointCloud(Target_syn_loc_in_R);
        Target_syn_loc_in_R_pointCloud_Denoise = pcdenoise(Target_syn_loc_in_R_pointCloud,"NumNeighbors",20,"Threshold",2.5);
        Target_syn_loc_in_R_Denoise=Target_syn_loc_in_R_pointCloud_Denoise.Location;
        Target_Table.syn_loc_R{i}=Target_syn_loc_in_R;
        Target_Table.syn_loc_R_Denoise{i}=Target_syn_loc_in_R_Denoise;
    else
        Target_Table.syn_loc_R{i}=[];
        Target_Table.syn_loc_R_Denoise{i}=[];
    end
    
    Target_syn_loc_in_L_Shift=Target_syn_loc_in_L; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% center shift
    Target_syn_loc_in_L_Shift(:,1)=Target_syn_loc_in_L_Shift(:,1)-Ref_x_mean_L;
    Target_syn_loc_in_L_Shift(:,2)=Target_syn_loc_in_L_Shift(:,2)-Ref_y_mean_L;
    Target_syn_loc_in_L_Shift(:,3)=Target_syn_loc_in_L_Shift(:,3)-Ref_z_mean_L;
    %%% pca axis
    Target_syn_loc_in_L_Shift_PCA=Target_syn_loc_in_L_Shift*coeffs_L;
    Target_Table.syn_loc_L_PCA{i}=Target_syn_loc_in_L_Shift_PCA;

    Target_syn_loc_in_R_Shift=Target_syn_loc_in_R; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% center shift
    Target_syn_loc_in_R_Shift(:,1)=Target_syn_loc_in_R_Shift(:,1)-Ref_x_mean_R;
    Target_syn_loc_in_R_Shift(:,2)=Target_syn_loc_in_R_Shift(:,2)-Ref_y_mean_R;
    Target_syn_loc_in_R_Shift(:,3)=Target_syn_loc_in_R_Shift(:,3)-Ref_z_mean_R;
    %%% pca axis
    Target_syn_loc_in_R_Shift_PCA=Target_syn_loc_in_R_Shift*coeffs_R;
    Target_Table.syn_loc_R_PCA{i}=Target_syn_loc_in_R_Shift_PCA;



    % % Given point

    z_hat_R=paraboloid_surface_R(Target_syn_loc_in_R_Shift_PCA(:,1),Target_syn_loc_in_R_Shift_PCA(:,2));
    z_distribution_R=Target_syn_loc_in_R_Shift_PCA(:,3)-z_hat_R;

    Target_Table.syn_zDist_from_Parabola_R{i}=z_distribution_R;

    z_hat_L=paraboloid_surface_L(Target_syn_loc_in_L_Shift_PCA(:,1),Target_syn_loc_in_L_Shift_PCA(:,2));
    z_distribution_L=Target_syn_loc_in_L_Shift_PCA(:,3)-z_hat_L;

    Target_Table.syn_zDist_from_Parabola_L{i}=z_distribution_L;
end

%%

edges = -7e-5:1e-6:5e-5;
synapseBox_z_R=zeros(size(edges,2)-1,size(Target_Table,1));
synapseBox_z_L=zeros(size(edges,2)-1,size(Target_Table,1));

for i=1:1:size(Target_Table,1)
    temp=Target_Table.syn_zDist_from_Parabola_R{i};
    [counts, ~] = histcounts(temp, 'BinEdges', edges);
    syn_histcounts = counts / max(counts);
    synapseBox_z_R(:,i)=syn_histcounts;

    temp=Target_Table.syn_zDist_from_Parabola_L{i};
    [counts, ~] = histcounts(temp, 'BinEdges', edges);
    syn_histcounts = counts / max(counts);
    synapseBox_z_L(:,i)=syn_histcounts;
end
%%
figure(5); set(gcf,'Color','w')
imagesc(synapseBox_z_R)
set(gca,'YTick',0.5:1:size(synapseBox_z_R,1)+0.5,'YTickLabel',edges,'XTick',1:1:size(Target_Table,1),'XTickLabel',Target_Table.cell_type,'TickDir','out','Box','off','YDir','normal')
xlabel( 'Neurons');
ylabel( 'Innvervaiton depth (μm)');

colormap(customColormapOut)
% ylim([-5e-5 3e-5])
title('R')
% grid on
xtickangle(45)
%%
figure(6); set(gcf,'Color','w')
imagesc(synapseBox_z_L)
set(gca,'YTick',0.5:1:size(synapseBox_z_L,1)+0.5,'YTickLabel',edges,'XTick',1:1:size(Target_Table,1),'XTickLabel',Target_Table.cell_type,'TickDir','out','Box','off','YDir','normal')
xlabel( 'Neurons');
ylabel( 'Innvervaiton depth (μm)');
colormap(customColormapOut)
% ylim([-5e-5 3e-5])
title('L')
% grid on
xtickangle(45)
%% 하나씩 경우 체크

% f3=figure(3);set(gcf,'Color','w'); hold on;
% [X,Y] = meshgrid(min(transformedData_R(:,1)):5e-7:max(transformedData_R(:,1)),min(transformedData_R(:,2)):5e-7:max(transformedData_R(:,2)));
% Z = paraboloid_surface(X,Y);
% surf(X,Y,Z,'EdgeColor','none')
% scatter3(transformedData_R(:,1),transformedData_R(:,2),transformedData_R(:,3)); hold on;
%
% FaceAlhpaValue=0.7;
% if strcmpi(Ref_Where,'Me')
%     ME_R_Vertices_PCA=ME_R_Vertices;
%     ME_R_Vertices_PCA(:,1)=ME_R_Vertices_PCA(:,1)-Ref_x_mean;
%     ME_R_Vertices_PCA(:,2)=ME_R_Vertices_PCA(:,2)-Ref_y_mean;
%     ME_R_Vertices_PCA(:,3)=ME_R_Vertices_PCA(:,3)-Ref_z_mean;
%     ME_R_Vertices_PCA=ME_R_Vertices_PCA*coeffs_R;
%     trisurf(ME_R_Faces,ME_R_Vertices_PCA(:,1),ME_R_Vertices_PCA(:,2),ME_R_Vertices_PCA(:,3),'FaceAlpha',0.2,'EdgeColor','None');
% elseif strcmpi(Ref_Where,'Lo')
%     LO_R_Vertices_PCA=LO_R_Vertices;
%     LO_R_Vertices_PCA(:,1)=LO_R_Vertices_PCA(:,1)-Ref_x_mean;
%     LO_R_Vertices_PCA(:,2)=LO_R_Vertices_PCA(:,2)-Ref_y_mean;
%     LO_R_Vertices_PCA(:,3)=LO_R_Vertices_PCA(:,3)-Ref_z_mean;
%     LO_R_Vertices_PCA=LO_R_Vertices_PCA*coeffs_R;
%     trisurf(LO_R_Faces,LO_R_Vertices_PCA(:,1),LO_R_Vertices_PCA(:,2),LO_R_Vertices_PCA(:,3),'FaceAlpha',0.2,'EdgeColor','None');
% elseif strcmpi(Ref_Where,'Lop')
%     LOP_R_Vertices_PCA=LOP_R_Vertices;
%     LOP_R_Vertices_PCA(:,1)=LOP_R_Vertices_PCA(:,1)-Ref_x_mean;
%     LOP_R_Vertices_PCA(:,2)=LOP_R_Vertices_PCA(:,2)-Ref_y_mean;
%     LOP_R_Vertices_PCA(:,3)=LOP_R_Vertices_PCA(:,3)-Ref_z_mean;
%     LOP_R_Vertices_PCA=LOP_R_Vertices_PCA*coeffs_R;
%     trisurf(LOP_R_Faces,LOP_R_Vertices_PCA(:,1),LOP_R_Vertices_PCA(:,2),LOP_R_Vertices_PCA(:,3),'FaceAlpha',0.2,'EdgeColor','None');
% else
%     msg = 'Lop or Lo or Me';
%     error(msg)
% end
% i=14;
% WantToScatter=Target_Table.syn_loc_R_PCA{11,1};
% scatter3(WantToScatter(:,1),WantToScatter(:,2),WantToScatter(:,3)); hold off;
% grid on;
% axis("equal")
% xlabel( 'PC1', 'Interpreter', 'none' );
% ylabel( 'PC2', 'Interpreter', 'none' );
% zlabel( 'PC3', 'Interpreter', 'none' );
%%

maxLayer=[];
% Iterate through each row
for i = 1:size(synapseBox_z_R, 2)
    % Find the column index with the maximum value in the current row
    [~, max_col] = max(synapseBox_z_R(:,i));
    maxLayer=[maxLayer max_col];
end
LayerIdx=[];
for i=1:1:size(synapseBox_z_R,1)
    idx=find(maxLayer==i);
    %%% sorting
    Values=synapseBox_z_R(i,idx);
    [B,ValuesIndex] = sort(Values,'descend');
    LayerIdx=[LayerIdx   idx(ValuesIndex)];
end
% figure(2); set(gcf,'Color','w')
% imagesc(synapseBox_z_R(:,LayerIdx))
% set(gca,'YTick',0.5:1:size(synapseBox_z_R,1)+0.5,'YTickLabel',edges,'XTick',1:1:size(Target_Neuron,1),'XTickLabel',Target_Neuron(LayerIdx),'TickDir','out','Box','off','YDir','normal')
%
% xlabel( 'Neurons');
% ylabel( 'Innvervaiton depth (μm)');
% colormap(customColormapOut)
% % %% density plot 그려보자
%
% densityBox=synapseBox_z(:,LayerIdx);
% figure(2); set(gcf,'Color','w')
% imagesc(densityBox)
% set(gca,'YTick',0.5:1:15.5,'YTickLabel',-20:5:50,'XTick',1:1:size(WantToSee,1),'XTickLabel',WantToSee(LayerIdx),'TickDir','out','Box','off')
%
% for i=1:1:size(densityBox,2)
%     plot(densityBox(:,i),'Color','#808080'); hold on;
%     plot(densityBox(:,i)); hold on;
%
% end
% plot(sum(densityBox,2)/max(sum(densityBox,2)),'Color','r','LineWidth',2); hold on;


%% 분산 그래프 그리기
figure(7); set(gcf,'Color','w'); hold on;

for i=1:1:size(Target_Table,1)
    pos=Target_Table.syn_loc_R_PCA{i};

    % --- 입력: pos (n x 3), 각 행은 시냅스의 3D 위치 ---
    n = size(pos, 1);
    numBoots = 1000;  % 샘플 수 기반 반복 수 자동 설정
    bootMeans = zeros(numBoots, 3);  % X, Y, Z 평균 저장

    % --- 부트스트랩 반복 (균등 샘플링) ---
    for iBoot = 1:numBoots
        idx = randsample(n, n, true);  % 균등 확률로 복원추출
        bootSample = pos(idx, :);
        bootMeans(iBoot, :) = mean(bootSample, 1);  % x, y, z 평균
    end

    % --- 신뢰구간 계산 ---
    CI_x = prctile(bootMeans(:,1), [2.5, 97.5]);
    CI_y = prctile(bootMeans(:,2), [2.5, 97.5]);
    CI_z = prctile(bootMeans(:,3), [2.5, 97.5]);

    mean_values = mean(bootMeans);  % 최종 평균 위치

    % --- 시각화 ---
    line_color = [0.8500 0.3250 0.0980]	;  % 파란색 계열

    % X축 신뢰구간 선
    line([CI_x(1), CI_x(2)], ...
        [mean_values(2), mean_values(2)], ...
        [mean_values(3), mean_values(3)], ...
        'Color', line_color, 'LineWidth', 1.5);

    % Y축 신뢰구간 선
    line([mean_values(1), mean_values(1)], ...
        [CI_y(1), CI_y(2)], ...
        [mean_values(3), mean_values(3)], ...
        'Color', line_color, 'LineWidth', 1.5);

    % Z축 신뢰구간 선
    line([mean_values(1), mean_values(1)], ...
        [mean_values(2), mean_values(2)], ...
        [CI_z(1), CI_z(2)], ...
        'Color', line_color, 'LineWidth', 1.5);

    % % 평균 위치에 마커 표시
    % plot3(mean_values(1), mean_values(2), mean_values(3), 'o', ...
    %     'MarkerSize', 6, 'MarkerFaceColor', line_color, 'MarkerEdgeColor', 'k');

    % 뉴런 이름 표시 (i는 바깥 루프에서 정의되어 있어야 함)
    text(mean_values(1), mean_values(2), mean_values(3), ...
        Target_Table.cell_type{i}, ...
        'FontSize', 12, 'HorizontalAlignment', 'center');

end

% 그래프 설정
xlabel( 'PC1', 'Interpreter', 'none' );
ylabel( 'PC2', 'Interpreter', 'none' );
zlabel( 'PC3', 'Interpreter', 'none' );
% title('3D Data Distribution with Mean and Variance');
grid off;
hold off;
% axis equal
% set(gca,'TickDir','out','ZDir','reverse')
set(gca,'TickDir','out')

% view([0, 0]);  % xZ 평면에서 보기 (X축을 기준으로 보는 시점 13 xz
view([90, 0]);  % YZ 평면에서 보기 (X축을 기준으로 보는 시점 23 yz
