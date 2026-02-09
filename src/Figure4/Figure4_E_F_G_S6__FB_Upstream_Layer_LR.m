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


opt = detectImportOptions('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\connections_no_threshold.csv');
opt = setvartype(opt,'pre_root_id','int64');
opt = setvartype(opt,'post_root_id','int64');
FAFBConnections = readtable('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\connections_no_threshold.csv',opt);

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
% colorbar
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
Target_What='dendrite';    %%%%%%%% 이거 주의...!! upstream이어서 무조건 dend
Target_Where='Lo';

%% Ref and Target root ids;
Ref_Neuron={'LT1a'}; %% Dm6 dendrite LT1a dendrite T5a axon
Ref_What='dendrite';
Ref_Where='Lo';



%%
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

reverse='False';
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
view([-41.264686231577031,27.673914691032877]) %%%%%%lop

ax1 = gca;  % 현재 figure(1)의 축 가져오기
x_limits = xlim(ax1);
y_limits = ylim(ax1);
z_limits = zlim(ax1);

figure(2);set(gcf,'Color','w');
hold on;
[X,Y] = meshgrid(min(transformedData_L(:,1)):5e-7:max(transformedData_L(:,1)),min(transformedData_L(:,2)):5e-7:max(transformedData_L(:,2)));
Z = paraboloid_surface_L(X,Y);
surf(X,Y,Z,'EdgeColor','none')
scatter3(transformedData_L(:,1),transformedData_L(:,2),transformedData_L(:,3),20,[0.8500 0.3250 0.0980],'filled');
% scatter3(transformedData_R(:,1),transformedData_R(:,2),transformedData_R(:,3),20,[0 0.4470 0.7410]	,'filled'); hold on;

grid on;

FaceAlhpaValue=0.7;
if strcmpi(Ref_Where,'Me')
    ME_L_Vertices_PCA=ME_L_Vertices;
    ME_L_Vertices_PCA(:,1)=ME_L_Vertices_PCA(:,1)-Ref_x_mean_L;
    ME_L_Vertices_PCA(:,2)=ME_L_Vertices_PCA(:,2)-Ref_y_mean_L;
    ME_L_Vertices_PCA(:,3)=ME_L_Vertices_PCA(:,3)-Ref_z_mean_L;
    ME_L_Vertices_PCA=ME_L_Vertices_PCA*coeffs_L;
    trisurf(ME_L_Faces,ME_L_Vertices_PCA(:,1),ME_L_Vertices_PCA(:,2),ME_L_Vertices_PCA(:,3),'FaceAlpha',0.2,'EdgeColor','None');
elseif strcmpi(Ref_Where,'Lo')
    LO_L_Vertices_PCA=LO_L_Vertices;
    LO_L_Vertices_PCA(:,1)=LO_L_Vertices_PCA(:,1)-Ref_x_mean_L;
    LO_L_Vertices_PCA(:,2)=LO_L_Vertices_PCA(:,2)-Ref_y_mean_L;
    LO_L_Vertices_PCA(:,3)=LO_L_Vertices_PCA(:,3)-Ref_z_mean_L;
    LO_L_Vertices_PCA=LO_L_Vertices_PCA*coeffs_L;
    trisurf(LO_L_Faces,LO_L_Vertices_PCA(:,1),LO_L_Vertices_PCA(:,2),LO_L_Vertices_PCA(:,3),'FaceAlpha',0.2,'EdgeColor','None');
elseif strcmpi(Ref_Where,'Lop')
    LOP_L_Vertices_PCA=LOP_L_Vertices;
    LOP_L_Vertices_PCA(:,1)=LOP_L_Vertices_PCA(:,1)-Ref_x_mean_L;
    LOP_L_Vertices_PCA(:,2)=LOP_L_Vertices_PCA(:,2)-Ref_y_mean_L;
    LOP_L_Vertices_PCA(:,3)=LOP_L_Vertices_PCA(:,3)-Ref_z_mean_L;
    LOP_L_Vertices_PCA=LOP_L_Vertices_PCA*coeffs_L;
    trisurf(LOP_L_Faces,LOP_L_Vertices_PCA(:,1),LOP_L_Vertices_PCA(:,2),LOP_L_Vertices_PCA(:,3),'FaceAlpha',0.2,'EdgeColor','None');
else
    msg = 'Lop or Lo or Me';
    error(msg)
end

axis("equal")
xlabel( 'PC1', 'Interpreter', 'none' );
ylabel( 'PC2', 'Interpreter', 'none' );
zlabel( 'PC3', 'Interpreter', 'none' );
view([-41.264686231577031,27.673914691032877]) %%%%%%lop

ax1 = gca;  % 현재 figure(1)의 축 가져오기
x_limits = xlim(ax1);
y_limits = ylim(ax1);
z_limits = zlim(ax1);
%% Target Synapse위치
edges = -7e-5:1e-6:5e-5;

Target_Upstream_Thr=0;
for i=1:1:size(Target_Table,1)
    Target_root_ids=Target_Table.root_ids{i};
    [upstreamNeurons]=seeConnection_root_id_NoOptic(Target_root_ids,FAFBConnections,FAFBConsolidated_type);

    syn_zDist_from_Parabola_one_FB_R=zeros(size(edges,2)-1,size(upstreamNeurons,1));
    syn_zDist_from_Parabola_one_FB_L=zeros(size(edges,2)-1,size(upstreamNeurons,1));
    upstreamNeuronsCell=num2cell(upstreamNeurons);
    for j=1:1:size(upstreamNeurons,1)
        current_upstream_root_ids=upstreamNeurons(j,1);

        upstream_syn_loc=[];

        if strcmpi(Target_What,'axon')
            Target_Syn_pre_root_idx=find(ismember(FAFB_synapse_coordinates.pre_root_id,current_upstream_root_ids));
            upstream_syn_loc=FAFB_synapse_coordinates(Target_Syn_pre_root_idx,1:3);
            

        elseif strcmpi(Target_What,'dendrite')
            Target_Syn_post_root_idx=find(ismember(FAFB_synapse_coordinates.post_root_id,current_upstream_root_ids));
            upstream_syn_loc=FAFB_synapse_coordinates(Target_Syn_post_root_idx,1:3);
        else
            msg = 'axon or dendrite';
            error(msg)
        end
        if isempty(upstream_syn_loc)
            syn_zDist_from_Parabola_one_FB_L(:,j)=nan(1,size(edges,2)-1);
            syn_zDist_from_Parabola_one_FB_R(:,j)=nan(1,size(edges,2)-1);
            continue;
        end

        upstream_syn_loc=upstream_syn_loc{:,:};
        if strcmpi(Target_Where,'Lop')
            upstream_syn_loc_in_R=upstream_syn_loc(intriangulation(LOP_R_Vertices,LOP_R_Faces,upstream_syn_loc),:);
            upstream_syn_loc_in_L=upstream_syn_loc(intriangulation(LOP_L_Vertices,LOP_L_Faces,upstream_syn_loc),:);
        elseif strcmpi(Target_Where,'Lo')
            upstream_syn_loc_in_R=upstream_syn_loc(intriangulation(LO_R_Vertices,LO_R_Faces,upstream_syn_loc),:);
            upstream_syn_loc_in_L=upstream_syn_loc(intriangulation(LO_L_Vertices,LO_L_Faces,upstream_syn_loc),:);
        elseif strcmpi(Target_Where,'Me')
            upstream_syn_loc_in_R=upstream_syn_loc(intriangulation(ME_R_Vertices,ME_R_Faces,upstream_syn_loc),:);
            upstream_syn_loc_in_L=upstream_syn_loc(intriangulation(ME_L_Vertices,ME_L_Faces,upstream_syn_loc),:);
        else
            msg = 'Lop or Lo or Me';
            error(msg)
        end

        counts_R=zeros(1,size(edges,2)-1);
        counts_L=zeros(1,size(edges,2)-1);

        if ~isempty(upstream_syn_loc_in_R)
            upstream_syn_loc_in_R_pointCloud = pointCloud(upstream_syn_loc_in_R);
            upstream_syn_loc_in_R_pointCloud_Denoise = pcdenoise(upstream_syn_loc_in_R_pointCloud,"NumNeighbors",20,"Threshold",2.5);
            upstream_syn_loc_in_R_Denoise=upstream_syn_loc_in_R_pointCloud_Denoise.Location;

            % upstream_syn_loc_in_R_shift=upstream_syn_loc_in_R_Denoise; %%%%%%%%%%%%%%%
            upstream_syn_loc_in_R_shift=upstream_syn_loc_in_R; %%%%%%%%%%%%%%%

            %%% center shift
            upstream_syn_loc_in_R_shift(:,1)=upstream_syn_loc_in_R_shift(:,1)-Ref_x_mean_R;
            upstream_syn_loc_in_R_shift(:,2)=upstream_syn_loc_in_R_shift(:,2)-Ref_y_mean_R;
            upstream_syn_loc_in_R_shift(:,3)=upstream_syn_loc_in_R_shift(:,3)-Ref_z_mean_R;
            %%% pca axis
            upstream_syn_loc_in_R_shift_PCA=upstream_syn_loc_in_R_shift*coeffs_R;

            z_hat_R=paraboloid_surface_R(upstream_syn_loc_in_R_shift_PCA(:,1),upstream_syn_loc_in_R_shift_PCA(:,2));
            [counts_R, ~] = histcounts(upstream_syn_loc_in_R_shift_PCA(:,3)-z_hat_R, 'BinEdges', edges);
            upstreamNeuronsCell{j,4}=upstream_syn_loc_in_R_shift_PCA;

        end

        if ~isempty(upstream_syn_loc_in_L)
            upstream_syn_loc_in_L_pointCloud = pointCloud(upstream_syn_loc_in_L);
            upstream_syn_loc_in_L_pointCloud_Denoise = pcdenoise(upstream_syn_loc_in_L_pointCloud,"NumNeighbors",20,"Threshold",2.5);
            upstream_syn_loc_in_L_Denoise=upstream_syn_loc_in_L_pointCloud_Denoise.Location;

            % upstream_syn_loc_in_L_shift=upstream_syn_loc_in_L_Denoise; %%%%%%%%%%%%%%%%%%%%%%%%
            upstream_syn_loc_in_L_shift=upstream_syn_loc_in_L; %%%%%%%%%%%%%%%%%%%%%%%%

            %%% center shift
            upstream_syn_loc_in_L_shift(:,1)=upstream_syn_loc_in_L_shift(:,1)-Ref_x_mean_L;
            upstream_syn_loc_in_L_shift(:,2)=upstream_syn_loc_in_L_shift(:,2)-Ref_y_mean_L;
            upstream_syn_loc_in_L_shift(:,3)=upstream_syn_loc_in_L_shift(:,3)-Ref_z_mean_L;
            %%% pca axis
            upstream_syn_loc_in_L_shift_PCA=upstream_syn_loc_in_L_shift*coeffs_L;

            z_hat_L=paraboloid_surface_L(upstream_syn_loc_in_L_shift_PCA(:,1),upstream_syn_loc_in_L_shift_PCA(:,2));
            [counts_L, ~] = histcounts(upstream_syn_loc_in_L_shift_PCA(:,3)-z_hat_L, 'BinEdges', edges);
            % if strcmpi(Target_Where,'Me')
            %    [counts_L, ~] = histcounts(-upstream_syn_loc_in_L_shift_PCA(:,3)+z_hat_L, 'BinEdges', edges);
            % end
            upstreamNeuronsCell{j,3}=upstream_syn_loc_in_L_shift_PCA;
        end


        counts_sum=sum([counts_L counts_R]);

        counts_L=counts_L/counts_sum*double(upstreamNeurons(j,2))/size(upstream_syn_loc,1)*size(upstream_syn_loc_in_L,1);
        counts_R=counts_R/counts_sum*double(upstreamNeurons(j,2))/size(upstream_syn_loc,1)*size(upstream_syn_loc_in_R,1);

        % Upstream_Target{j,9}=syn_histcounts*Upstream_Target{j,2}/size(Upstream_Target{j,6},1)*size(Upstream_Target{j,7},1);
        % Upstream_Target{j,10}=syn_histcounts*Upstream_Target{j,2}/size(Upstream_Target{j,6},1)*(size(Upstream_Target{j,7},1)+size(Upstream_Target{j,8},1));
        % Upstream_Target{j,10}=syn_histcounts*웨이트/전체시냅스갯수*(왼쪽+오른쪽);

        syn_zDist_from_Parabola_one_FB_L(:,j)=counts_L';
        syn_zDist_from_Parabola_one_FB_R(:,j)=counts_R';
    end
    syn_zDist_from_Parabola_one_FB_All_L{i}=syn_zDist_from_Parabola_one_FB_L;
    syn_zDist_from_Parabola_one_FB_All_R{i}=syn_zDist_from_Parabola_one_FB_R;

    Target_Table.upstream_synapses{i}=upstreamNeuronsCell;

    Target_Table.syn_zDist_from_Parabola_L{i}=sum(syn_zDist_from_Parabola_one_FB_L,2,'omitnan')/max([sum(syn_zDist_from_Parabola_one_FB_L,2,'omitnan');sum(syn_zDist_from_Parabola_one_FB_R,2,'omitnan')]);
    Target_Table.syn_zDist_from_Parabola_R{i}=sum(syn_zDist_from_Parabola_one_FB_R,2,'omitnan')/max([sum(syn_zDist_from_Parabola_one_FB_L,2,'omitnan');sum(syn_zDist_from_Parabola_one_FB_R,2,'omitnan')]);
end


%%
% figure(1);set(gcf,'Color','w');
% temp_L=syn_zDist_from_Parabola_one_FB_All_L{i};
% imagesc(syn_zDist_from_Parabola_one_FB_All_L{i})
% figure(2);set(gcf,'Color','w');
% temp_R=syn_zDist_from_Parabola_one_FB_All_R{i};
% imagesc(syn_zDist_from_Parabola_one_FB_All_R{i})
%%
synapseBox_z_Target_L=zeros(size(edges,2)-1,size(Target_Table,1));

for i=1:1:size(Target_Table,1)
    if ~isempty(Target_Table.syn_zDist_from_Parabola_L{i})
        synapseBox_z_Target_L(:,i)=Target_Table.syn_zDist_from_Parabola_L{i};
    end
end
%%
figure(3); set(gcf,'Color','w')
imagesc(synapseBox_z_Target_L)
set(gca,'YTick',0.5:1:size(synapseBox_z_Target_L,1)+0.5,'YTickLabel',edges,'XTick',1:1:size(Target_Table,1),'XTickLabel',Target_Table.cell_type,'TickDir','out','Box','off','YDir','normal')
xlabel( 'Neurons');
ylabel( 'Innvervaiton depth (μm)');
colormap(customColormapIn_contra)
xtickangle(45)
title('L')
clim([0 1])

%%
synapseBox_z_Target_R=zeros(size(edges,2)-1,size(Target_Table,1));

for i=1:1:size(Target_Table,1)
    if ~isempty(Target_Table.syn_zDist_from_Parabola_R{i})
        synapseBox_z_Target_R(:,i)=Target_Table.syn_zDist_from_Parabola_R{i};
    end
end
%%
figure(4); set(gcf,'Color','w')
imagesc(synapseBox_z_Target_R)
set(gca,'YTick',0.5:1:size(synapseBox_z_Target_R,1)+0.5,'YTickLabel',edges,'XTick',1:1:size(Target_Table,1),'XTickLabel',Target_Table.cell_type,'TickDir','out','Box','off','YDir','normal')
xlabel( 'Neurons');
ylabel( 'Innvervaiton depth (μm)');
colormap(customColormapIn)
xtickangle(45)
title('R')
clim([0 1])

%% 분산 그래프 그리기 L
figure(5); set(gcf,'Color','w'); hold on;

for i=1:1:size(Target_Table,1)
    Upstream_syn=Target_Table.upstream_synapses{i,1} ;
    if size(Upstream_syn,2)==2
        continue;
    end
    pos = {};  % 결과 저장용
    w=[];
    for j = 1:size(Upstream_syn,1)
        current = Upstream_syn{j,3};
        if ~isempty(current)
            pos = [pos; current];  % 빈 셀이 아닌 경우만 이어붙이기
            w=[w; double(Upstream_syn{j,2})];
        end
    end
    if isempty(w)
        continue;
    end
    % ===== 1. 시냅스 위치 통합 및 가중치 확장 =====
    all_pos = [];  % 모든 시냅스 위치 (n x 3)
    all_w = [];    % 각 시냅스에 대응되는 weight (n x 1)

    for j = 1:length(pos)
        n_syn = size(pos{j}, 1);       % 해당 뉴런의 시냅스 개수
        all_pos = [all_pos; pos{j}];   % 위치 추가
        all_w = [all_w; repmat(w(j), n_syn, 1)];  % 같은 weight 반복 추가
    end

    p = all_w / sum(all_w);
    n = size(all_pos, 1);  % 총 시냅스 수
    numBoots = 1000;
    bootMeans = zeros(numBoots, 3);

    TotalSampled_Pos=[];
    for j = 1:numBoots
        idx = randsample(1:n, n, true, p);  % 확률 기반 샘플링
        sampled_pos = all_pos(idx, :);
        TotalSampled_Pos=[TotalSampled_Pos; sampled_pos];
        sampled_w = all_w(idx);
    
        % 가중 평균 계산
        bootMeans(j, 1) = mean(sampled_pos(:,1));  % x
        bootMeans(j, 2) = mean(sampled_pos(:,2));  % y
        bootMeans(j, 3) = mean(sampled_pos(:,3));  % z
    end
    z_hat_L=paraboloid_surface_L(TotalSampled_Pos(:,1),TotalSampled_Pos(:,2));
    Target_Table.syn_zDist_from_Parabola_L_Raw{i}=TotalSampled_Pos(:,3)-z_hat_L;
    if strcmpi(Target_Where,'Me')
        Target_Table.syn_zDist_from_Parabola_L_Raw{i}=-TotalSampled_Pos(:,3)+z_hat_L;
    end
    CI_x = prctile(bootMeans(:,1), [2.5, 97.5]);
    CI_y = prctile(bootMeans(:,2), [2.5, 97.5]);
    CI_z = prctile(bootMeans(:,3), [2.5, 97.5]);

    % mean 위치 (부트스트랩 평균 위치)
    mean_values = mean(bootMeans);  % 1x3 벡터 [mean_x, mean_y, mean_z]

    % CI_x, CI_y, CI_z는 부트스트랩에서 이미 구한 95% 신뢰구간

    % 색상 설정 (optional)
    line_color = [0.4660 0.6740 0.1880]	;

    % 신뢰구간 라인 (X축)
    line([CI_x(1), CI_x(2)], ...
        [mean_values(2), mean_values(2)], ...
        [mean_values(3), mean_values(3)], ...
        'Color', line_color, 'LineWidth', 1);

    % 신뢰구간 라인 (Y축)
    line([mean_values(1), mean_values(1)], ...
        [CI_y(1), CI_y(2)], ...
        [mean_values(3), mean_values(3)], ...
        'Color', line_color, 'LineWidth', 1);

    % 신뢰구간 라인 (Z축)
    line([mean_values(1), mean_values(1)], ...
        [mean_values(2), mean_values(2)], ...
        [CI_z(1), CI_z(2)], ...
        'Color', line_color, 'LineWidth', 1);

    % 뉴런 이름 텍스트 추가 (i는 루프 외부에서 지정되어 있어야 함)
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
view([0, 0]);  % xZ 평면에서 보기 (X축을 기준으로 보는 시점 13 xz
% view([90, 0]);  % YZ 평면에서 보기 (X축을 기준으로 보는 시점 23 yz

%% 분산 그래프 그리기 R
figure(6); set(gcf,'Color','w'); hold on;

for i=1:1:size(Target_Table,1)
    Upstream_syn=Target_Table.upstream_synapses{i,1} ;
    if size(Upstream_syn,2)<=3
        continue;
    end
    pos = {};  % 결과 저장용
    w=[];
    for j = 1:size(Upstream_syn,1)
        current = Upstream_syn{j,4};
        if ~isempty(current)
            pos = [pos; current];  % 빈 셀이 아닌 경우만 이어붙이기
            w=[w; double(Upstream_syn{j,2})];
        end
    end
    if isempty(w)
        continue;
    end
    % ===== 1. 시냅스 위치 통합 및 가중치 확장 =====
    all_pos = [];  % 모든 시냅스 위치 (n x 3)
    all_w = [];    % 각 시냅스에 대응되는 weight (n x 1)

    for j = 1:length(pos)
        n_syn = size(pos{j}, 1);       % 해당 뉴런의 시냅스 개수
        all_pos = [all_pos; pos{j}];   % 위치 추가
        all_w = [all_w; repmat(w(j), n_syn, 1)];  % 같은 weight 반복 추가
    end

    p = all_w / sum(all_w);
    n = size(all_pos, 1);  % 총 시냅스 수
    numBoots = 1000;
    bootMeans = zeros(numBoots, 3);
    TotalSampled_Pos=[];
    for j = 1:numBoots
        idx = randsample(1:n, n, true, p);  % 확률 기반 샘플링
        sampled_pos = all_pos(idx, :);
        TotalSampled_Pos=[TotalSampled_Pos; sampled_pos];
        sampled_w = all_w(idx);
    
        % 가중 평균 계산
        bootMeans(j, 1) = mean(sampled_pos(:,1));  % x
        bootMeans(j, 2) = mean(sampled_pos(:,2));  % y
        bootMeans(j, 3) = mean(sampled_pos(:,3));  % z
    end
    z_hat_R=paraboloid_surface_L(TotalSampled_Pos(:,1),TotalSampled_Pos(:,2));
    Target_Table.syn_zDist_from_Parabola_R_Raw{i}=TotalSampled_Pos(:,3)-z_hat_R;

    CI_x = prctile(bootMeans(:,1), [2.5, 97.5]);
    CI_y = prctile(bootMeans(:,2), [2.5, 97.5]);
    CI_z = prctile(bootMeans(:,3), [2.5, 97.5]);

    % mean 위치 (부트스트랩 평균 위치)
    mean_values = mean(bootMeans);  % 1x3 벡터 [mean_x, mean_y, mean_z]

    % CI_x, CI_y, CI_z는 부트스트랩에서 이미 구한 95% 신뢰구간

    % 색상 설정 (optional)
    line_color = [0 0.4470 0.7410]	;

    % 신뢰구간 라인 (X축)
    line([CI_x(1), CI_x(2)], ...
        [mean_values(2), mean_values(2)], ...
        [mean_values(3), mean_values(3)], ...
        'Color', line_color, 'LineWidth', 1);

    % 신뢰구간 라인 (Y축)
    line([mean_values(1), mean_values(1)], ...
        [CI_y(1), CI_y(2)], ...
        [mean_values(3), mean_values(3)], ...
        'Color', line_color, 'LineWidth', 1);

    % 신뢰구간 라인 (Z축)
    line([mean_values(1), mean_values(1)], ...
        [mean_values(2), mean_values(2)], ...
        [CI_z(1), CI_z(2)], ...
        'Color', line_color, 'LineWidth', 1);

    % 뉴런 이름 텍스트 추가 (i는 루프 외부에서 지정되어 있어야 함)
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
view([0, 0]);  % xZ 평면에서 보기 (X축을 기준으로 보는 시점 13 xz
% view([90, 0]);  % YZ 평면에서 보기 (X축을 기준으로 보는 시점 23 yz
