% 데이터 로드
clear all; close all; clc;

opt = detectImportOptions('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\Figure5_BL\Mi1_columns_sungsooKim.csv');
opt = setvartype(opt,'mi1_root_id','int64');
Mi1_Columns = readtable('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\Figure5_BL\Mi1_columns_sungsooKim.csv',opt);
Mi1_Columns = rmmissing(Mi1_Columns);

Me_R_Faces = readmatrix('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\Neuropil\Me_R_faces.csv') + 1;
Me_R_Vertices = readmatrix('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\Neuropil\Me_R_vertices.csv');
Me_L_Faces = readmatrix('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\Neuropil\Me_L_faces.csv') + 1;
Me_L_Vertices = readmatrix('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\Neuropil\Me_L_vertices.csv');

AMe_R_Faces = readmatrix('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\Neuropil\AMe_R_faces.csv') + 1;
AMe_R_Vertices = readmatrix('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\Neuropil\AMe_R_vertices.csv');
AMe_L_Faces = readmatrix('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\Neuropil\AMe_L_faces.csv') + 1;
AMe_L_Vertices = readmatrix('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\Neuropil\AMe_L_vertices.csv');

opt = detectImportOptions('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\fafb_v783_princeton_synapse_table.csv');
opt = setvartype(opt,'pre_root_id_720575940','int64');
opt = setvartype(opt,'post_root_id_720575940','int64');
FAFB_synapse_coordinates = readtable('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\fafb_v783_princeton_synapse_table.csv',opt);
FAFB_synapse_coordinates.Properties.VariableNames(11) = "pre_root_id";
FAFB_synapse_coordinates.Properties.VariableNames(12) = "post_root_id";
FAFB_synapse_coordinates.pre_root_id=FAFB_synapse_coordinates.pre_root_id+int64(720575940000000000);
FAFB_synapse_coordinates.post_root_id=FAFB_synapse_coordinates.post_root_id+int64(720575940000000000);
%% MI1 neurons 복셀들 넣어놓기

swc_location='C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\sk_lod1_783_healed\';
for i=1:1:size(Mi1_Columns,1)
    swc_names=[swc_location num2str(Mi1_Columns.mi1_root_id(i)) '.swc'];
    swc_position=SWC_XYZ(swc_names);
    Mi1_Columns.voxel_position{i}=swc_position;

    currentVoxels=table2array(swc_position);

    currentVoxels = currentVoxels((intriangulation(Me_R_Vertices,Me_R_Faces,currentVoxels)...
        |intriangulation(AMe_R_Vertices,AMe_R_Faces,currentVoxels)...
        |intriangulation(Me_L_Vertices,Me_L_Faces,currentVoxels)...
        |intriangulation(AMe_L_Vertices,AMe_L_Faces,currentVoxels)),:);
    
    Mi1_Columns.voxel_position_Medulla{i}=table(currentVoxels(:,1),currentVoxels(:,2),currentVoxels(:,3),'VariableNames',{'x','y','z'});
end

%%
Me_R_Faces = readmatrix('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\Neuropil\Me_R_faces.csv');
Me_R_Vertices=readmatrix('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\Neuropil\Me_R_vertices.csv');
Me_R_Faces=Me_R_Faces+1;

Me_L_Faces = readmatrix('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\Neuropil\Me_L_faces.csv');
Me_L_Vertices=readmatrix('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\Neuropil\Me_L_vertices.csv');
Me_L_Faces=Me_L_Faces+1;
Mi1_Columns.Properties.VariableNames(1) = "root_id";

%% T4a voxel들 중 Lobula에 있는 것만..
%% synapse 위치
zeroIdx=[];
for i=1:1:size(Mi1_Columns,1)
    current_root_id=Mi1_Columns.root_id(i);

    out_syn_idx=find(ismember(FAFB_synapse_coordinates.pre_root_id,current_root_id));
    out_syn_loc=FAFB_synapse_coordinates(out_syn_idx,1:3);

    if isempty(out_syn_loc)
        zeroIdx=[zeroIdx i];
        continue;
    end
    out_syn_temp=table2array(out_syn_loc);

    out_syn_loc=out_syn_loc((intriangulation(Me_R_Vertices,Me_R_Faces,out_syn_temp)|intriangulation(Me_L_Vertices,Me_L_Faces,out_syn_temp)),:); % LO LR 바뀐듯
    if isempty(out_syn_loc)
        zeroIdx=[zeroIdx i];
        continue;
    end

    Mi1_Columns.out_syn_loc{i}=out_syn_loc;
end
Mi1_Columns(zeroIdx,:)=[];

%%
out_syn_loc_R=[];
out_syn_loc_L=[];
for i=1:1:size(Mi1_Columns,1)
    if isempty(Mi1_Columns.out_syn_loc{i})
        continue;
    elseif strcmp(Mi1_Columns.hemisphere{i},'R')
        out_syn_loc_R=[out_syn_loc_R;table2array(Mi1_Columns.out_syn_loc{i})];
    elseif strcmp(Mi1_Columns.hemisphere{i},'L')
        out_syn_loc_L=[out_syn_loc_L;table2array(Mi1_Columns.out_syn_loc{i})];
    end
end
%%
figure(1);set(gcf,'Color','w');hold on;
scatter3(out_syn_loc_R(:,1),out_syn_loc_R(:,2),out_syn_loc_R(:,3))
out_syn_loc_R_pointCloud = pointCloud(out_syn_loc_R);
out_syn_loc_R_pointCloud_Denoise = pcdenoise(out_syn_loc_R_pointCloud,"NumNeighbors",20,"Threshold",1);
out_syn_loc_R_pointCloud_Denoise=out_syn_loc_R_pointCloud_Denoise.Location;

scatter3(out_syn_loc_R_pointCloud_Denoise(:,1),out_syn_loc_R_pointCloud_Denoise(:,2),out_syn_loc_R_pointCloud_Denoise(:,3))
%%
figure(2);set(gcf,'Color','w');hold on;
scatter3(out_syn_loc_L(:,1),out_syn_loc_L(:,2),out_syn_loc_L(:,3))
out_syn_loc_L_pointCloud = pointCloud(out_syn_loc_L);
out_syn_loc_L_pointCloud_Denoise = pcdenoise(out_syn_loc_L_pointCloud,"NumNeighbors",20,"Threshold",1);
out_syn_loc_L_pointCloud_Denoise=out_syn_loc_L_pointCloud_Denoise.Location;

scatter3(out_syn_loc_L_pointCloud_Denoise(:,1),out_syn_loc_L_pointCloud_Denoise(:,2),out_syn_loc_L_pointCloud_Denoise(:,3))

%%
zeroIdx=[];
for i=1:1:size(Mi1_Columns,1)

    if isempty(Mi1_Columns.out_syn_loc{i})
        continue
    end
    current_out_syn_loc_Denoise=table2array(Mi1_Columns.out_syn_loc{i});
    idx=ismember(current_out_syn_loc_Denoise,out_syn_loc_R_pointCloud_Denoise,'rows')|ismember(current_out_syn_loc_Denoise,out_syn_loc_L_pointCloud_Denoise,'rows');
    if sum(idx)==0
        zeroIdx=[zeroIdx i];
        continue;
    end
    Mi1_Columns.out_syn_loc_denoise{i}=table(current_out_syn_loc_Denoise(idx,1),current_out_syn_loc_Denoise(idx,2),current_out_syn_loc_Denoise(idx,3),'VariableNames',{'x','y','z'});
end
Mi1_Columns(zeroIdx,:)=[];