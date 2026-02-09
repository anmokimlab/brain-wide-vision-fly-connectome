clear all; close all; clc

opt = detectImportOptions('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\consolidated_cell_types.csv');
opt = setvartype(opt,'root_id','int64');
FAFBConsolidated_type = readtable('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\consolidated_cell_types.csv',opt);


opt = detectImportOptions('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\connections_no_threshold.csv');
opt = setvartype(opt,'pre_root_id','int64');
opt = setvartype(opt,'post_root_id','int64');
FAFBConnections = readtable('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\connections_no_threshold.csv',opt);


opt = detectImportOptions('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\fafb_v783_princeton_synapse_table.csv');
opt = setvartype(opt,'pre_root_id_720575940','int64');
opt = setvartype(opt,'post_root_id_720575940','int64');
FAFB_synapse_coordinates = readtable('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\fafb_v783_princeton_synapse_table.csv',opt);
FAFB_synapse_coordinates.Properties.VariableNames(11) = "pre_root_id";
FAFB_synapse_coordinates.Properties.VariableNames(12) = "post_root_id";
FAFB_synapse_coordinates.pre_root_id=FAFB_synapse_coordinates.pre_root_id+int64(720575940000000000);
FAFB_synapse_coordinates.post_root_id=FAFB_synapse_coordinates.post_root_id+int64(720575940000000000);

%%
T4a_Columns=table(FAFBConsolidated_type.root_id(strcmp(FAFBConsolidated_type.primary_type,'T4a'),:),'VariableNames',{'root_id'});
load('Mi1_voxels.mat')

%%
zeroIdx=[];

for i=1:1:size(T4a_Columns,1)
    current_T4a=T4a_Columns.root_id(i);
    upstreamConnections=FAFBConnections(FAFBConnections.post_root_id==current_T4a,:);
    upstreamConnectionsMi1=upstreamConnections(ismember(upstreamConnections.pre_root_id,Mi1_Columns.root_id),:);
    upstreamMi1=unique(upstreamConnectionsMi1.pre_root_id);
    upstreamMi1=table(upstreamMi1,'VariableNames',{'root_id'});
    if isempty(upstreamMi1)
        zeroIdx=[zeroIdx i];
        continue;

    end
    for j=1:1:size(upstreamMi1,1)
        idx_connection=upstreamConnectionsMi1.pre_root_id==upstreamMi1.root_id(j);
        upstreamMi1.syn_count(j)=sum(upstreamConnectionsMi1.syn_count(idx_connection));

        idx_column=Mi1_Columns.root_id==upstreamMi1.root_id(j);
        upstreamMi1.theta(j)=Mi1_Columns.theta(idx_column);
        upstreamMi1.phi(j)=Mi1_Columns.phi(idx_column);
        upstreamMi1.hemisphere(j)=Mi1_Columns.hemisphere(idx_column);
    end
    T4a_Columns.upstreamMi1{i}=upstreamMi1;
    T4a_Columns.hemisphere(i)=upstreamMi1.hemisphere(1);
    T4a_Columns.theta{i}=[upstreamMi1.theta upstreamMi1.syn_count/sum(upstreamMi1.syn_count)];
    T4a_Columns.phi{i}=[upstreamMi1.phi upstreamMi1.syn_count/sum(upstreamMi1.syn_count)];
end
T4a_Columns(zeroIdx,:)=[];

%%
zeroIdx=[];

swc_location='C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\sk_lod1_783_healed\';
for i=1:1:size(T4a_Columns,1)
    swc_names=[swc_location num2str(T4a_Columns.root_id(i)) '.swc'];
    [swc_position]=SWC_XYZ(swc_names);
    T4a_Columns.voxel_position{i}=swc_position;
    if isempty(swc_position)
        zeroIdx=[zeroIdx i];
    end
end
T4a_Columns(zeroIdx,:)=[];
%%
LoP_R_Faces = readmatrix('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\Neuropil\LoP_R_faces.csv');
LoP_R_Vertices=readmatrix('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\Neuropil\LoP_R_vertices.csv');
LoP_R_Faces=LoP_R_Faces+1;

LoP_L_Faces = readmatrix('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\Neuropil\LoP_L_faces.csv');
LoP_L_Vertices=readmatrix('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\Neuropil\LoP_L_vertices.csv');
LoP_L_Faces=LoP_L_Faces+1;

%% T4a voxel들 중 Lobula에 있는 것만..
zeroIdx=[];
for i=1:1:size(T4a_Columns,1)
    currentVoxels=T4a_Columns.voxel_position{i};
    currentVoxels_array=table2array(currentVoxels);

    VoxelsInLobulaPlate=currentVoxels(intriangulation(LoP_L_Vertices,LoP_L_Faces,currentVoxels_array)|intriangulation(LoP_R_Vertices,LoP_R_Faces,currentVoxels_array),:); % LO LR 바뀐듯
    T4a_Columns.voxel_position_LobulaPlate{i}=VoxelsInLobulaPlate;
    if isempty(VoxelsInLobulaPlate)
        zeroIdx=[zeroIdx i];
    end
end
T4a_Columns(zeroIdx,:)=[];

%% synapse 위치
zeroIdx=[];
for i=1:1:size(T4a_Columns,1)
    current_root_id=T4a_Columns.root_id(i);

    out_syn_idx=find(ismember(FAFB_synapse_coordinates.pre_root_id,current_root_id));
    out_syn_loc=FAFB_synapse_coordinates(out_syn_idx,1:3);

    if isempty(out_syn_loc)
        zeroIdx=[zeroIdx i];
        continue;
    end
    out_syn_temp=table2array(out_syn_loc);

    out_syn_loc=out_syn_loc((intriangulation(LoP_R_Vertices,LoP_R_Faces,out_syn_temp)|intriangulation(LoP_L_Vertices,LoP_L_Faces,out_syn_temp)),:); % LO LR 바뀐듯
    if isempty(out_syn_loc)
        zeroIdx=[zeroIdx i];
        continue;
    end

    T4a_Columns.out_syn_loc{i}=out_syn_loc;
end
T4a_Columns(zeroIdx,:)=[];

%%
out_syn_loc_R=[];
out_syn_loc_L=[];
for i=1:1:size(T4a_Columns,1)
    if isempty(T4a_Columns.out_syn_loc{i})
        continue;
    elseif strcmp(T4a_Columns.hemisphere{i},'R')
        out_syn_loc_R=[out_syn_loc_R;table2array(T4a_Columns.out_syn_loc{i})];
    elseif strcmp(T4a_Columns.hemisphere{i},'L')
        out_syn_loc_L=[out_syn_loc_L;table2array(T4a_Columns.out_syn_loc{i})];
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
for i=1:1:size(T4a_Columns,1)

    if isempty(T4a_Columns.out_syn_loc{i})
        continue
    end
    current_out_syn_loc_Denoise=table2array(T4a_Columns.out_syn_loc{i});
    idx=ismember(current_out_syn_loc_Denoise,out_syn_loc_R_pointCloud_Denoise,'rows')|ismember(current_out_syn_loc_Denoise,out_syn_loc_L_pointCloud_Denoise,'rows');
    if sum(idx)==0
        zeroIdx=[zeroIdx i];
        continue;
    end
    T4a_Columns.out_syn_loc_denoise{i}=table(current_out_syn_loc_Denoise(idx,1),current_out_syn_loc_Denoise(idx,2),current_out_syn_loc_Denoise(idx,3),'VariableNames',{'x','y','z'});
end
T4a_Columns(zeroIdx,:)=[];
%% ===================== 설정 =====================
N_ANC = 20;        % 조상 단계를 최대 100개
zeroIdx = [];       % out_syn 없음 행 모으기
swc_location='C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\sk_lod1_783_healed\';

%%% ===================== 메인 루프 =====================
for i = 1:height(T4a_Columns)
    %%% 0) 출력 시냅스 없으면 스킵 (나중에 행 삭제)

    swcPath=[swc_location num2str(T4a_Columns.root_id(i)) '.swc'];

    % 1) 입력 데이터 꺼내기
    A = table2array(T4a_Columns.voxel_position_LobulaPlate{i,1});  % [M x 3] 복셀 좌표
    B = table2array(T4a_Columns.out_syn_loc_denoise{i,1});             % [K x 3] 시냅스 좌표

    % 방어: 비어있으면 빈 저장 후 continue
    if isempty(A) || isempty(B) || ~isfile(swcPath)
        T4a_Columns.ancestors_outsyn{i} = table([],[],[], 'VariableNames',{'x','y','z'});
        warning('Row %d: A/B empty or SWC missing -> empty result saved.', i);
        zeroIdx=[zeroIdx i];
        continue
    end

    % 2) SWC 로드 (id,type,x,y,z,r,parent)
    swcData =SWC_XYZ_full(swcPath);

    % 3) 각 시냅스에서 "가장 가까운 복셀" 인덱스(시드 복셀) 선택
    seedVoxelIdx = nearestVoxelSeeds(A, B);     % [K x 1]
    seedVoxels   = unique(A(seedVoxelIdx, :), 'rows', 'stable');  % 좌표 기준 중복 제거

    %%% 4) 시드 복셀 -> 가장 가까운 SWC 노드 매핑 (툴박스 없이 pdist2 사용 권장)
    swcXYZ = swcData{:, {'x','y','z'}};                 % [Ms x 3]
    if isempty(swcXYZ) || isempty(seedVoxels)
        T4a_Columns.ancestors_outsyn{i} = table([],[],[], 'VariableNames',{'x','y','z'});
        warning('Row %d: empty swcXYZ or seedVoxels -> empty result saved.', i);
        continue
    end

    % knnsearch 대신 pdist2로 호환성↑ (Statistics Toolbox 없어도 됨)
    Dsv = pdist2(seedVoxels, swcXYZ);                   % [S x Ms]
    [~, idxSeedNode] = min(Dsv, [], 2);                 % [S x 1], 각 시드 복셀의 최근접 SWC row 인덱스

    %%% 5) id -> row 매핑 준비 (parent는 "id"를 가리킴)
    if ismember('id', swcData.Properties.VariableNames)
        ids = swcData.id;                               % 진짜 SWC id
        [id2row, hasMap] = build_id2row(ids);
    else
        % id가 없으면 row==id로 가정(안전성↓). parent가 실제 row 번호여야만 동작.
        ids = (1:height(swcData))';
        id2row = []; hasMap = false;
    end
    parentId = swcData.parent_node;                           % parent는 "id" 값이어야 함 (루트 -1)

    %%% 6) 각 시드 노드에서 조상 N 단계 수집 후 유니크
    allAncRows = [];
    for s = 1:numel(idxSeedNode)
        ancRows = trace_ancestors(idxSeedNode(s), parentId, ids, id2row, hasMap, N_ANC);
        if ~isempty(ancRows)
            allAncRows = [allAncRows; ancRows(:)]; %#ok<AGROW>
        end
    end
    allAncRows = unique(allAncRows, 'stable');

    %% 7) 좌표 추출 후 테이블 저장
    V = swcData{allAncRows, {'x','y','z'}};
    V = unique(V, 'rows', 'stable');                     % 좌표 유니크(안전)
    T4a_Columns.ancestors_outsyn{i} = array2table(V, 'VariableNames', {'x','y','z'});

end


%% ===================== 보조 함수들 =====================
function [swc_position]=SWC_XYZ_full(swcFilePath)
% Function to visualize a neuron from an SWC file in MATLAB
% Input:
%   swcFilePath - string, path to the SWC file

% Read the SWC file
data = readSWCFile(swcFilePath);
% end_idx=data(:,2)==6;
x = data(:, 3);             % X coordinates
y = data(:, 4);             % Y coordinates
z = data(:, 5);             % Z coordinates
parent_node=data(:,7);
% end_point= data(end_idx,[3,4,5]);
% end_point=table(end_point(:,1),end_point(:,2),end_point(:,3),'VariableNames',{'x','y','z'});
swc_position=table(x,y,z,parent_node);

end

function data = readSWCFile(swcFilePath)
% Reads an SWC file and returns the data as a matrix
% Input:
%   swcFilePath - string, path to the SWC file
% Output:
%   data - matrix, parsed SWC data

% Open the file
fid = fopen(swcFilePath, 'r');
if fid == -1
    error('Unable to open SWC file: %s', swcFilePath);
end

% Read the file, skipping comments
data = [];
while ~feof(fid)
    line = fgetl(fid);
    if isempty(line) || startsWith(line, '#')
        continue;
    end
    data = [data; sscanf(line, '%f')']; %#ok<AGROW>
end

% Close the file
fclose(fid);
end

function seedIdx = nearestVoxelSeeds(A, B)
% 각 시냅스 B(k,:)에 가장 가까운 복셀 A의 인덱스 반환 (중복 허용)
if isempty(B) || isempty(A)
    seedIdx = zeros(0,1);
    return
end
D = pdist2(B, A);            % [K x M]
[~, seedIdx] = min(D, [], 2);
end

function [id2row, hasMap] = build_id2row(ids)
% SWC id가 1..M 연속이면 맵 불필요, 아니면 Map 구축
    M = numel(ids);
    if all(ismember(1:M, ids)) && all(ismember(ids, 1:M))
        id2row = []; hasMap = false;      % id==row 간주 가능
    else
        K = strcat("i", string(ids));     % 문자열 키
        V = num2cell(1:M);
        id2row = containers.Map(K, V);
        hasMap = true;
    end
end

function ancRows = trace_ancestors(seedRow, parentCol, ids, id2row, hasMap, N)
% seedRow: 시작 row 인덱스 (knnsearch/pdist2 결과)
% parentCol: 각 row의 parent "id" (루트는 -1)
% ids: 각 row의 "id"
% 반환: row 인덱스 벡터(시드 포함), 최대 N+1개
    M = numel(parentCol);
    if seedRow < 1 || seedRow > M
        ancRows = [];
        return
    end
    ancRows = seedRow;
    curRow = seedRow;
    for t = 1:N
        pId = parentCol(curRow);          % parent는 "id"
        if pId < 0, break; end            % -1이면 루트
        if hasMap
            key = "i" + string(pId);
            if ~isKey(id2row, key), break; end
            pRow = id2row(key);
        else
            pRow = pId;                   % id==row 가정 (fallback)
            if pRow < 1 || pRow > M || ids(pRow) ~= pId
                break;                    % 불일치 시 중단
            end
        end
        ancRows(end+1) = pRow; %#ok<AGROW>
        curRow = pRow;
    end
end