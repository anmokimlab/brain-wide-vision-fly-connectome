%%
clear all; close all; clc;

load('BD_Neuron_In_Out.mat')

opt = detectImportOptions('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\consolidated_cell_types.csv');
opt = setvartype(opt,'root_id','int64');
FAFB_consolidated_cell_types = readtable('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\consolidated_cell_types.csv',opt);

opt = detectImportOptions('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\classification.csv');
opt = setvartype(opt,'root_id','int64');
FAFB_classfication = readtable('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\classification.csv',opt);

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

%%
LC9_Total_in=Wantsee.InNeuronTypes{1};
LC9_Total_in = cell2table(LC9_Total_in(:,[1 2]), 'VariableNames', {'Type','Synapse'});

LC9_Total_out=Wantsee.OutNeuronTypes{1};
LC9_Total_out = cell2table(LC9_Total_out(:,[1 2]), 'VariableNames', {'Type','Synapse'});



LC9_in_optic=Wantsee.unique_InOptic_Types{1};
LC9_in_optic = cell2table(LC9_in_optic, 'VariableNames', {'Type','Synapse','Number'});
LC9_in_central=Wantsee.unique_InCentral_Types{1};
LC9_in_central = cell2table(LC9_in_central, 'VariableNames', {'Type','Synapse','Number'});
LC9_out_optic=Wantsee.unique_OutOptic_Types{1};
LC9_out_optic = cell2table(LC9_out_optic, 'VariableNames', {'Type','Synapse','Number'});
LC9_out_central=Wantsee.unique_OutCentral_Types{1};
LC9_out_central = cell2table(LC9_out_central, 'VariableNames', {'Type','Synapse','Number'});


%% In LC9
Indata = [sum(LC9_in_central.Synapse),sum(LC9_in_optic.Synapse)];
figure;set(gcf,'Color','w')
piechart(Indata,Names={'Central','Optic'})
title('LC9 Input');

in_Thr=1000;
LC9_Total_in_Thr=LC9_Total_in;
LC9_Total_in_Thr(LC9_Total_in_Thr.Synapse<in_Thr,:)=[];

for i=1:1:size(LC9_Total_in_Thr,1)
    current_types=LC9_Total_in_Thr.Type{i};
    current_types_root_ids=FAFB_consolidated_cell_types.root_id(strcmp(FAFB_consolidated_cell_types.primary_type,current_types));
    superclasses=FAFB_classfication.super_class(ismember(FAFB_classfication.root_id,current_types_root_ids));
    [u, ~, idx] = unique(superclasses);           % 고유 문자열 u와 index
    counts = histcounts(idx, 1:numel(u)+1); % 각 문자열 빈도
    [~, maxIdx] = max(counts);         % 최빈값 인덱스
    mostCommon_superclasses = superclasses{maxIdx};            % 최빈 문자열
    LC9_Total_in_Thr.superclass{i}=mostCommon_superclasses;
end

LC9_Total_in_Thr = sortrows(LC9_Total_in_Thr,"superclass","ascend");
LC9_Total_in_Thr = sortrows(LC9_Total_in_Thr,"Type","ascend");
LC9_Total_in_Thr = sortrows(LC9_Total_in_Thr,"superclass","ascend");
% p=piechart(LC9_in_central_Thr,"Synapse","Type",'LabelStyle','name');
% 
% %% Pie chart with custom colors based on superclass

% 1. LC9_in_central_Thr 테이블의 순서에 맞는 색상 행렬(pieColors) 생성
% 테이블의 행 개수만큼 3열(RGB)짜리 행렬을 0으로 초기화
numSlices = height(LC9_Total_in_Thr);
pieColors = zeros(numSlices, 3);

% for 루프를 사용해 테이블의 각 행을 순회
for i = 1:1:numSlices
    % 현재 행의 superclass를 가져옴
    current_superclass = LC9_Total_in_Thr.superclass{i};
    
    % uniqueSuperclass 목록에서 현재 superclass와 일치하는 위치(인덱스)를 찾음
    % matches 함수는 strcmp와 유사하게 문자열을 비교하여 논리형 배열을 반환
    color_index = matches(uniqueSuperclass, current_superclass);
    
    % 찾은 인덱스를 사용해 미리 정의된 colors 행렬에서 해당 색상(RGB 값)을 가져옴
    % pieColors의 현재 행에 색상 값을 할당
    pieColors(i, :) = colors(color_index, :);
end
figure;set(gcf,'Color','w')
% 2. 파이 차트를 그리고, 차트 객체를 'p' 변수에 저장
p = piechart(LC9_Total_in_Thr, "Synapse", "Type", 'LabelStyle', 'name');
title(sprintf('LC9 Inputs Total (synapse > %d)', in_Thr));

% 3. 생성된 pieColors 행렬을 파이 차트 객체의 FaceColor 속성에 할당
p.ColorOrder = pieColors;
%%
Outdata = [sum(LC9_out_central.Synapse),sum(LC9_out_optic.Synapse)];
figure;set(gcf,'Color','w')
piechart(Outdata,Names={'Central','Optic'})
title('LC9 Output');

out_Thr=800;
LC9_Total_out_Thr=LC9_Total_out;
LC9_Total_out_Thr(LC9_Total_out_Thr.Synapse<out_Thr,:)=[];

for i=1:1:size(LC9_Total_out_Thr,1)
    current_types=LC9_Total_out_Thr.Type{i};
    current_types_root_ids=FAFB_consolidated_cell_types.root_id(strcmp(FAFB_consolidated_cell_types.primary_type,current_types));
    superclasses=FAFB_classfication.super_class(ismember(FAFB_classfication.root_id,current_types_root_ids));
    [u, ~, idx] = unique(superclasses);           % 고유 문자열 u와 index
    counts = histcounts(idx, 1:numel(u)+1); % 각 문자열 빈도
    [~, maxIdx] = max(counts);         % 최빈값 인덱스
    mostCommon_superclasses = superclasses{maxIdx};            % 최빈 문자열
    LC9_Total_out_Thr.superclass{i}=mostCommon_superclasses;
end

LC9_Total_out_Thr = sortrows(LC9_Total_out_Thr,"superclass","ascend");
LC9_Total_out_Thr = sortrows(LC9_Total_out_Thr,"Type","ascend");
LC9_Total_out_Thr = sortrows(LC9_Total_out_Thr,"superclass","ascend");
% p=piechart(LC9_in_central_Thr,"Synapse","Type",'LabelStyle','name');
% 
% %% Pie chart with custom colors based on superclass

% 1. LC9_in_central_Thr 테이블의 순서에 맞는 색상 행렬(pieColors) 생성
% 테이블의 행 개수만큼 3열(RGB)짜리 행렬을 0으로 초기화
numSlices = height(LC9_Total_out_Thr);
pieColors = zeros(numSlices, 3);

% for 루프를 사용해 테이블의 각 행을 순회
for i = 1:1:numSlices
    % 현재 행의 superclass를 가져옴
    current_superclass = LC9_Total_out_Thr.superclass{i};
    
    % uniqueSuperclass 목록에서 현재 superclass와 일치하는 위치(인덱스)를 찾음
    % matches 함수는 strcmp와 유사하게 문자열을 비교하여 논리형 배열을 반환
    color_index = matches(uniqueSuperclass, current_superclass);
    
    % 찾은 인덱스를 사용해 미리 정의된 colors 행렬에서 해당 색상(RGB 값)을 가져옴
    % pieColors의 현재 행에 색상 값을 할당
    pieColors(i, :) = colors(color_index, :);
end
figure;set(gcf,'Color','w')
% 2. 파이 차트를 그리고, 차트 객체를 'p' 변수에 저장
p = piechart(LC9_Total_out_Thr, "Synapse", "Type", 'LabelStyle', 'name');
title(sprintf('LC9 Outputs Total (synapse > %d)', out_Thr));

% 3. 생성된 pieColors 행렬을 파이 차트 객체의 FaceColor 속성에 할당
p.ColorOrder = pieColors;

%% in_central LC9
in_central_Thr=800;
LC9_in_central_Thr=LC9_in_central;
LC9_in_central_Thr(LC9_in_central_Thr.Synapse<in_central_Thr,:)=[];

for i=1:1:size(LC9_in_central_Thr,1)
    current_types=LC9_in_central_Thr.Type{i};
    current_types_root_ids=FAFB_consolidated_cell_types.root_id(strcmp(FAFB_consolidated_cell_types.primary_type,current_types));
    superclasses=FAFB_classfication.super_class(ismember(FAFB_classfication.root_id,current_types_root_ids));
    [u, ~, idx] = unique(superclasses);           % 고유 문자열 u와 index
    counts = histcounts(idx, 1:numel(u)+1); % 각 문자열 빈도
    [~, maxIdx] = max(counts);         % 최빈값 인덱스
    mostCommon_superclasses = superclasses{maxIdx};            % 최빈 문자열
    LC9_in_central_Thr.superclass{i}=mostCommon_superclasses;
end

LC9_in_central_Thr = sortrows(LC9_in_central_Thr,"superclass","ascend");
LC9_in_central_Thr = sortrows(LC9_in_central_Thr,"Type","ascend");
LC9_in_central_Thr = sortrows(LC9_in_central_Thr,"superclass","ascend");
% p=piechart(LC9_in_central_Thr,"Synapse","Type",'LabelStyle','name');
% 
% %% Pie chart with custom colors based on superclass

% 1. LC9_in_central_Thr 테이블의 순서에 맞는 색상 행렬(pieColors) 생성
% 테이블의 행 개수만큼 3열(RGB)짜리 행렬을 0으로 초기화
numSlices = height(LC9_in_central_Thr);
pieColors = zeros(numSlices, 3);

% for 루프를 사용해 테이블의 각 행을 순회
for i = 1:1:numSlices
    % 현재 행의 superclass를 가져옴
    current_superclass = LC9_in_central_Thr.superclass{i};
    
    % uniqueSuperclass 목록에서 현재 superclass와 일치하는 위치(인덱스)를 찾음
    % matches 함수는 strcmp와 유사하게 문자열을 비교하여 논리형 배열을 반환
    color_index = matches(uniqueSuperclass, current_superclass);
    
    % 찾은 인덱스를 사용해 미리 정의된 colors 행렬에서 해당 색상(RGB 값)을 가져옴
    % pieColors의 현재 행에 색상 값을 할당
    pieColors(i, :) = colors(color_index, :);
end
figure;set(gcf,'Color','w')
% 2. 파이 차트를 그리고, 차트 객체를 'p' 변수에 저장
p = piechart(LC9_in_central_Thr, "Synapse", "Type", 'LabelStyle', 'name');
title(sprintf('LC9 Inputs from Central brains (synapse > %d)', in_central_Thr));

% 3. 생성된 pieColors 행렬을 파이 차트 객체의 FaceColor 속성에 할당
p.ColorOrder = pieColors;

%% in_optic LC9

in_optic_Thr=1100;
LC9_in_optic_Thr=LC9_in_optic;
LC9_in_optic_Thr(LC9_in_optic_Thr.Synapse<in_optic_Thr,:)=[];

for i=1:1:size(LC9_in_optic_Thr,1)
    current_types=LC9_in_optic_Thr.Type{i};
    current_types_root_ids=FAFB_consolidated_cell_types.root_id(strcmp(FAFB_consolidated_cell_types.primary_type,current_types));
    superclasses=FAFB_classfication.super_class(ismember(FAFB_classfication.root_id,current_types_root_ids));
    [u, ~, idx] = unique(superclasses);           % 고유 문자열 u와 index
    counts = histcounts(idx, 1:numel(u)+1); % 각 문자열 빈도
    [~, maxIdx] = max(counts);         % 최빈값 인덱스
    mostCommon_superclasses = superclasses{maxIdx};            % 최빈 문자열
    LC9_in_optic_Thr.superclass{i}=mostCommon_superclasses;
end

LC9_in_optic_Thr = sortrows(LC9_in_optic_Thr,"superclass","ascend");
LC9_in_optic_Thr = sortrows(LC9_in_optic_Thr,"Type","ascend");
LC9_in_optic_Thr = sortrows(LC9_in_optic_Thr,"superclass","ascend");
% p=piechart(LC9_in_optic_Thr,"Synapse","Type",'LabelStyle','name');
% 
% %% Pie chart with custom colors based on superclass

% 1. LC9_in_optic_Thr 테이블의 순서에 맞는 색상 행렬(pieColors) 생성
% 테이블의 행 개수만큼 3열(RGB)짜리 행렬을 0으로 초기화
numSlices = height(LC9_in_optic_Thr);
pieColors = zeros(numSlices, 3);

% for 루프를 사용해 테이블의 각 행을 순회
for i = 1:1:numSlices
    % 현재 행의 superclass를 가져옴
    current_superclass = LC9_in_optic_Thr.superclass{i};
    
    % uniqueSuperclass 목록에서 현재 superclass와 일치하는 위치(인덱스)를 찾음
    % matches 함수는 strcmp와 유사하게 문자열을 비교하여 논리형 배열을 반환
    color_index = matches(uniqueSuperclass, current_superclass);
    
    % 찾은 인덱스를 사용해 미리 정의된 colors 행렬에서 해당 색상(RGB 값)을 가져옴
    % pieColors의 현재 행에 색상 값을 할당
    pieColors(i, :) = colors(color_index, :);
end
figure;set(gcf,'Color','w')

% 2. 파이 차트를 그리고, 차트 객체를 'p' 변수에 저장
p = piechart(LC9_in_optic_Thr, "Synapse", "Type", 'LabelStyle', 'name');
title(sprintf('LC9 Inputs from Optic lobes (synapse > %d)', in_optic_Thr));

% 3. 생성된 pieColors 행렬을 파이 차트 객체의 FaceColor 속성에 할당
p.ColorOrder = pieColors;

%% out_central
out_central_Thr=700;
LC9_out_central_Thr=LC9_out_central;
LC9_out_central_Thr(LC9_out_central_Thr.Synapse<out_central_Thr,:)=[];

for i=1:1:size(LC9_out_central_Thr,1)
    current_types=LC9_out_central_Thr.Type{i};
    current_types_root_ids=FAFB_consolidated_cell_types.root_id(strcmp(FAFB_consolidated_cell_types.primary_type,current_types));
    superclasses=FAFB_classfication.super_class(ismember(FAFB_classfication.root_id,current_types_root_ids));
    [u, ~, idx] = unique(superclasses);           % 고유 문자열 u와 outdex
    counts = histcounts(idx, 1:numel(u)+1); % 각 문자열 빈도
    [~, maxIdx] = max(counts);         % 최빈값 인덱스
    mostCommon_superclasses = superclasses{maxIdx};            % 최빈 문자열
    LC9_out_central_Thr.superclass{i}=mostCommon_superclasses;
end

LC9_out_central_Thr = sortrows(LC9_out_central_Thr,"superclass","ascend");
LC9_out_central_Thr = sortrows(LC9_out_central_Thr,"Type","ascend");
LC9_out_central_Thr = sortrows(LC9_out_central_Thr,"superclass","ascend");
% 1. LC9_out_central_Thr 테이블의 순서에 맞는 색상 행렬(pieColors) 생성
% 테이블의 행 개수만큼 3열(RGB)짜리 행렬을 0으로 초기화
numSlices = height(LC9_out_central_Thr);
pieColors = zeros(numSlices, 3);

% for 루프를 사용해 테이블의 각 행을 순회
for i = 1:1:numSlices
    % 현재 행의 superclass를 가져옴
    current_superclass = LC9_out_central_Thr.superclass{i};
    
    % uniqueSuperclass 목록에서 현재 superclass와 일치하는 위치(인덱스)를 찾음
    % matches 함수는 strcmp와 유사하게 문자열을 비교하여 논리형 배열을 반환
    color_index = matches(uniqueSuperclass, current_superclass);
    
    % 찾은 인덱스를 사용해 미리 정의된 colors 행렬에서 해당 색상(RGB 값)을 가져옴
    % pieColors의 현재 행에 색상 값을 할당
    pieColors(i, :) = colors(color_index, :);
end
figure;set(gcf,'Color','w')

% 2. 파이 차트를 그리고, 차트 객체를 'p' 변수에 저장
p = piechart(LC9_out_central_Thr, "Synapse", "Type", 'LabelStyle', 'name');
title(sprintf('LC9 Outputs to Central brains (synapse > %d)', out_central_Thr));

% 3. 생성된 pieColors 행렬을 파이 차트 객체의 FaceColor 속성에 할당
p.ColorOrder = pieColors;

%% out_optic
out_optic_Thr=500;
LC9_out_optic_Thr=LC9_out_optic;
LC9_out_optic_Thr(LC9_out_optic_Thr.Synapse<out_optic_Thr,:)=[];

for i=1:1:size(LC9_out_optic_Thr,1)
    current_types=LC9_out_optic_Thr.Type{i};
    current_types_root_ids=FAFB_consolidated_cell_types.root_id(strcmp(FAFB_consolidated_cell_types.primary_type,current_types));
    superclasses=FAFB_classfication.super_class(ismember(FAFB_classfication.root_id,current_types_root_ids));
    [u, ~, idx] = unique(superclasses);           % 고유 문자열 u와 index
    counts = histcounts(idx, 1:numel(u)+1); % 각 문자열 빈도
    [~, maxIdx] = max(counts);         % 최빈값 인덱스
    mostCommon_superclasses = superclasses{maxIdx};            % 최빈 문자열
    LC9_out_optic_Thr.superclass{i}=mostCommon_superclasses;
end

LC9_out_optic_Thr = sortrows(LC9_out_optic_Thr,"superclass","ascend");
LC9_out_optic_Thr = sortrows(LC9_out_optic_Thr,"Type","ascend");
LC9_out_optic_Thr = sortrows(LC9_out_optic_Thr,"superclass","ascend");
% p=piechart(LC9_in_optic_Thr,"Synapse","Type",'LabelStyle','name');
% 
% %% Pie chart with custom colors based on superclass

% 1. LC9_in_optic_Thr 테이블의 순서에 맞는 색상 행렬(pieColors) 생성
% 테이블의 행 개수만큼 3열(RGB)짜리 행렬을 0으로 초기화
numSlices = height(LC9_out_optic_Thr);
pieColors = zeros(numSlices, 3);

% for 루프를 사용해 테이블의 각 행을 순회
for i = 1:1:numSlices
    % 현재 행의 superclass를 가져옴
    current_superclass = LC9_out_optic_Thr.superclass{i};
    
    % uniqueSuperclass 목록에서 현재 superclass와 일치하는 위치(인덱스)를 찾음
    % matches 함수는 strcmp와 유사하게 문자열을 비교하여 논리형 배열을 반환
    color_index = matches(uniqueSuperclass, current_superclass);
    
    % 찾은 인덱스를 사용해 미리 정의된 colors 행렬에서 해당 색상(RGB 값)을 가져옴
    % pieColors의 현재 행에 색상 값을 할당
    pieColors(i, :) = colors(color_index, :);
end
figure;set(gcf,'Color','w')

% 2. 파이 차트를 그리고, 차트 객체를 'p' 변수에 저장
p = piechart(LC9_out_optic_Thr, "Synapse", "Type", 'LabelStyle', 'name');
title(sprintf('LC9 Outputs to Optic lobes (synapse > %d)', out_optic_Thr));

% 3. 생성된 pieColors 행렬을 파이 차트 객체의 FaceColor 속성에 할당
p.ColorOrder = pieColors;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
LT43_Total_in=sum(Wantsee.InConnections{2}.syn_count);
LT43_Total_out=sum(Wantsee.OutConnections{2}.syn_count);
LT43_in_optic=Wantsee.unique_InOptic_Types{2};
LT43_in_optic = cell2table(LT43_in_optic, 'VariableNames', {'Type','Synapse','Number'});
LT43_in_central=Wantsee.unique_InCentral_Types{2};
LT43_in_central = cell2table(LT43_in_central, 'VariableNames', {'Type','Synapse','Number'});
LT43_out_optic=Wantsee.unique_OutOptic_Types{2};
LT43_out_optic = cell2table(LT43_out_optic, 'VariableNames', {'Type','Synapse','Number'});
LT43_out_central=Wantsee.unique_OutCentral_Types{2};
LT43_out_central = cell2table(LT43_out_central, 'VariableNames', {'Type','Synapse','Number'});


%%
Indata = [sum(LT43_in_central.Synapse),sum(LT43_in_optic.Synapse)];
figure;set(gcf,'Color','w')
piechart(Indata,Names={'Central','Optic'})
title('LT43 Input');
%%
Outdata = [sum(LT43_out_central.Synapse),sum(LT43_out_optic.Synapse)];
figure;set(gcf,'Color','w')
piechart(Outdata,Names={'Central','Optic'})
title('LT43 Output');
%% === (추가) LT43: 전체 In/Out 타입별 piechart ===
LT43_Total_in = Wantsee.InNeuronTypes{2};
LT43_Total_in = cell2table(LT43_Total_in(:,[1 2]), 'VariableNames', {'Type','Synapse'});

LT43_Total_out = Wantsee.OutNeuronTypes{2};
LT43_Total_out = cell2table(LT43_Total_out(:,[1 2]), 'VariableNames', {'Type','Synapse'});

% 전체 Input 타입별 pie (임계값 직접 조정)
in_Thr = 40;   % <-- 네가 조정
LT43_Total_in_Thr = LT43_Total_in;
LT43_Total_in_Thr(LT43_Total_in_Thr.Synapse < in_Thr, :) = [];

for i=1:height(LT43_Total_in_Thr)
    current_types = LT43_Total_in_Thr.Type{i};
    current_types_root_ids = FAFB_consolidated_cell_types.root_id(strcmp(FAFB_consolidated_cell_types.primary_type,current_types));
    superclasses = FAFB_classfication.super_class(ismember(FAFB_classfication.root_id,current_types_root_ids));
    [u,~,idx] = unique(superclasses); counts = histcounts(idx,1:numel(u)+1);
    [~,maxIdx] = max(counts);
    LT43_Total_in_Thr.superclass{i} = superclasses{maxIdx};
end

LT43_Total_in_Thr = sortrows(LT43_Total_in_Thr,"superclass","ascend");
LT43_Total_in_Thr = sortrows(LT43_Total_in_Thr,"Type","ascend");
LT43_Total_in_Thr = sortrows(LT43_Total_in_Thr,"superclass","ascend");

numSlices = height(LT43_Total_in_Thr);
pieColors = zeros(numSlices,3);
for i=1:numSlices
    color_index = matches(uniqueSuperclass, LT43_Total_in_Thr.superclass{i});
    pieColors(i,:) = colors(color_index,:);
end
figure; set(gcf,'Color','w')
p = piechart(LT43_Total_in_Thr, "Synapse", "Type", 'LabelStyle','name');
title(sprintf('LT43 Inputs Total (synapse > %d)', in_Thr));
p.ColorOrder = pieColors;

% 전체 Output 타입별 pie (임계값 직접 조정)
out_Thr = 80;  % <-- 네가 조정
LT43_Total_out_Thr = LT43_Total_out;
LT43_Total_out_Thr(LT43_Total_out_Thr.Synapse < out_Thr, :) = [];

for i=1:height(LT43_Total_out_Thr)
    current_types = LT43_Total_out_Thr.Type{i};
    current_types_root_ids = FAFB_consolidated_cell_types.root_id(strcmp(FAFB_consolidated_cell_types.primary_type,current_types));
    superclasses = FAFB_classfication.super_class(ismember(FAFB_classfication.root_id,current_types_root_ids));
    [u,~,idx] = unique(superclasses); counts = histcounts(idx,1:numel(u)+1);
    [~,maxIdx] = max(counts);
    LT43_Total_out_Thr.superclass{i} = superclasses{maxIdx};
end

LT43_Total_out_Thr = sortrows(LT43_Total_out_Thr,"superclass","ascend");
LT43_Total_out_Thr = sortrows(LT43_Total_out_Thr,"Type","ascend");
LT43_Total_out_Thr = sortrows(LT43_Total_out_Thr,"superclass","ascend");

numSlices = height(LT43_Total_out_Thr);
pieColors = zeros(numSlices,3);
for i=1:numSlices
    color_index = matches(uniqueSuperclass, LT43_Total_out_Thr.superclass{i});
    pieColors(i,:) = colors(color_index,:);
end
figure; set(gcf,'Color','w')
p = piechart(LT43_Total_out_Thr, "Synapse", "Type", 'LabelStyle','name');
title(sprintf('LT43 Outputs Total (synapse > %d)', out_Thr));
p.ColorOrder = pieColors;

%% in_central LT43
in_central_Thr=20;
LT43_in_central_Thr=LT43_in_central;
LT43_in_central_Thr(LT43_in_central_Thr.Synapse<in_central_Thr,:)=[];

for i=1:1:size(LT43_in_central_Thr,1)
    current_types=LT43_in_central_Thr.Type{i};
    current_types_root_ids=FAFB_consolidated_cell_types.root_id(strcmp(FAFB_consolidated_cell_types.primary_type,current_types));
    superclasses=FAFB_classfication.super_class(ismember(FAFB_classfication.root_id,current_types_root_ids));
    [u, ~, idx] = unique(superclasses);           % 고유 문자열 u와 index
    counts = histcounts(idx, 1:numel(u)+1); % 각 문자열 빈도
    [~, maxIdx] = max(counts);         % 최빈값 인덱스
    mostCommon_superclasses = superclasses{maxIdx};            % 최빈 문자열
    LT43_in_central_Thr.superclass{i}=mostCommon_superclasses;
end

LT43_in_central_Thr = sortrows(LT43_in_central_Thr,"superclass","ascend");
LT43_in_central_Thr = sortrows(LT43_in_central_Thr,"Type","ascend");
LT43_in_central_Thr = sortrows(LT43_in_central_Thr,"superclass","ascend");
% p=piechart(LT43_in_central_Thr,"Synapse","Type",'LabelStyle','name');
% 
% %% Pie chart with custom colors based on superclass

% 1. LT43_in_central_Thr 테이블의 순서에 맞는 색상 행렬(pieColors) 생성
% 테이블의 행 개수만큼 3열(RGB)짜리 행렬을 0으로 초기화
numSlices = height(LT43_in_central_Thr);
pieColors = zeros(numSlices, 3);

% for 루프를 사용해 테이블의 각 행을 순회
for i = 1:1:numSlices
    % 현재 행의 superclass를 가져옴
    current_superclass = LT43_in_central_Thr.superclass{i};
    
    % uniqueSuperclass 목록에서 현재 superclass와 일치하는 위치(인덱스)를 찾음
    % matches 함수는 strcmp와 유사하게 문자열을 비교하여 논리형 배열을 반환
    color_index = matches(uniqueSuperclass, current_superclass);
    
    % 찾은 인덱스를 사용해 미리 정의된 colors 행렬에서 해당 색상(RGB 값)을 가져옴
    % pieColors의 현재 행에 색상 값을 할당
    pieColors(i, :) = colors(color_index, :);
end
figure;set(gcf,'Color','w')
% 2. 파이 차트를 그리고, 차트 객체를 'p' 변수에 저장
p = piechart(LT43_in_central_Thr, "Synapse", "Type", 'LabelStyle', 'name');
title(sprintf('LT43 Inputs from Central brains (synapse > %d)', in_central_Thr));

% 3. 생성된 pieColors 행렬을 파이 차트 객체의 FaceColor 속성에 할당
p.ColorOrder = pieColors;

%% in_optic LT43

in_optic_Thr=50;
LT43_in_optic_Thr=LT43_in_optic;
LT43_in_optic_Thr(LT43_in_optic_Thr.Synapse<in_optic_Thr,:)=[];

for i=1:1:size(LT43_in_optic_Thr,1)
    current_types=LT43_in_optic_Thr.Type{i};
    current_types_root_ids=FAFB_consolidated_cell_types.root_id(strcmp(FAFB_consolidated_cell_types.primary_type,current_types));
    superclasses=FAFB_classfication.super_class(ismember(FAFB_classfication.root_id,current_types_root_ids));
    [u, ~, idx] = unique(superclasses);           % 고유 문자열 u와 index
    counts = histcounts(idx, 1:numel(u)+1); % 각 문자열 빈도
    [~, maxIdx] = max(counts);         % 최빈값 인덱스
    mostCommon_superclasses = superclasses{maxIdx};            % 최빈 문자열
    LT43_in_optic_Thr.superclass{i}=mostCommon_superclasses;
end

LT43_in_optic_Thr = sortrows(LT43_in_optic_Thr,"superclass","ascend");
LT43_in_optic_Thr = sortrows(LT43_in_optic_Thr,"Type","ascend");
LT43_in_optic_Thr = sortrows(LT43_in_optic_Thr,"superclass","ascend");
% p=piechart(LT43_in_optic_Thr,"Synapse","Type",'LabelStyle','name');
% 
% %% Pie chart with custom colors based on superclass

% 1. LT43_in_optic_Thr 테이블의 순서에 맞는 색상 행렬(pieColors) 생성
% 테이블의 행 개수만큼 3열(RGB)짜리 행렬을 0으로 초기화
numSlices = height(LT43_in_optic_Thr);
pieColors = zeros(numSlices, 3);

% for 루프를 사용해 테이블의 각 행을 순회
for i = 1:1:numSlices
    % 현재 행의 superclass를 가져옴
    current_superclass = LT43_in_optic_Thr.superclass{i};
    
    % uniqueSuperclass 목록에서 현재 superclass와 일치하는 위치(인덱스)를 찾음
    % matches 함수는 strcmp와 유사하게 문자열을 비교하여 논리형 배열을 반환
    color_index = matches(uniqueSuperclass, current_superclass);
    
    % 찾은 인덱스를 사용해 미리 정의된 colors 행렬에서 해당 색상(RGB 값)을 가져옴
    % pieColors의 현재 행에 색상 값을 할당
    pieColors(i, :) = colors(color_index, :);
end
figure;set(gcf,'Color','w')

% 2. 파이 차트를 그리고, 차트 객체를 'p' 변수에 저장
p = piechart(LT43_in_optic_Thr, "Synapse", "Type", 'LabelStyle', 'name');
title(sprintf('LT43 Inputs from Optic lobes (synapse > %d)', in_optic_Thr));

% 3. 생성된 pieColors 행렬을 파이 차트 객체의 FaceColor 속성에 할당
p.ColorOrder = pieColors;

%% out_central
out_central_Thr=45;
LT43_out_central_Thr=LT43_out_central;
LT43_out_central_Thr(LT43_out_central_Thr.Synapse<out_central_Thr,:)=[];

for i=1:1:size(LT43_out_central_Thr,1)
    current_types=LT43_out_central_Thr.Type{i};
    current_types_root_ids=FAFB_consolidated_cell_types.root_id(strcmp(FAFB_consolidated_cell_types.primary_type,current_types));
    superclasses=FAFB_classfication.super_class(ismember(FAFB_classfication.root_id,current_types_root_ids));
    [u, ~, idx] = unique(superclasses);           % 고유 문자열 u와 outdex
    counts = histcounts(idx, 1:numel(u)+1); % 각 문자열 빈도
    [~, maxIdx] = max(counts);         % 최빈값 인덱스
    mostCommon_superclasses = superclasses{maxIdx};            % 최빈 문자열
    LT43_out_central_Thr.superclass{i}=mostCommon_superclasses;
end

LT43_out_central_Thr = sortrows(LT43_out_central_Thr,"superclass","ascend");
LT43_out_central_Thr = sortrows(LT43_out_central_Thr,"Type","ascend");
LT43_out_central_Thr = sortrows(LT43_out_central_Thr,"superclass","ascend");
% 1. LT43_out_central_Thr 테이블의 순서에 맞는 색상 행렬(pieColors) 생성
% 테이블의 행 개수만큼 3열(RGB)짜리 행렬을 0으로 초기화
numSlices = height(LT43_out_central_Thr);
pieColors = zeros(numSlices, 3);

% for 루프를 사용해 테이블의 각 행을 순회
for i = 1:1:numSlices
    % 현재 행의 superclass를 가져옴
    current_superclass = LT43_out_central_Thr.superclass{i};
    
    % uniqueSuperclass 목록에서 현재 superclass와 일치하는 위치(인덱스)를 찾음
    % matches 함수는 strcmp와 유사하게 문자열을 비교하여 논리형 배열을 반환
    color_index = matches(uniqueSuperclass, current_superclass);
    
    % 찾은 인덱스를 사용해 미리 정의된 colors 행렬에서 해당 색상(RGB 값)을 가져옴
    % pieColors의 현재 행에 색상 값을 할당
    pieColors(i, :) = colors(color_index, :);
end
figure;set(gcf,'Color','w')

% 2. 파이 차트를 그리고, 차트 객체를 'p' 변수에 저장
p = piechart(LT43_out_central_Thr, "Synapse", "Type", 'LabelStyle', 'name');
title(sprintf('LT43 Outputs to Central Neurons (synapse > %d)', out_central_Thr));

% 3. 생성된 pieColors 행렬을 파이 차트 객체의 FaceColor 속성에 할당
p.ColorOrder = pieColors;

%% out_optic
out_optic_Thr=60;
LT43_out_optic_Thr=LT43_out_optic;
LT43_out_optic_Thr(LT43_out_optic_Thr.Synapse<out_optic_Thr,:)=[];

for i=1:1:size(LT43_out_optic_Thr,1)
    current_types=LT43_out_optic_Thr.Type{i};
    current_types_root_ids=FAFB_consolidated_cell_types.root_id(strcmp(FAFB_consolidated_cell_types.primary_type,current_types));
    superclasses=FAFB_classfication.super_class(ismember(FAFB_classfication.root_id,current_types_root_ids));
    [u, ~, idx] = unique(superclasses);           % 고유 문자열 u와 index
    counts = histcounts(idx, 1:numel(u)+1); % 각 문자열 빈도
    [~, maxIdx] = max(counts);         % 최빈값 인덱스
    mostCommon_superclasses = superclasses{maxIdx};            % 최빈 문자열
    LT43_out_optic_Thr.superclass{i}=mostCommon_superclasses;
end

LT43_out_optic_Thr = sortrows(LT43_out_optic_Thr,"superclass","ascend");
LT43_out_optic_Thr = sortrows(LT43_out_optic_Thr,"Type","ascend");
LT43_out_optic_Thr = sortrows(LT43_out_optic_Thr,"superclass","ascend");
% p=piechart(LT43_in_optic_Thr,"Synapse","Type",'LabelStyle','name');
% 
% %% Pie chart with custom colors based on superclass

% 1. LT43_in_optic_Thr 테이블의 순서에 맞는 색상 행렬(pieColors) 생성
% 테이블의 행 개수만큼 3열(RGB)짜리 행렬을 0으로 초기화
numSlices = height(LT43_out_optic_Thr);
pieColors = zeros(numSlices, 3);

% for 루프를 사용해 테이블의 각 행을 순회
for i = 1:1:numSlices
    % 현재 행의 superclass를 가져옴
    current_superclass = LT43_out_optic_Thr.superclass{i};
    
    % uniqueSuperclass 목록에서 현재 superclass와 일치하는 위치(인덱스)를 찾음
    % matches 함수는 strcmp와 유사하게 문자열을 비교하여 논리형 배열을 반환
    color_index = matches(uniqueSuperclass, current_superclass);
    
    % 찾은 인덱스를 사용해 미리 정의된 colors 행렬에서 해당 색상(RGB 값)을 가져옴
    % pieColors의 현재 행에 색상 값을 할당
    pieColors(i, :) = colors(color_index, :);
end
figure;set(gcf,'Color','w')

% 2. 파이 차트를 그리고, 차트 객체를 'p' 변수에 저장
p = piechart(LT43_out_optic_Thr, "Synapse", "Type", 'LabelStyle', 'name');
title(sprintf('LT43 Outputs to Optic lobes (synapse > %d)', out_optic_Thr));

% 3. 생성된 pieColors 행렬을 파이 차트 객체의 FaceColor 속성에 할당
p.ColorOrder = pieColors;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% =========================
%% LT52
%% =========================
LT52_Total_in  = sum(Wantsee.InConnections{3}.syn_count);
LT52_Total_out = sum(Wantsee.OutConnections{3}.syn_count);
LT52_in_optic   = Wantsee.unique_InOptic_Types{3};
LT52_in_optic   = cell2table(LT52_in_optic,   'VariableNames', {'Type','Synapse','Number'});
LT52_in_central = Wantsee.unique_InCentral_Types{3};
LT52_in_central = cell2table(LT52_in_central, 'VariableNames', {'Type','Synapse','Number'});
LT52_out_optic  = Wantsee.unique_OutOptic_Types{3};
LT52_out_optic  = cell2table(LT52_out_optic,  'VariableNames', {'Type','Synapse','Number'});
LT52_out_central= Wantsee.unique_OutCentral_Types{3};
LT52_out_central= cell2table(LT52_out_central,'VariableNames', {'Type','Synapse','Number'});

% 전체 입력/출력 파이
Indata = [sum(LT52_in_central.Synapse), sum(LT52_in_optic.Synapse)];
figure; set(gcf,'Color','w'); piechart(Indata, Names={'Central','Optic'}); title('LT52 Input');

Outdata = [sum(LT52_out_central.Synapse), sum(LT52_out_optic.Synapse)];
figure; set(gcf,'Color','w'); piechart(Outdata, Names={'Central','Optic'}); title('LT52 Output');
%% === (추가) LT52: 전체 In/Out 타입별 piechart ===
LT52_Total_in_tbl = Wantsee.InNeuronTypes{3};
LT52_Total_in_tbl = cell2table(LT52_Total_in_tbl(:,[1 2]), 'VariableNames', {'Type','Synapse'});

LT52_Total_out_tbl = Wantsee.OutNeuronTypes{3};
LT52_Total_out_tbl = cell2table(LT52_Total_out_tbl(:,[1 2]), 'VariableNames', {'Type','Synapse'});

% 전체 Input 타입별 pie
in_Thr = 700; % <-- 네가 조정
LT52_Total_in_Thr = LT52_Total_in_tbl;
LT52_Total_in_Thr(LT52_Total_in_Thr.Synapse < in_Thr, :) = [];

for i=1:height(LT52_Total_in_Thr)
    current_types = LT52_Total_in_Thr.Type{i};
    current_types_root_ids = FAFB_consolidated_cell_types.root_id(strcmp(FAFB_consolidated_cell_types.primary_type,current_types));
    superclasses = FAFB_classfication.super_class(ismember(FAFB_classfication.root_id,current_types_root_ids));
    [u,~,idx] = unique(superclasses); counts = histcounts(idx,1:numel(u)+1);
    [~,maxIdx] = max(counts);
    LT52_Total_in_Thr.superclass{i} = superclasses{maxIdx};
end

LT52_Total_in_Thr = sortrows(LT52_Total_in_Thr,"superclass","ascend");
LT52_Total_in_Thr = sortrows(LT52_Total_in_Thr,"Type","ascend");
LT52_Total_in_Thr = sortrows(LT52_Total_in_Thr,"superclass","ascend");

numSlices = height(LT52_Total_in_Thr);
pieColors = zeros(numSlices,3);
for i=1:numSlices
    color_index = matches(uniqueSuperclass, LT52_Total_in_Thr.superclass{i});
    pieColors(i,:) = colors(color_index,:);
end
figure; set(gcf,'Color','w')
p = piechart(LT52_Total_in_Thr, "Synapse", "Type", 'LabelStyle','name');
title(sprintf('LT52 Inputs Total (synapse > %d)', in_Thr));
p.ColorOrder = pieColors;

% 전체 Output 타입별 pie
out_Thr = 500; % <-- 네가 조정
LT52_Total_out_Thr = LT52_Total_out_tbl;
LT52_Total_out_Thr(LT52_Total_out_Thr.Synapse < out_Thr, :) = [];

for i=1:height(LT52_Total_out_Thr)
    current_types = LT52_Total_out_Thr.Type{i};
    current_types_root_ids = FAFB_consolidated_cell_types.root_id(strcmp(FAFB_consolidated_cell_types.primary_type,current_types));
    superclasses = FAFB_classfication.super_class(ismember(FAFB_classfication.root_id,current_types_root_ids));
    [u,~,idx] = unique(superclasses); counts = histcounts(idx,1:numel(u)+1);
    [~,maxIdx] = max(counts);
    LT52_Total_out_Thr.superclass{i} = superclasses{maxIdx};
end

LT52_Total_out_Thr = sortrows(LT52_Total_out_Thr,"superclass","ascend");
LT52_Total_out_Thr = sortrows(LT52_Total_out_Thr,"Type","ascend");
LT52_Total_out_Thr = sortrows(LT52_Total_out_Thr,"superclass","ascend");

numSlices = height(LT52_Total_out_Thr);
pieColors = zeros(numSlices,3);
for i=1:numSlices
    color_index = matches(uniqueSuperclass, LT52_Total_out_Thr.superclass{i});
    pieColors(i,:) = colors(color_index,:);
end
figure; set(gcf,'Color','w')
p = piechart(LT52_Total_out_Thr, "Synapse", "Type", 'LabelStyle','name');
title(sprintf('LT52 Outputs Total (synapse > %d)', out_Thr));
p.ColorOrder = pieColors;

%% in_central LT52
in_central_Thr = 160;  % <- 직접 조정
LT52_in_central_Thr = LT52_in_central;
LT52_in_central_Thr(LT52_in_central_Thr.Synapse < in_central_Thr, :) = [];

for i=1:size(LT52_in_central_Thr,1)
    current_types = LT52_in_central_Thr.Type{i};
    current_types_root_ids = FAFB_consolidated_cell_types.root_id(strcmp(FAFB_consolidated_cell_types.primary_type,current_types));
    superclasses = FAFB_classfication.super_class(ismember(FAFB_classfication.root_id,current_types_root_ids));
    [u,~,idx] = unique(superclasses); counts = histcounts(idx,1:numel(u)+1);
    [~,maxIdx] = max(counts); LT52_in_central_Thr.superclass{i} = superclasses{maxIdx};
end
LT52_in_central_Thr = sortrows(LT52_in_central_Thr,"superclass","ascend");
LT52_in_central_Thr = sortrows(LT52_in_central_Thr,"Type","ascend");
LT52_in_central_Thr = sortrows(LT52_in_central_Thr,"superclass","ascend");
numSlices = height(LT52_in_central_Thr); pieColors = zeros(numSlices,3);
for i=1:numSlices
    color_index = matches(uniqueSuperclass, LT52_in_central_Thr.superclass{i});
    pieColors(i,:) = colors(color_index,:);
end
figure; set(gcf,'Color','w');
p = piechart(LT52_in_central_Thr, "Synapse", "Type", 'LabelStyle','name');
title(sprintf('LT52 Inputs from Central brains (synapse > %d)', in_central_Thr));
p.ColorOrder = pieColors;

%% in_optic LT52
in_optic_Thr = 700;  % <- 직접 조정
LT52_in_optic_Thr = LT52_in_optic;
LT52_in_optic_Thr(LT52_in_optic_Thr.Synapse < in_optic_Thr, :) = [];

for i=1:size(LT52_in_optic_Thr,1)
    current_types = LT52_in_optic_Thr.Type{i};
    current_types_root_ids = FAFB_consolidated_cell_types.root_id(strcmp(FAFB_consolidated_cell_types.primary_type,current_types));
    superclasses = FAFB_classfication.super_class(ismember(FAFB_classfication.root_id,current_types_root_ids));
    [u,~,idx] = unique(superclasses); counts = histcounts(idx,1:numel(u)+1);
    [~,maxIdx] = max(counts); LT52_in_optic_Thr.superclass{i} = superclasses{maxIdx};
end
LT52_in_optic_Thr = sortrows(LT52_in_optic_Thr,"superclass","ascend");
LT52_in_optic_Thr = sortrows(LT52_in_optic_Thr,"Type","ascend");
LT52_in_optic_Thr = sortrows(LT52_in_optic_Thr,"superclass","ascend");
numSlices = height(LT52_in_optic_Thr); pieColors = zeros(numSlices,3);
for i=1:numSlices
    color_index = matches(uniqueSuperclass, LT52_in_optic_Thr.superclass{i});
    pieColors(i,:) = colors(color_index,:);
end
figure; set(gcf,'Color','w');
p = piechart(LT52_in_optic_Thr, "Synapse", "Type", 'LabelStyle','name');
title(sprintf('LT52 Inputs from Optic lobes (synapse > %d)', in_optic_Thr));
p.ColorOrder = pieColors;

%% out_central LT52
out_central_Thr = 300; % <- 직접 조정
LT52_out_central_Thr = LT52_out_central;
LT52_out_central_Thr(LT52_out_central_Thr.Synapse < out_central_Thr, :) = [];

for i=1:size(LT52_out_central_Thr,1)
    current_types = LT52_out_central_Thr.Type{i};
    current_types_root_ids = FAFB_consolidated_cell_types.root_id(strcmp(FAFB_consolidated_cell_types.primary_type,current_types));
    superclasses = FAFB_classfication.super_class(ismember(FAFB_classfication.root_id,current_types_root_ids));
    [u,~,idx] = unique(superclasses); counts = histcounts(idx,1:numel(u)+1);
    [~,maxIdx] = max(counts); LT52_out_central_Thr.superclass{i} = superclasses{maxIdx};
end
LT52_out_central_Thr = sortrows(LT52_out_central_Thr,"superclass","ascend");
LT52_out_central_Thr = sortrows(LT52_out_central_Thr,"Type","ascend");
LT52_out_central_Thr = sortrows(LT52_out_central_Thr,"superclass","ascend");
numSlices = height(LT52_out_central_Thr); pieColors = zeros(numSlices,3);
for i=1:numSlices
    color_index = matches(uniqueSuperclass, LT52_out_central_Thr.superclass{i});
    pieColors(i,:) = colors(color_index,:);
end
figure; set(gcf,'Color','w');
p = piechart(LT52_out_central_Thr, "Synapse", "Type", 'LabelStyle','name');
title(sprintf('LT52 Outputs to Central Neurons (synapse > %d)', out_central_Thr));
p.ColorOrder = pieColors;

%% out_optic LT52
out_optic_Thr = 400; % <- 직접 조정
LT52_out_optic_Thr = LT52_out_optic;
LT52_out_optic_Thr(LT52_out_optic_Thr.Synapse < out_optic_Thr, :) = [];

for i=1:size(LT52_out_optic_Thr,1)
    current_types = LT52_out_optic_Thr.Type{i};
    current_types_root_ids = FAFB_consolidated_cell_types.root_id(strcmp(FAFB_consolidated_cell_types.primary_type,current_types));
    superclasses = FAFB_classfication.super_class(ismember(FAFB_classfication.root_id,current_types_root_ids));
    [u,~,idx] = unique(superclasses); counts = histcounts(idx,1:numel(u)+1);
    [~,maxIdx] = max(counts); LT52_out_optic_Thr.superclass{i} = superclasses{maxIdx};
end
LT52_out_optic_Thr = sortrows(LT52_out_optic_Thr,"superclass","ascend");
LT52_out_optic_Thr = sortrows(LT52_out_optic_Thr,"Type","ascend");
LT52_out_optic_Thr = sortrows(LT52_out_optic_Thr,"superclass","ascend");
numSlices = height(LT52_out_optic_Thr); pieColors = zeros(numSlices,3);
for i=1:numSlices
    color_index = matches(uniqueSuperclass, LT52_out_optic_Thr.superclass{i});
    pieColors(i,:) = colors(color_index,:);
end
figure; set(gcf,'Color','w');
p = piechart(LT52_out_optic_Thr, "Synapse", "Type", 'LabelStyle','name');
title(sprintf('LT52 Outputs to Optic lobes (synapse > %d)', out_optic_Thr));
p.ColorOrder = pieColors;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% =========================
%% aMe1
%% =========================
aMe1_Total_in  = sum(Wantsee.InConnections{4}.syn_count);
aMe1_Total_out = sum(Wantsee.OutConnections{4}.syn_count);
aMe1_in_optic   = Wantsee.unique_InOptic_Types{4};
aMe1_in_optic   = cell2table(aMe1_in_optic,   'VariableNames', {'Type','Synapse','Number'});
aMe1_in_central = Wantsee.unique_InCentral_Types{4};
aMe1_in_central = cell2table(aMe1_in_central, 'VariableNames', {'Type','Synapse','Number'});
aMe1_out_optic  = Wantsee.unique_OutOptic_Types{4};
aMe1_out_optic  = cell2table(aMe1_out_optic,  'VariableNames', {'Type','Synapse','Number'});
aMe1_out_central= Wantsee.unique_OutCentral_Types{4};
aMe1_out_central= cell2table(aMe1_out_central,'VariableNames', {'Type','Synapse','Number'});

Indata = [sum(aMe1_in_central.Synapse), sum(aMe1_in_optic.Synapse)];
figure; set(gcf,'Color','w'); piechart(Indata, Names={'Central','Optic'}); title('aMe1 Input');

Outdata = [sum(aMe1_out_central.Synapse), sum(aMe1_out_optic.Synapse)];
figure; set(gcf,'Color','w'); piechart(Outdata, Names={'Central','Optic'}); title('aMe1 Output');
%% === (추가) aMe1: 전체 In/Out 타입별 piechart ===
aMe1_Total_in_tbl = Wantsee.InNeuronTypes{4};
aMe1_Total_in_tbl = cell2table(aMe1_Total_in_tbl(:,[1 2]), 'VariableNames', {'Type','Synapse'});

aMe1_Total_out_tbl = Wantsee.OutNeuronTypes{4};
aMe1_Total_out_tbl = cell2table(aMe1_Total_out_tbl(:,[1 2]), 'VariableNames', {'Type','Synapse'});

% 전체 Input 타입별 pie
in_Thr = 60; % <-- 네가 조정
aMe1_Total_in_Thr = aMe1_Total_in_tbl;
aMe1_Total_in_Thr(aMe1_Total_in_Thr.Synapse < in_Thr, :) = [];

for i=1:height(aMe1_Total_in_Thr)
    current_types = aMe1_Total_in_Thr.Type{i};
    current_types_root_ids = FAFB_consolidated_cell_types.root_id(strcmp(FAFB_consolidated_cell_types.primary_type,current_types));
    superclasses = FAFB_classfication.super_class(ismember(FAFB_classfication.root_id,current_types_root_ids));
    if ~isempty(superclasses)
        [u,~,idx] = unique(superclasses); counts = histcounts(idx,1:numel(u)+1);
        [~,maxIdx] = max(counts);
        aMe1_Total_in_Thr.superclass{i} = superclasses{maxIdx};
    else
        current_types_root_ids = str_to_int64_exact(current_types);   % int64로 정확 변환

        superclasses=FAFB_classfication.super_class(FAFB_classfication.root_id==current_types_root_ids);
        aMe1_Total_in_Thr.superclass{i}=superclasses{1};
    end
end

aMe1_Total_in_Thr = sortrows(aMe1_Total_in_Thr,"superclass","ascend");
aMe1_Total_in_Thr = sortrows(aMe1_Total_in_Thr,"Type","ascend");
aMe1_Total_in_Thr = sortrows(aMe1_Total_in_Thr,"superclass","ascend");

numSlices = height(aMe1_Total_in_Thr);
pieColors = zeros(numSlices,3);
for i=1:numSlices
    color_index = matches(uniqueSuperclass, aMe1_Total_in_Thr.superclass{i});
    pieColors(i,:) = colors(color_index,:);
end
figure; set(gcf,'Color','w')
p = piechart(aMe1_Total_in_Thr, "Synapse", "Type", 'LabelStyle','name');
title(sprintf('aMe1 Inputs Total (synapse > %d)', in_Thr));
p.ColorOrder = pieColors;

% 전체 Output 타입별 pie
out_Thr = 50; % <-- 네가 조정
aMe1_Total_out_Thr = aMe1_Total_out_tbl;
aMe1_Total_out_Thr(aMe1_Total_out_Thr.Synapse < out_Thr, :) = [];

for i=1:height(aMe1_Total_out_Thr)
    current_types = aMe1_Total_out_Thr.Type{i};
    current_types_root_ids = FAFB_consolidated_cell_types.root_id(strcmp(FAFB_consolidated_cell_types.primary_type,current_types));
    superclasses = FAFB_classfication.super_class(ismember(FAFB_classfication.root_id,current_types_root_ids));
    if ~isempty(superclasses)
        [u,~,idx] = unique(superclasses); counts = histcounts(idx,1:numel(u)+1);
        [~,maxIdx] = max(counts);
        aMe1_Total_out_Thr.superclass{i} = superclasses{maxIdx};
    else
        current_types_root_ids = str_to_int64_exact(current_types);   % int64로 정확 변환
        superclasses=FAFB_classfication.super_class(FAFB_classfication.root_id==current_types_root_ids);
        aMe1_Total_out_Thr.superclass{i}=superclasses{1};
    end
end

aMe1_Total_out_Thr = sortrows(aMe1_Total_out_Thr,"superclass","ascend");
aMe1_Total_out_Thr = sortrows(aMe1_Total_out_Thr,"Type","ascend");
aMe1_Total_out_Thr = sortrows(aMe1_Total_out_Thr,"superclass","ascend");

numSlices = height(aMe1_Total_out_Thr);
pieColors = zeros(numSlices,3);
for i=1:numSlices
    color_index = matches(uniqueSuperclass, aMe1_Total_out_Thr.superclass{i});
    pieColors(i,:) = colors(color_index,:);
end
figure; set(gcf,'Color','w')
p = piechart(aMe1_Total_out_Thr, "Synapse", "Type", 'LabelStyle','name');
title(sprintf('aMe1 Outputs Total (synapse > %d)', out_Thr));
p.ColorOrder = pieColors;

%% in_central aMe1
in_central_Thr = 30; % <- 직접 조정
aMe1_in_central_Thr = aMe1_in_central;
aMe1_in_central_Thr(aMe1_in_central_Thr.Synapse < in_central_Thr, :) = [];

for i=1:size(aMe1_in_central_Thr,1)
    current_types = aMe1_in_central_Thr.Type{i};
    current_types_root_ids = FAFB_consolidated_cell_types.root_id(strcmp(FAFB_consolidated_cell_types.primary_type,current_types));
    superclasses = FAFB_classfication.super_class(ismember(FAFB_classfication.root_id,current_types_root_ids));
    [u,~,idx] = unique(superclasses); counts = histcounts(idx,1:numel(u)+1);
    [~,maxIdx] = max(counts); aMe1_in_central_Thr.superclass{i} = superclasses{maxIdx};
end
aMe1_in_central_Thr = sortrows(aMe1_in_central_Thr,"superclass","ascend");
aMe1_in_central_Thr = sortrows(aMe1_in_central_Thr,"Type","ascend");
aMe1_in_central_Thr = sortrows(aMe1_in_central_Thr,"superclass","ascend");
numSlices = height(aMe1_in_central_Thr); pieColors = zeros(numSlices,3);
for i=1:numSlices
    color_index = matches(uniqueSuperclass, aMe1_in_central_Thr.superclass{i});
    pieColors(i,:) = colors(color_index,:);
end
figure; set(gcf,'Color','w');
p = piechart(aMe1_in_central_Thr, "Synapse", "Type", 'LabelStyle','name');
title(sprintf('aMe1 Inputs from Central brains (synapse > %d)', in_central_Thr));
p.ColorOrder = pieColors;

%% in_optic aMe1
in_optic_Thr =40; % <- 직접 조정
aMe1_in_optic_Thr = aMe1_in_optic;
aMe1_in_optic_Thr(aMe1_in_optic_Thr.Synapse < in_optic_Thr, :) = [];

for i=1:size(aMe1_in_optic_Thr,1)
    current_types = aMe1_in_optic_Thr.Type{i};
    current_types_root_ids = FAFB_consolidated_cell_types.root_id(strcmp(FAFB_consolidated_cell_types.primary_type,current_types));
    superclasses = FAFB_classfication.super_class(ismember(FAFB_classfication.root_id,current_types_root_ids));
    [u,~,idx] = unique(superclasses); counts = histcounts(idx,1:numel(u)+1);
    [~,maxIdx] = max(counts); aMe1_in_optic_Thr.superclass{i} = superclasses{maxIdx};
end
aMe1_in_optic_Thr = sortrows(aMe1_in_optic_Thr,"superclass","ascend");
aMe1_in_optic_Thr = sortrows(aMe1_in_optic_Thr,"Type","ascend");
aMe1_in_optic_Thr = sortrows(aMe1_in_optic_Thr,"superclass","ascend");
numSlices = height(aMe1_in_optic_Thr); pieColors = zeros(numSlices,3);
for i=1:numSlices
    color_index = matches(uniqueSuperclass, aMe1_in_optic_Thr.superclass{i});
    pieColors(i,:) = colors(color_index,:);
end
figure; set(gcf,'Color','w');
p = piechart(aMe1_in_optic_Thr, "Synapse", "Type", 'LabelStyle','name');
title(sprintf('aMe1 Inputs from Optic lobes (synapse > %d)', in_optic_Thr));
p.ColorOrder = pieColors;

%% out_central aMe1
out_central_Thr = 20; % <- 직접 조정
aMe1_out_central_Thr = aMe1_out_central;
aMe1_out_central_Thr(aMe1_out_central_Thr.Synapse < out_central_Thr, :) = [];

for i=1:size(aMe1_out_central_Thr,1)
    current_types = aMe1_out_central_Thr.Type{i};
    current_types_root_ids = FAFB_consolidated_cell_types.root_id(strcmp(FAFB_consolidated_cell_types.primary_type,current_types));
    superclasses = FAFB_classfication.super_class(ismember(FAFB_classfication.root_id,current_types_root_ids));
    [u,~,idx] = unique(superclasses); counts = histcounts(idx,1:numel(u)+1);
    [~,maxIdx] = max(counts); aMe1_out_central_Thr.superclass{i} = superclasses{maxIdx};
end
aMe1_out_central_Thr = sortrows(aMe1_out_central_Thr,"superclass","ascend");
aMe1_out_central_Thr = sortrows(aMe1_out_central_Thr,"Type","ascend");
aMe1_out_central_Thr = sortrows(aMe1_out_central_Thr,"superclass","ascend");
numSlices = height(aMe1_out_central_Thr); pieColors = zeros(numSlices,3);
for i=1:numSlices
    color_index = matches(uniqueSuperclass, aMe1_out_central_Thr.superclass{i});
    pieColors(i,:) = colors(color_index,:);
end
figure; set(gcf,'Color','w');
p = piechart(aMe1_out_central_Thr, "Synapse", "Type", 'LabelStyle','name');
title(sprintf('aMe1 Outputs to Central Neurons (synapse > %d)', out_central_Thr));
p.ColorOrder = pieColors;

%% out_optic aMe1
out_optic_Thr = 45; % <- 직접 조정
aMe1_out_optic_Thr = aMe1_out_optic;
aMe1_out_optic_Thr(aMe1_out_optic_Thr.Synapse < out_optic_Thr, :) = [];

for i=1:size(aMe1_out_optic_Thr,1)
    current_types = aMe1_out_optic_Thr.Type{i};
    current_types_root_ids = FAFB_consolidated_cell_types.root_id(strcmp(FAFB_consolidated_cell_types.primary_type,current_types));
    superclasses = FAFB_classfication.super_class(ismember(FAFB_classfication.root_id,current_types_root_ids));
    [u,~,idx] = unique(superclasses); counts = histcounts(idx,1:numel(u)+1);
    [~,maxIdx] = max(counts); aMe1_out_optic_Thr.superclass{i} = superclasses{maxIdx};
end
aMe1_out_optic_Thr = sortrows(aMe1_out_optic_Thr,"superclass","ascend");
aMe1_out_optic_Thr = sortrows(aMe1_out_optic_Thr,"Type","ascend");
aMe1_out_optic_Thr = sortrows(aMe1_out_optic_Thr,"superclass","ascend");
numSlices = height(aMe1_out_optic_Thr); pieColors = zeros(numSlices,3);
for i=1:numSlices
    color_index = matches(uniqueSuperclass, aMe1_out_optic_Thr.superclass{i});
    pieColors(i,:) = colors(color_index,:);
end
figure; set(gcf,'Color','w');
p = piechart(aMe1_out_optic_Thr, "Synapse", "Type", 'LabelStyle','name');
title(sprintf('aMe1 Outputs to Optic lobes (synapse > %d)', out_optic_Thr));
p.ColorOrder = pieColors;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% =========================
%% cLM01
%% =========================
cLM01_Total_in  = sum(Wantsee.InConnections{5}.syn_count);
cLM01_Total_out = sum(Wantsee.OutConnections{5}.syn_count);
cLM01_in_optic   = Wantsee.unique_InOptic_Types{5};
cLM01_in_optic   = cell2table(cLM01_in_optic,   'VariableNames', {'Type','Synapse','Number'});
cLM01_in_central = Wantsee.unique_InCentral_Types{5};
cLM01_in_central = cell2table(cLM01_in_central, 'VariableNames', {'Type','Synapse','Number'});
cLM01_out_optic  = Wantsee.unique_OutOptic_Types{5};
cLM01_out_optic  = cell2table(cLM01_out_optic,  'VariableNames', {'Type','Synapse','Number'});
cLM01_out_central= Wantsee.unique_OutCentral_Types{5};
cLM01_out_central= cell2table(cLM01_out_central,'VariableNames', {'Type','Synapse','Number'});

Indata = [sum(cLM01_in_central.Synapse), sum(cLM01_in_optic.Synapse)];
figure; set(gcf,'Color','w'); piechart(Indata, Names={'Central','Optic'}); title('cLM01 Input');

Outdata = [sum(cLM01_out_central.Synapse), sum(cLM01_out_optic.Synapse)];
figure; set(gcf,'Color','w'); piechart(Outdata, Names={'Central','Optic'}); title('cLM01 Output');


%% === (추가) cLM01: 전체 In/Out 타입별 piechart ===
cLM01_Total_in_tbl = Wantsee.InNeuronTypes{5};
cLM01_Total_in_tbl = cell2table(cLM01_Total_in_tbl(:,[1 2]), 'VariableNames', {'Type','Synapse'});

cLM01_Total_out_tbl = Wantsee.OutNeuronTypes{5};
cLM01_Total_out_tbl = cell2table(cLM01_Total_out_tbl(:,[1 2]), 'VariableNames', {'Type','Synapse'});

% 전체 Input 타입별 pie
in_Thr = 15; % <-- 네가 조정
cLM01_Total_in_Thr = cLM01_Total_in_tbl;
cLM01_Total_in_Thr(cLM01_Total_in_Thr.Synapse < in_Thr, :) = [];

for i=1:height(cLM01_Total_in_Thr)
    current_types = cLM01_Total_in_Thr.Type{i};
    current_types_root_ids = FAFB_consolidated_cell_types.root_id(strcmp(FAFB_consolidated_cell_types.primary_type,current_types));
    superclasses = FAFB_classfication.super_class(ismember(FAFB_classfication.root_id,current_types_root_ids));
    [u,~,idx] = unique(superclasses); counts = histcounts(idx,1:numel(u)+1);
    [~,maxIdx] = max(counts);
    cLM01_Total_in_Thr.superclass{i} = superclasses{maxIdx};
end

cLM01_Total_in_Thr = sortrows(cLM01_Total_in_Thr,"superclass","ascend");
cLM01_Total_in_Thr = sortrows(cLM01_Total_in_Thr,"Type","ascend");
cLM01_Total_in_Thr = sortrows(cLM01_Total_in_Thr,"superclass","ascend");

numSlices = height(cLM01_Total_in_Thr);
pieColors = zeros(numSlices,3);
for i=1:numSlices
    color_index = matches(uniqueSuperclass, cLM01_Total_in_Thr.superclass{i});
    pieColors(i,:) = colors(color_index,:);
end
figure; set(gcf,'Color','w')
p = piechart(cLM01_Total_in_Thr, "Synapse", "Type", 'LabelStyle','name');
title(sprintf('cLM01 Inputs Total (synapse > %d)', in_Thr));
p.ColorOrder = pieColors;

% 전체 Output 타입별 pie
out_Thr = 35; % <-- 네가 조정
cLM01_Total_out_Thr = cLM01_Total_out_tbl;
cLM01_Total_out_Thr(cLM01_Total_out_Thr.Synapse < out_Thr, :) = [];

for i=1:height(cLM01_Total_out_Thr)
    current_types = cLM01_Total_out_Thr.Type{i};
    current_types_root_ids = FAFB_consolidated_cell_types.root_id(strcmp(FAFB_consolidated_cell_types.primary_type,current_types));
    superclasses = FAFB_classfication.super_class(ismember(FAFB_classfication.root_id,current_types_root_ids));
    [u,~,idx] = unique(superclasses); counts = histcounts(idx,1:numel(u)+1);
    [~,maxIdx] = max(counts);
    cLM01_Total_out_Thr.superclass{i} = superclasses{maxIdx};
end

cLM01_Total_out_Thr = sortrows(cLM01_Total_out_Thr,"superclass","ascend");
cLM01_Total_out_Thr = sortrows(cLM01_Total_out_Thr,"Type","ascend");
cLM01_Total_out_Thr = sortrows(cLM01_Total_out_Thr,"superclass","ascend");

numSlices = height(cLM01_Total_out_Thr);
pieColors = zeros(numSlices,3);
for i=1:numSlices
    color_index = matches(uniqueSuperclass, cLM01_Total_out_Thr.superclass{i});
    pieColors(i,:) = colors(color_index,:);
end
figure; set(gcf,'Color','w')
p = piechart(cLM01_Total_out_Thr, "Synapse", "Type", 'LabelStyle','name');
title(sprintf('cLM01 Outputs Total (synapse > %d)', out_Thr));
p.ColorOrder = pieColors;

%% in_central cLM01
in_central_Thr = 6; % <- 직접 조정
cLM01_in_central_Thr = cLM01_in_central;
cLM01_in_central_Thr(cLM01_in_central_Thr.Synapse < in_central_Thr, :) = [];

for i=1:size(cLM01_in_central_Thr,1)
    current_types = cLM01_in_central_Thr.Type{i};
    current_types_root_ids = FAFB_consolidated_cell_types.root_id(strcmp(FAFB_consolidated_cell_types.primary_type,current_types));
    superclasses = FAFB_classfication.super_class(ismember(FAFB_classfication.root_id,current_types_root_ids));
    [u,~,idx] = unique(superclasses); counts = histcounts(idx,1:numel(u)+1);
    [~,maxIdx] = max(counts); cLM01_in_central_Thr.superclass{i} = superclasses{maxIdx};
end
cLM01_in_central_Thr = sortrows(cLM01_in_central_Thr,"superclass","ascend");
cLM01_in_central_Thr = sortrows(cLM01_in_central_Thr,"Type","ascend");
cLM01_in_central_Thr = sortrows(cLM01_in_central_Thr,"superclass","ascend");
numSlices = height(cLM01_in_central_Thr); pieColors = zeros(numSlices,3);
for i=1:numSlices
    color_index = matches(uniqueSuperclass, cLM01_in_central_Thr.superclass{i});
    pieColors(i,:) = colors(color_index,:);
end
figure; set(gcf,'Color','w');
p = piechart(cLM01_in_central_Thr, "Synapse", "Type", 'LabelStyle','name');
title(sprintf('cLM01 Inputs from Central brains (synapse > %d)', in_central_Thr));
p.ColorOrder = pieColors;

%% in_optic cLM01
in_optic_Thr = 15; % <- 직접 조정
cLM01_in_optic_Thr = cLM01_in_optic;
cLM01_in_optic_Thr(cLM01_in_optic_Thr.Synapse < in_optic_Thr, :) = [];

for i=1:size(cLM01_in_optic_Thr,1)
    current_types = cLM01_in_optic_Thr.Type{i};
    current_types_root_ids = FAFB_consolidated_cell_types.root_id(strcmp(FAFB_consolidated_cell_types.primary_type,current_types));
    superclasses = FAFB_classfication.super_class(ismember(FAFB_classfication.root_id,current_types_root_ids));
    [u,~,idx] = unique(superclasses); counts = histcounts(idx,1:numel(u)+1);
    [~,maxIdx] = max(counts); cLM01_in_optic_Thr.superclass{i} = superclasses{maxIdx};
end
cLM01_in_optic_Thr = sortrows(cLM01_in_optic_Thr,"superclass","ascend");
cLM01_in_optic_Thr = sortrows(cLM01_in_optic_Thr,"Type","ascend");
cLM01_in_optic_Thr = sortrows(cLM01_in_optic_Thr,"superclass","ascend");
numSlices = height(cLM01_in_optic_Thr); pieColors = zeros(numSlices,3);
for i=1:numSlices
    color_index = matches(uniqueSuperclass, cLM01_in_optic_Thr.superclass{i});
    pieColors(i,:) = colors(color_index,:);
end
figure; set(gcf,'Color','w');
p = piechart(cLM01_in_optic_Thr, "Synapse", "Type", 'LabelStyle','name');
title(sprintf('cLM01 Inputs from Optic lobes (synapse > %d)', in_optic_Thr));
p.ColorOrder = pieColors;

%% out_central cLM01
out_central_Thr = 8; % <- 직접 조정
cLM01_out_central_Thr = cLM01_out_central;
cLM01_out_central_Thr(cLM01_out_central_Thr.Synapse < out_central_Thr, :) = [];

for i=1:size(cLM01_out_central_Thr,1)
    current_types = cLM01_out_central_Thr.Type{i};
    current_types_root_ids = FAFB_consolidated_cell_types.root_id(strcmp(FAFB_consolidated_cell_types.primary_type,current_types));
    superclasses = FAFB_classfication.super_class(ismember(FAFB_classfication.root_id,current_types_root_ids));
    [u,~,idx] = unique(superclasses); counts = histcounts(idx,1:numel(u)+1);
    [~,maxIdx] = max(counts); cLM01_out_central_Thr.superclass{i} = superclasses{maxIdx};
end
cLM01_out_central_Thr = sortrows(cLM01_out_central_Thr,"superclass","ascend");
cLM01_out_central_Thr = sortrows(cLM01_out_central_Thr,"Type","ascend");
cLM01_out_central_Thr = sortrows(cLM01_out_central_Thr,"superclass","ascend");
numSlices = height(cLM01_out_central_Thr); pieColors = zeros(numSlices,3);
for i=1:numSlices
    color_index = matches(uniqueSuperclass, cLM01_out_central_Thr.superclass{i});
    pieColors(i,:) = colors(color_index,:);
end
figure; set(gcf,'Color','w');
p = piechart(cLM01_out_central_Thr, "Synapse", "Type", 'LabelStyle','name');
title(sprintf('cLM01 Outputs to Central Neurons (synapse > %d)', out_central_Thr));
p.ColorOrder = pieColors;

%% out_optic cLM01
out_optic_Thr = 40; % <- 직접 조정
cLM01_out_optic_Thr = cLM01_out_optic;
cLM01_out_optic_Thr(cLM01_out_optic_Thr.Synapse < out_optic_Thr, :) = [];

for i=1:size(cLM01_out_optic_Thr,1)
    current_types = cLM01_out_optic_Thr.Type{i};
    current_types_root_ids = FAFB_consolidated_cell_types.root_id(strcmp(FAFB_consolidated_cell_types.primary_type,current_types));
    superclasses = FAFB_classfication.super_class(ismember(FAFB_classfication.root_id,current_types_root_ids));
    [u,~,idx] = unique(superclasses); counts = histcounts(idx,1:numel(u)+1);
    [~,maxIdx] = max(counts); cLM01_out_optic_Thr.superclass{i} = superclasses{maxIdx};
end
cLM01_out_optic_Thr = sortrows(cLM01_out_optic_Thr,"superclass","ascend");
cLM01_out_optic_Thr = sortrows(cLM01_out_optic_Thr,"Type","ascend");
cLM01_out_optic_Thr = sortrows(cLM01_out_optic_Thr,"superclass","ascend");
numSlices = height(cLM01_out_optic_Thr); pieColors = zeros(numSlices,3);
for i=1:numSlices
    color_index = matches(uniqueSuperclass, cLM01_out_optic_Thr.superclass{i});
    pieColors(i,:) = colors(color_index,:);
end
figure; set(gcf,'Color','w');
p = piechart(cLM01_out_optic_Thr, "Synapse", "Type", 'LabelStyle','name');
title(sprintf('cLM01 Outputs to Optic lobes (synapse > %d)', out_optic_Thr));
p.ColorOrder = pieColors;


function n = str_to_int64_exact(s)
% 정확한 int64 파서 (double 미사용)
% 사용법: n = str_to_int64_exact('720575940624333706')

    % ---- 입력 정리 ----
    if isstring(s), s = char(s); end
    s = strtrim(s);
    if isempty(s), error('빈 문자열입니다.'); end

    % 부호 처리 (여기 예시는 양수지만 일반 처리)
    sign = 1; i = 1;
    if s(1) == '+' || s(1) == '-'
        if s(1) == '-', sign = -1; end
        i = 2;
        if i > numel(s), error('숫자가 없습니다.'); end
    end

    % ---- int64 누적 계산 ----
    n   = int64(0);
    ten = int64(10);
    z   = int64('0');

    for k = i:numel(s)
        d = int64(s(k)) - z;            % 현재 자리의 숫자 (int64)
        if d < 0 || d > 9
            error('숫자가 아닌 문자가 %d번째에 있습니다.', k);
        end
        n = n * ten + d;               % 모두 int64 연산 → 정밀도 손실 없음
    end

    if sign < 0, n = -n; end
end
