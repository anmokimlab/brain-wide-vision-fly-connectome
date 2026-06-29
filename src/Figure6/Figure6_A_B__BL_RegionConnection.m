clear all; close all; clc;

load('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\Figure1_Overall\Bilateral_Neurons_Thr0.mat')

%% Codex 데이터 가져오기
opt = detectImportOptions('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\neurons.csv');
opt = setvartype(opt,'root_id','int64');
FAFBNeurons = readtable('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\neurons.csv',opt);

opt = detectImportOptions('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\connections_no_threshold.csv');
opt = setvartype(opt,'pre_root_id','int64');
opt = setvartype(opt,'post_root_id','int64');
FAFBConnections = readtable('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\connections_no_threshold.csv',opt);

opt = detectImportOptions('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\consolidated_cell_types.csv');
opt = setvartype(opt,'root_id','int64');
FAFBTypes = readtable('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\consolidated_cell_types.csv',opt);

% load('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_v2\CodexData\FAFBConnectome_my.mat')
NT_Types=unique(FAFBNeurons.nt_type);
NT_Types(1,:)=[];
NT_Types{end+1,1}='Unidentified';

FAFBNeuropils=unique(FAFBConnections.neuropil);
FAFBNeuropil_OpticLobeRight={'LA_R';'AME_R';'ME_R';'LO_R';'LOP_R'};
FAFBNeuropil_OpticLobeLeft={'LA_L';'AME_L';'ME_L';'LO_L';'LOP_L'};
FAFBNeuropil_Central=FAFBNeuropils(~ismember(FAFBNeuropils,FAFBNeuropil_OpticLobeRight));
FAFBNeuropil_Central=FAFBNeuropil_Central(~ismember(FAFBNeuropil_Central,FAFBNeuropil_OpticLobeLeft));
FAFBNeuropil_Central=FAFBNeuropil_Central(~ismember(FAFBNeuropil_Central,'UNASGD'));
% FAFBNeuropil_Central = strrep(FAFBNeuropil_Central,'_',' ');
%%

for i=1:1:size(Bilateral_R_NPIs,1)
    idx=FAFBNeurons.root_id==Bilateral_R_NPIs.root_id(i);
    if ~any(idx)
        Bilateral_R_NPIs.nt_type{i}='Unidentified';
    elseif isempty(FAFBNeurons.nt_type{idx})
        Bilateral_R_NPIs.nt_type{i}='Unidentified';
    else
        Bilateral_R_NPIs.nt_type{i}=FAFBNeurons.nt_type{idx};
    end
end

for i=1:1:size(NT_Types,1)
    idx=strcmp(Bilateral_R_NPIs.nt_type,NT_Types{i,1});
    NT_Types{i,2}=sum(idx);
end
%%  From Optic lobe to Optic lobe
BL_Connection=cell(length(FAFBNeuropils),length(FAFBNeuropils));

for i=1:1:size(Bilateral_R_NPIs,1)
    current_root_id=Bilateral_R_NPIs.root_id(i);
    PostConnections=FAFBConnections(FAFBConnections.post_root_id==current_root_id,:);
    PreConnections=FAFBConnections(FAFBConnections.pre_root_id==current_root_id,:);
    [OutNeuropils, ~, ic]=unique(PreConnections.neuropil);
    for j=1:1:size(OutNeuropils,1)
        idx=ic==j;
        OutNeuropils{j,2}=sum(PreConnections.syn_count(idx));
    end
    [InNeuropils, ~, ic]=unique(PostConnections.neuropil);
    for j=1:1:size(InNeuropils,1)
        idx=ic==j;
        InNeuropils{j,2}=sum(PostConnections.syn_count(idx));
    end

    for j=1:1:size(InNeuropils,1)
        idx_In=find(strcmp(FAFBNeuropils,InNeuropils{j,1}));
        syn_ratio_In=InNeuropils{j,2}/sum(cell2mat(InNeuropils(:,2)));
        for k=1:1:size(OutNeuropils,1)
            idx_Out=find(strcmp(FAFBNeuropils,OutNeuropils{k,1}));

            temp=BL_Connection{idx_In,idx_Out};
            temp{end+1,1}=current_root_id;
            temp{end,2}=OutNeuropils{k,2}*syn_ratio_In;
            temp{end,3}=Bilateral_R_NPIs.nt_type{i};
            BL_Connection{idx_In,idx_Out}=temp;

        end
    end
end


%% AME 를 ME에 합치기
% RightFB._NPIs.root_id가 유니크하다는 걸 전제로
aMe_R_idx=find(strcmp(FAFBNeuropils,'AME_R'));
Me_R_idx=find(strcmp(FAFBNeuropils,'ME_R'));
aMe_L_idx=find(strcmp(FAFBNeuropils,'AME_L'));
Me_L_idx=find(strcmp(FAFBNeuropils,'ME_L'));

for i=1:1:size(BL_Connection,2)
    Me_R_temp=BL_Connection{Me_R_idx,i};
    AMe_R_Temp=BL_Connection{aMe_R_idx,i};
    if isempty(Me_R_temp)||isempty(AMe_R_Temp)
        BL_Connection{Me_R_idx,i}=[Me_R_temp;AMe_R_Temp];
    else
        [MergedAMeIdx,MergedMeIdx]=(ismember(cell2mat(AMe_R_Temp(:,1)),cell2mat(Me_R_temp(:,1))));
        MergedAMeIdx=find(MergedAMeIdx);
        MergedMeIdx(MergedMeIdx==0)=[];
        if ~isempty(MergedMeIdx)
            for j=1:1:size(MergedMeIdx,1)
                Me_R_temp{MergedMeIdx(j),2}= Me_R_temp{MergedMeIdx(j),2}+AMe_R_Temp{MergedAMeIdx(j),2};
            end
        end
        % --- [수정] AMe_R에만 있는 고유 뉴런 찾아서 합치기 ---
        % 'ismember'를 다시 사용하여 Me_R에 '없는' 뉴런을 찾습니다.
        NonMergedAMeIdx = ~ismember(cell2mat(AMe_R_Temp(:,1)), cell2mat(Me_R_temp(:,1)));

        % 그 고유 뉴런 데이터(NonMergedAMeIdx)를 Me_R_temp의 맨 아래에 붙입니다.
        Me_R_temp = [Me_R_temp; AMe_R_Temp(NonMergedAMeIdx, :)];
        % --- 수정 끝 ---

        % 이제 'Me_R_temp'는 [뉴런A, 10], [뉴런B, 25], [뉴런C, 30]을 모두 가집니다.
        BL_Connection{Me_R_idx,i}=Me_R_temp;
    end
end

for i=1:1:size(BL_Connection,1)
    Me_L_temp=BL_Connection{i,Me_L_idx};
    AMe_L_Temp=BL_Connection{i,aMe_L_idx};
    if isempty(Me_L_temp)||isempty(AMe_L_Temp)
        BL_Connection{i,Me_L_idx}=[Me_L_temp;AMe_L_Temp];
    else
        [MergedAMeIdx,MergedMeIdx]=(ismember(cell2mat(AMe_L_Temp(:,1)),cell2mat(Me_L_temp(:,1))));
        MergedAMeIdx=find(MergedAMeIdx);
        MergedMeIdx(MergedMeIdx==0)=[];
        if ~isempty(MergedMeIdx)
            for j=1:1:size(MergedMeIdx,1)
                Me_L_temp{MergedMeIdx(j),2}= Me_L_temp{MergedMeIdx(j),2}+AMe_L_Temp{MergedAMeIdx(j),2}; % 비율이니까 그냥 더해도 됨
            end
        end
        % --- [수정] AMe_R에만 있는 고유 뉴런 찾아서 합치기 ---
        % 'ismember'를 다시 사용하여 Me_R에 '없는' 뉴런을 찾습니다.
        NonMergedAMeIdx = ~ismember(cell2mat(AMe_L_Temp(:,1)), cell2mat(Me_L_temp(:,1)));

        % 그 고유 뉴런 데이터(NonMergedAMeIdx)를 Me_L_temp의 맨 아래에 붙입니다.
        Me_L_temp = [Me_L_temp; AMe_L_Temp(NonMergedAMeIdx, :)];
        % --- 수정 끝 ---

        % 이제 'Me_L_temp'는 [뉴런A, 10], [뉴런B, 25], [뉴런C, 30]을 모두 가집니다.
        BL_Connection{i,Me_L_idx}=Me_L_temp;
    end
end
%% 버블
%%% Me R 55 Lo R 45 Lop R 43
idx_RightOpticLobes=[55 45 43];
idx_LeftOpticLobes=[54 44 42];

BL_Connection_Matrix=zeros(length(idx_RightOpticLobes),length(idx_LeftOpticLobes));

for i=1:1:length(idx_RightOpticLobes)
    for j=1:1:length(idx_LeftOpticLobes)
        currentConnection=BL_Connection{idx_RightOpticLobes(i),idx_LeftOpticLobes(j)};
        if isempty(currentConnection)
            BL_Connection_Matrix(i,j)=0;
        else
            BL_Connection_Matrix(i,j)=size(currentConnection,1);
        end
    end
end


figure(1);set(gcf,'Color','w')
NeuronN_Thr=0;
temp_matrix=(BL_Connection_Matrix);
[YData, XData]=find(temp_matrix>=NeuronN_Thr);
Bino_Connection_Matrix_bubble=temp_matrix(temp_matrix>=NeuronN_Thr);
XData=XData';
YData=YData';
Bino_Connection_Matrix_bubble=Bino_Connection_Matrix_bubble';
bubblechart(XData,YData,Bino_Connection_Matrix_bubble);
bubblelegend('Neuron','Location','eastoutside')
set(gca,'XTick',1:1:size(BL_Connection_Matrix,2),'XTickLabels',{'Me R','Lo R','LoP R'},...
    'YTick',1:1:size(BL_Connection_Matrix,1),'YTickLabels',{'Me R','Lo R','LoP R'},'Box','off','TickDir','out')
title('Right Bino Neuron')
%% sending neuron names...
temp=sum(BL_Connection_Matrix,2);
[~,SortingIndex_Neuron]=sort(temp,'descend');

%%% Me R 55 Lo R 45 Lop R 43
idx_RightOpticLobes=[55 45 43];
idx_LeftOpticLobes=[54 44 42];
BL_Connection_NeuronsNames=cell(length(idx_RightOpticLobes),length(idx_LeftOpticLobes));
BL_Connection_NeuronsNamesType=cell(length(idx_RightOpticLobes),length(idx_LeftOpticLobes));

for i=1:1:length(idx_RightOpticLobes)
    for j=1:1:length(idx_LeftOpticLobes)
        currentConnection=BL_Connection{idx_RightOpticLobes(i),idx_LeftOpticLobes(j)};
        if isempty(currentConnection)
            continue;
        else
            for k=1:1:size(currentConnection,1)
                idx_Types=FAFBTypes.root_id==currentConnection{k,1};
                currentConnection{k,4}=FAFBTypes.primary_type{idx_Types};
            end
            [uniquecurrentConnectionType,~,ic] = unique(currentConnection(:,4));
            a_counts = accumarray(ic,1);
            uniquecurrentConnectionType(:,2)=num2cell(a_counts);
            for k=1:1:size(uniquecurrentConnectionType,1)
                idx=ic==k;
                uniquecurrentConnectionType{k,3}=sum(cell2mat(currentConnection(idx,2)));
            end
            uniquecurrentConnectionType=sortrows(uniquecurrentConnectionType,3,'descend');
            BL_Connection_NeuronsNames{i,j}=currentConnection;
            BL_Connection_NeuronsNamesType{i,j}=uniquecurrentConnectionType;
        end
    end
end
%%
ME_R_NumberNeurons=[BL_Connection_NeuronsNames{1, 1};BL_Connection_NeuronsNames{1, 2};BL_Connection_NeuronsNames{1, 3}];
ME_R_NumberNeurons=size(unique(cell2mat(ME_R_NumberNeurons(:,1))),1);
LO_R_NumberNeurons=[BL_Connection_NeuronsNames{2, 1};BL_Connection_NeuronsNames{2, 2};BL_Connection_NeuronsNames{2, 3}];
LO_R_NumberNeurons=size(unique(cell2mat(LO_R_NumberNeurons(:,1))),1);
LoP_R_NumberNeurons=[BL_Connection_NeuronsNames{3, 1};BL_Connection_NeuronsNames{3, 2};BL_Connection_NeuronsNames{3, 3}];
LoP_R_NumberNeurons=size(unique(cell2mat(LoP_R_NumberNeurons(:,1))),1);

ME_L_NumberNeurons=[BL_Connection_NeuronsNames{1, 1};BL_Connection_NeuronsNames{2, 1};BL_Connection_NeuronsNames{3, 1}];
ME_L_NumberNeurons=size(unique(cell2mat(ME_L_NumberNeurons(:,1))),1);
LO_L_NumberNeurons=[BL_Connection_NeuronsNames{1, 2};BL_Connection_NeuronsNames{2, 2};BL_Connection_NeuronsNames{3, 2}];
LO_L_NumberNeurons=size(unique(cell2mat(LO_L_NumberNeurons(:,1))),1);
LoP_L_NumberNeurons=[BL_Connection_NeuronsNames{1, 3};BL_Connection_NeuronsNames{2, 3};BL_Connection_NeuronsNames{3, 3}];
LoP_L_NumberNeurons=size(unique(cell2mat(LoP_L_NumberNeurons(:,1))),1);

BL_Connection_Matrix(4,1)=ME_L_NumberNeurons;
BL_Connection_Matrix(4,2)=LO_L_NumberNeurons;
BL_Connection_Matrix(4,3)=LoP_L_NumberNeurons;

BL_Connection_Matrix(1,4)=ME_R_NumberNeurons;
BL_Connection_Matrix(2,4)=LO_R_NumberNeurons;
BL_Connection_Matrix(3,4)=LoP_R_NumberNeurons;

figure(1);set(gcf,'Color','w')
NeuronN_Thr=0;
temp_matrix=(BL_Connection_Matrix);
[YData, XData]=find(temp_matrix>=NeuronN_Thr);
Bino_Connection_Matrix_bubble=temp_matrix(temp_matrix>=NeuronN_Thr);
XData=XData';
YData=YData';
Bino_Connection_Matrix_bubble=Bino_Connection_Matrix_bubble';
bubblechart(XData,YData,Bino_Connection_Matrix_bubble);
bubblelegend('Neuron','Location','eastoutside')
set(gca,'XTick',1:1:size(BL_Connection_Matrix,2),'XTickLabels',{'Me L','Lo L','LoP L'},...
    'YTick',1:1:size(BL_Connection_Matrix,1),'YTickLabels',{'Me R','Lo R','LoP R'},'Box','off','TickDir','out')
title('Right Bino Neuron')
print(gcf,'-depsc2','-vector','BL_Neuropil_dots.eps')

%% piechart로 해보기

NT_Me_R={};
NT_Me_L={};
NT_Lo_R={};
NT_Lo_L={};
NT_LoP_R={};
NT_LoP_L={};

for i=1:1:size(BL_Connection_NeuronsNames,1)
    for j=1:1:size(BL_Connection_NeuronsNames,2)
        if size(BL_Connection_NeuronsNames{i,j},1)<NeuronN_Thr
            BL_Connection_NeuronsNames{i,j}=[];
        end
    end
end

for i=1:1:size(BL_Connection_NeuronsNames,2)
    NT_Me_R=[NT_Me_R;BL_Connection_NeuronsNames{1,i}];
    NT_Lo_R=[NT_Lo_R;BL_Connection_NeuronsNames{2,i}];
    NT_LoP_R=[NT_LoP_R;BL_Connection_NeuronsNames{3,i}];
end

for i=1:1:size(BL_Connection_NeuronsNames,1)
    NT_Me_L=[NT_Me_L;BL_Connection_NeuronsNames{i,1}];
    NT_Lo_L=[NT_Lo_L;BL_Connection_NeuronsNames{i,2}];
    NT_LoP_L=[NT_LoP_L;BL_Connection_NeuronsNames{i,3}];
end
%%
%%%
NT_Me_R_unique=num2cell(unique(cell2mat(NT_Me_R(:,1))));

for i=1:1:size(NT_Me_R_unique,1)
    idx=find(cell2mat(NT_Me_R(:,1))==NT_Me_R_unique{i,1},1);
    NT_Me_R_unique{i,2}=NT_Me_R{idx,3};
end

[NT_Me_R,~,ic] = unique(NT_Me_R_unique(:,2));
a_counts = accumarray(ic,1);
NT_Me_R(:,2)=num2cell(a_counts);
%%%
NT_Lo_R_unique=num2cell(unique(cell2mat(NT_Lo_R(:,1))));

for i=1:1:size(NT_Lo_R_unique,1)
    idx=find(cell2mat(NT_Lo_R(:,1))==NT_Lo_R_unique{i,1},1);
    NT_Lo_R_unique{i,2}=NT_Lo_R{idx,3};
end

[NT_Lo_R,~,ic] = unique(NT_Lo_R_unique(:,2));
a_counts = accumarray(ic,1);
NT_Lo_R(:,2)=num2cell(a_counts);
%%%
NT_LoP_R_unique=num2cell(unique(cell2mat(NT_LoP_R(:,1))));

for i=1:1:size(NT_LoP_R_unique,1)
    idx=find(cell2mat(NT_LoP_R(:,1))==NT_LoP_R_unique{i,1},1);
    NT_LoP_R_unique{i,2}=NT_LoP_R{idx,3};
end

[NT_LoP_R,~,ic] = unique(NT_LoP_R_unique(:,2));
a_counts = accumarray(ic,1);
NT_LoP_R(:,2)=num2cell(a_counts);


%%%
NT_Me_L_unique=num2cell(unique(cell2mat(NT_Me_L(:,1))));

for i=1:1:size(NT_Me_L_unique,1)
    idx=find(cell2mat(NT_Me_L(:,1))==NT_Me_L_unique{i,1},1);
    NT_Me_L_unique{i,2}=NT_Me_L{idx,3};
end

[NT_Me_L,~,ic] = unique(NT_Me_L_unique(:,2));
a_counts = accumarray(ic,1);
NT_Me_L(:,2)=num2cell(a_counts);
%%%
NT_Lo_L_unique=num2cell(unique(cell2mat(NT_Lo_L(:,1))));

for i=1:1:size(NT_Lo_L_unique,1)
    idx=find(cell2mat(NT_Lo_L(:,1))==NT_Lo_L_unique{i,1},1);
    NT_Lo_L_unique{i,2}=NT_Lo_L{idx,3};
end

[NT_Lo_L,~,ic] = unique(NT_Lo_L_unique(:,2));
a_counts = accumarray(ic,1);
NT_Lo_L(:,2)=num2cell(a_counts);
%%%
NT_LoP_L_unique=num2cell(unique(cell2mat(NT_LoP_L(:,1))));

for i=1:1:size(NT_LoP_L_unique,1)
    idx=find(cell2mat(NT_LoP_L(:,1))==NT_LoP_L_unique{i,1},1);
    NT_LoP_L_unique{i,2}=NT_LoP_L{idx,3};
end

[NT_LoP_L,~,ic] = unique(NT_LoP_L_unique(:,2));
a_counts = accumarray(ic,1);
NT_LoP_L(:,2)=num2cell(a_counts);


%%

figure(2);set(gcf,'Color','w')

BubbleSize=[58 58 10 63 68 28 1];


b=bubblechart(1:1:length(BubbleSize),1:1:length(BubbleSize),BubbleSize);
bubblelegend('Neuron','Location','eastoutside')
% bubblesize([25*min(BubbleSize)/max(BubbleSize) 25])
bubblesize([5 25])

print(gcf,'-depsc2','-vector','BLRegionCircle.eps')

%%
figure(3);set(gcf,'Color','w')

x = ["LC14b" "LC14a1" "LC14a2" "MeMe e02" "MTe07" "MeMe e01" "MeMe e08" "MeMe e07" "l-LNv" "Others"];
y = [18 15 12 9 6 6 6 5 4 25];
bar(x,y)
set(gca,'Box','off','TickDir','out')
print(gcf,'-depsc2','-vector','BLTypeBar.eps')
