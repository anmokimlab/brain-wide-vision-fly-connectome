%% 데이터 가져오기
clear all; close all; clc
load('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\Figure1_Overall\Right_Neurons_Thr0.mat')
%%

[RightFB_type,~,ic]=unique(RightFB_NPIs.type);

RightFB_type=table(RightFB_type,'VariableNames',{'type'});

for i=1:1:size(RightFB_type,1)
    idx=ic==i;

    RightFB_type.root_id{i}=RightFB_NPIs.root_id(idx);
end
%%
% opt = detectImportOptions('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_v3\CodexData\connections_no_threshold.csv');
% opt = setvartype(opt,'pre_root_id','int64');
% opt = setvartype(opt,'post_root_id','int64');
% FAFBConnections = readtable('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_v3\CodexData\connections_no_threshold.csv',opt);

opt = detectImportOptions('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\connections_no_threshold.csv');
opt = setvartype(opt,'pre_root_id','int64');
opt = setvartype(opt,'post_root_id','int64');
FAFBConnections = readtable('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\connections_no_threshold.csv',opt);

opt = detectImportOptions('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\classification.csv');
opt = setvartype(opt,'root_id','int64');
FAFBClassification = readtable('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\classification.csv',opt);

opt = detectImportOptions('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\consolidated_cell_types.csv');
opt = setvartype(opt,'root_id','int64');
FAFBConsolidated_type = readtable('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\consolidated_cell_types.csv',opt);

FAFB_Superclass=unique(FAFBClassification.super_class);
%% FB의 타입별로 input neuron 들을 지정하기..
for i=1:1:size(RightFB_type,1)
    current_FB_root_id=RightFB_type.root_id{i};
    InConnections=FAFBConnections(ismember(FAFBConnections.post_root_id,current_FB_root_id),:);
    %%% 오직 central 에서의 입력만 보기
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    InConnections(ismember(InConnections.neuropil,{'LA_R';'AME_R';'ME_R';'LO_R';'LOP_R'}),:)=[];
    InConnections(ismember(InConnections.neuropil,{'LA_L';'AME_L';'ME_L';'LO_L';'LOP_L'}),:)=[];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% recursive connection 지우기
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % recursive_temp=RightFB_type.Recursive_neurons{i};
    % InConnections(ismember(InConnections.pre_root_id,unique(vertcat(recursive_temp{:}))),:)=[];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [Unique_inputs,~,ic]=unique(InConnections.pre_root_id);
    Unique_inputs=num2cell(Unique_inputs);
    for j=1:1:size(Unique_inputs,1)
        idx_syn=ic==j;
        Unique_inputs{j,2}=sum(InConnections.syn_count(idx_syn));
        idx_superclass=FAFBClassification.root_id==Unique_inputs{j,1};
        Unique_inputs{j,3}=FAFBClassification.super_class{idx_superclass};
    end
    In_SuperClass=cell2table(FAFB_Superclass,'VariableNames',{'superclass'});
    for j=1:1:size(In_SuperClass,1)
        idx_superclass=strcmp(Unique_inputs(:,3),In_SuperClass.superclass{j});
        In_SuperClass.syn_count{j}=sum(cell2mat(Unique_inputs(idx_superclass,2)));
        In_SuperClass.root_ids{j}=Unique_inputs(idx_superclass,[1 2]);
    end
    RightFB_type.input_superclass{i}=In_SuperClass;
  

    OutConnections=FAFBConnections(ismember(FAFBConnections.pre_root_id,current_FB_root_id),:);
    %%% recursive
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OutConnections(ismember(OutConnections.post_root_id,unique(vertcat(recursive_temp{:}))),:)=[];

    %%% 오직 central 에서의 입력만 보기
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    OutConnections(~ismember(OutConnections.neuropil,{'LA_R';'AME_R';'ME_R';'LO_R';'LOP_R'}),:)=[];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [Unique_outputs,~,ic]=unique(OutConnections.post_root_id);
    Unique_outputs=num2cell(Unique_outputs);
    for j=1:1:size(Unique_outputs,1)
        idx_syn=ic==j;
        Unique_outputs{j,2}=sum(OutConnections.syn_count(idx_syn));
        idx_superclass=FAFBClassification.root_id==Unique_outputs{j,1};
        Unique_outputs{j,3}=FAFBClassification.super_class{idx_superclass};
    end
    Out_SuperClass=cell2table(FAFB_Superclass,'VariableNames',{'superclass'});
    for j=1:1:size(Out_SuperClass,1)
        idx_superclass=strcmp(Unique_outputs(:,3),Out_SuperClass.superclass{j});
        Out_SuperClass.syn_count{j}=sum(cell2mat(Unique_outputs(idx_superclass,2)));
        Out_SuperClass.root_ids{j}=Unique_outputs(idx_superclass,[1 2]);
    end
    RightFB_type.output_superclass{i}=Out_SuperClass;
end

%% 그래프 그리기 direct 부터
% colors = [ 0.3010, 0.7450, 0.9330; %endocrine
%     0.2780, 0.6000, 0.8000; %motor
%     0.8500, 0.3250, 0.0980; %sensory
%     0.0000, 0.4470, 0.7410; %visual_projection
%     0.4660, 0.6740, 0.1880; %visual_centrifugalS
%     0.9290, 0.6940, 0.1250; % 'optic',...
%     0.4940, 0.1840, 0.5560;  % 'central',...
%     0.8500, 0.1500, 0.2000; % 'descending',...
%     0.6350, 0.5090, 0.2540 % 'ascending'
%     ];
colors=[ 0.6350, 0.5090, 0.2540;... %ascending
    0.4940, 0.1840, 0.5560;...%'central 
    0.8500, 0.1500, 0.2000;... %descending
    1.00,0.07,0.65;...%endocrine 
    0.8500, 0.3250, 0.0980; %motor 
    0.9290, 0.6940, 0.1250;... %optic
    0.3010, 0.7450, 0.9330;... %sensory 
    0.4660, 0.6740, 0.1880;... %visual_centifugal
    0.0000, 0.4470, 0.7410%visual_projection
    ];

%% ME, LO, LOP MULTI 로 나눠서 해보자
Me_FB_idx=[4 22 25 26 29 30 31 32 33 68 69 70 71 72 73 74 75 76 78 81 82 83 84 85 86 87 93];
Lo_FB_idx=[2 10 11 12 13 14 15 16 17 18 19 21 34 35 36 37 38 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 95];
Lop_FB_idx=[3 7 8 9 28 62 63 64 65 66];
Multi_FB_idx=[1 5 6 23 24 27 39 60 61 67 77 79 80 88 89 90 91 92 94 96];


Me_FB_input_superclass=zeros(size(FAFB_Superclass));
for i=1:1:size(Me_FB_idx,2)
    temp=cell2mat(RightFB_type.input_superclass{Me_FB_idx(i)}.syn_count);
    temp=temp/sum(temp);
    Me_FB_input_superclass=Me_FB_input_superclass +temp;
end
Me_FB_input_superclass=Me_FB_input_superclass./i;
Me_FB_input_superclass=[FAFB_Superclass num2cell(Me_FB_input_superclass)];

Lo_FB_input_superclass=zeros(size(FAFB_Superclass));
for i=1:1:size(Lo_FB_idx,2)
    temp=cell2mat(RightFB_type.input_superclass{Lo_FB_idx(i)}.syn_count);
    temp=temp/sum(temp);
    Lo_FB_input_superclass=Lo_FB_input_superclass +temp;
end
Lo_FB_input_superclass=Lo_FB_input_superclass./i;
Lo_FB_input_superclass=[FAFB_Superclass num2cell(Lo_FB_input_superclass)];

Lop_FB_input_superclass=zeros(size(FAFB_Superclass));
for i=1:1:size(Lop_FB_idx,2)
    temp=cell2mat(RightFB_type.input_superclass{Lop_FB_idx(i)}.syn_count);
    temp=temp/sum(temp);
    Lop_FB_input_superclass=Lop_FB_input_superclass +temp;
end
Lop_FB_input_superclass=Lop_FB_input_superclass./i;
Lop_FB_input_superclass=[FAFB_Superclass num2cell(Lop_FB_input_superclass)];

Multi_FB_input_superclass=zeros(size(FAFB_Superclass));
for i=1:1:size(Multi_FB_idx,2)
    temp=cell2mat(RightFB_type.input_superclass{Multi_FB_idx(i)}.syn_count);
    temp=temp/sum(temp);
    Multi_FB_input_superclass=Multi_FB_input_superclass +temp;
end
Multi_FB_input_superclass=Multi_FB_input_superclass./i;
Multi_FB_input_superclass=[FAFB_Superclass num2cell(Multi_FB_input_superclass)];

Me_FB_output_superclass=zeros(size(FAFB_Superclass));
for i=1:1:size(Me_FB_idx,2)
    temp=cell2mat(RightFB_type.output_superclass{Me_FB_idx(i)}.syn_count);
    temp=temp/sum(temp);
    Me_FB_output_superclass=Me_FB_output_superclass +temp;
end
Me_FB_output_superclass=Me_FB_output_superclass./i;
Me_FB_output_superclass=[FAFB_Superclass num2cell(Me_FB_output_superclass)];

Lo_FB_output_superclass=zeros(size(FAFB_Superclass));
for i=1:1:size(Lo_FB_idx,2)
    temp=cell2mat(RightFB_type.output_superclass{Lo_FB_idx(i)}.syn_count);
    temp=temp/sum(temp);
    Lo_FB_output_superclass=Lo_FB_output_superclass +temp;
end
Lo_FB_output_superclass=Lo_FB_output_superclass./i;
Lo_FB_output_superclass=[FAFB_Superclass num2cell(Lo_FB_output_superclass)];

Lop_FB_output_superclass=zeros(size(FAFB_Superclass));
for i=1:1:size(Lop_FB_idx,2)
    temp=cell2mat(RightFB_type.output_superclass{Lop_FB_idx(i)}.syn_count);
    temp=temp/sum(temp);
    Lop_FB_output_superclass=Lop_FB_output_superclass +temp;
end
Lop_FB_output_superclass=Lop_FB_output_superclass./i;
Lop_FB_output_superclass=[FAFB_Superclass num2cell(Lop_FB_output_superclass)];

Multi_FB_output_superclass=zeros(size(FAFB_Superclass));
for i=1:1:size(Multi_FB_idx,2)
    temp=cell2mat(RightFB_type.output_superclass{Multi_FB_idx(i)}.syn_count);
    temp=temp/sum(temp);
    Multi_FB_output_superclass=Multi_FB_output_superclass +temp;
end
Multi_FB_output_superclass=Multi_FB_output_superclass./i;
Multi_FB_output_superclass=[FAFB_Superclass num2cell(Multi_FB_output_superclass)];
%%
All_idx=[Me_FB_idx Lo_FB_idx Lop_FB_idx Multi_FB_idx];
%%
figure(1);set(gcf,'Color','w'); hold on;
for i=1:1:size(All_idx,2)
    values=cell2mat(RightFB_type.input_superclass{All_idx(i)}.syn_count);
    values=values/sum(values);
    b= bar(i,values,'stacked','EdgeColor','flat');


    for j = 1:length(b)
        b(j).FaceColor = 'flat';      % 'flat'으로 설정
        b(j).CData = colors(j, :); % 각 스택의 색 지정

    end
end
hold off;
set(gca,'TickDir','out','XTick',1:1:size(RightFB_type,1),'XTickLabel',RightFB_type.type(All_idx))
title('Input')
ylim([0 1])
% print(gcf,'-depsc2','-vector','All_FB_input_superclass.eps')
%%
figure(2);set(gcf,'Color','w'); hold on;
for i=1:1:size(All_idx,2)
    values=cell2mat(RightFB_type.output_superclass{All_idx(i)}.syn_count);
    values=values/sum(values);
    b= bar(i,values,'stacked','EdgeColor','flat');


    for j = 1:length(b)
        b(j).FaceColor = 'flat';      % 'flat'으로 설정
        b(j).CData = colors(j, :); % 각 스택의 색 지정

    end
end
hold off;
set(gca,'TickDir','out','XTick',1:1:size(RightFB_type,1),'XTickLabel',RightFB_type.type(All_idx))
ylim([0 1])
title('Output')
% print(gcf,'-depsc2','-vector','All_FB_output_superclass.eps')

%%
figure(3);set(gcf,'Color','w'); hold on;

b= bar(1,cell2mat(Me_FB_input_superclass(:,2)),'stacked','EdgeColor','flat');
for j = 1:length(b)
    b(j).FaceColor = 'flat';      % 'flat'으로 설정
    b(j).CData = colors(j, :); % 각 스택의 색 지정
end

b= bar(2,cell2mat(Lo_FB_input_superclass(:,2)),'stacked','EdgeColor','flat');
for j = 1:length(b)
    b(j).FaceColor = 'flat';      % 'flat'으로 설정
    b(j).CData = colors(j, :); % 각 스택의 색 지정
end

b= bar(3,cell2mat(Lop_FB_input_superclass(:,2)),'stacked','EdgeColor','flat');
for j = 1:length(b)
    b(j).FaceColor = 'flat';      % 'flat'으로 설정
    b(j).CData = colors(j, :); % 각 스택의 색 지정
end

b= bar(4,cell2mat(Multi_FB_input_superclass(:,2)),'stacked','EdgeColor','flat');
for j = 1:length(b)
    b(j).FaceColor = 'flat';      % 'flat'으로 설정
    b(j).CData = colors(j, :); % 각 스택의 색 지정
end
hold off;

set(gca,'TickDir','out','XTick',1:1:4,'XTickLabel',{'Me FB','Lo FB','Lop FB','Multi FB'})
ylim([0 1])
title('Input')
% print(gcf,'-depsc2','-vector','Neuropil_FB_input_superclass.eps')

%%
figure(4);set(gcf,'Color','w'); hold on;

b= bar(1,cell2mat(Me_FB_output_superclass(:,2)),'stacked','EdgeColor','flat');
for j = 1:length(b)
    b(j).FaceColor = 'flat';      % 'flat'으로 설정
    b(j).CData = colors(j, :); % 각 스택의 색 지정
end

b= bar(2,cell2mat(Lo_FB_output_superclass(:,2)),'stacked','EdgeColor','flat');
for j = 1:length(b)
    b(j).FaceColor = 'flat';      % 'flat'으로 설정
    b(j).CData = colors(j, :); % 각 스택의 색 지정
end

b= bar(3,cell2mat(Lop_FB_output_superclass(:,2)),'stacked','EdgeColor','flat');
for j = 1:length(b)
    b(j).FaceColor = 'flat';      % 'flat'으로 설정
    b(j).CData = colors(j, :); % 각 스택의 색 지정
end

b= bar(4,cell2mat(Multi_FB_output_superclass(:,2)),'stacked','EdgeColor','flat');
for j = 1:length(b)
    b(j).FaceColor = 'flat';      % 'flat'으로 설정
    b(j).CData = colors(j, :); % 각 스택의 색 지정
end
hold off;

set(gca,'TickDir','out','XTick',1:1:4,'XTickLabel',{'Me FB','Lo FB','Lop FB','Multi FB'})
ylim([0 1])
title('Output')

%% central 의 top5  input 을 보고 그중에서 다른 종류가 있으면 그걸로 편입
RightFB_type.input_superclass2=RightFB_type.input_superclass;
TopN=5;
for i=1:1:size(RightFB_type,1)
    centralsToFB=RightFB_type.input_superclass{i}.root_ids{2};
    for j=size(centralsToFB,1):-1:1
        inConnections_centralsToFB=FAFBConnections(FAFBConnections.post_root_id==centralsToFB{j,1},:);
        for k=1:1:size(inConnections_centralsToFB,1)
            idx_type=FAFBConsolidated_type.root_id==inConnections_centralsToFB.pre_root_id(k);
            if any(idx_type)
                inConnections_centralsToFB.pre_type{k}=FAFBConsolidated_type.primary_type{idx_type};
            else
                inConnections_centralsToFB.pre_type{k}='NoLabel';
            end
            idx_superclass=FAFBClassification.root_id==inConnections_centralsToFB.pre_root_id(k);
            inConnections_centralsToFB.pre_superclass{k}=FAFBClassification.super_class{idx_superclass};
        end
        if isempty(inConnections_centralsToFB)
            continue;
        end
        noLabel_idx= strcmp(inConnections_centralsToFB.pre_type,'NoLabel');

        inConnections_centralsToFB.pre_superclass(noLabel_idx)={'NoLabel'};

        [unique_inputType_centralsToFB,~,ic]=unique(inConnections_centralsToFB.pre_type);
        for k=1:1:size(unique_inputType_centralsToFB,1)
            idx=ic==k;
            unique_inputType_centralsToFB{k,2}=sum(inConnections_centralsToFB.syn_count(idx));
            [unique_inputType_centralsToFB_superclass,~,icc]=unique(inConnections_centralsToFB.pre_superclass(idx));
            if length(unique_inputType_centralsToFB_superclass)~=1
                temp=inConnections_centralsToFB(idx,:);
                for l=1:1:size(unique_inputType_centralsToFB_superclass,1)
                    idxx=icc==l;
                    unique_inputType_centralsToFB_superclass{l,2}=sum(temp.syn_count(idxx));
                end
                [~,idx_max]=max(cell2mat(unique_inputType_centralsToFB_superclass(:,2)));
                unique_inputType_centralsToFB{k,3}=unique_inputType_centralsToFB_superclass{idx_max,1};
            else
                unique_inputType_centralsToFB(k,3)=unique_inputType_centralsToFB_superclass;
            end
        end
        unique_inputType_centralsToFB=sortrows(unique_inputType_centralsToFB,2,'descend');
        TopN_inputType_centralsToFB=unique_inputType_centralsToFB(1:1:min(TopN,size(unique_inputType_centralsToFB,1)),:);
        superClass_Input_centralsToFB=cell2table(FAFB_Superclass,'VariableNames',{'superclass'});
        for k=1:1:size(superClass_Input_centralsToFB,1)
            idx_superclass=strcmp(TopN_inputType_centralsToFB(:,3),superClass_Input_centralsToFB.superclass{k});
            superClass_Input_centralsToFB.syn_count{k}=sum(cell2mat(TopN_inputType_centralsToFB(idx_superclass,2)));
        end
        %max value로
        % % [~,MaxSuperClass_idx]=max(cell2mat(superClass_Input_centralsToFB.syn_count));
        % if MaxSuperClass_idx~=2
        %     RightFB_type.input_superclass2{i}.syn_count{MaxSuperClass_idx}=...
        %         RightFB_type.input_superclass2{i}.syn_count{MaxSuperClass_idx}+...
        %         RightFB_type.input_superclass2{i}.root_ids{2,1}{j,2};
        %     RightFB_type.input_superclass2{i}.syn_count{2}=...
        %         RightFB_type.input_superclass2{i}.syn_count{2}-...
        %         RightFB_type.input_superclass2{i}.root_ids{2,1}{j,2};
        %     RightFB_type.input_superclass2{i}.root_ids{2,1}(j, :)=[];
        % end
        % central 빼고 비율로
        superClass_Input_centralsToFB.syn_count{2}=0;
        if sum(cell2mat(superClass_Input_centralsToFB.syn_count))==0
            continue;
        else
            RatioSuperClass=cell2mat(superClass_Input_centralsToFB.syn_count)/sum(cell2mat(superClass_Input_centralsToFB.syn_count));
            for k=1:1:size(RatioSuperClass,1)
                if RatioSuperClass(k)==0
                    continue;
                else

                    RightFB_type.input_superclass2{i}.syn_count{k}=...
                        RightFB_type.input_superclass2{i}.syn_count{k}+...
                        RightFB_type.input_superclass2{i}.root_ids{2,1}{j,2}*RatioSuperClass(k);
                    RightFB_type.input_superclass2{i}.syn_count{2}=...
                        RightFB_type.input_superclass2{i}.syn_count{2}-...
                        RightFB_type.input_superclass2{i}.root_ids{2,1}{j,2}*RatioSuperClass(k);
                end
            end

            RightFB_type.input_superclass2{i}.root_ids{2,1}(j, :)=[];

        end
    end

end

%%
figure(5);set(gcf,'Color','w'); hold on;
for i=1:1:size(All_idx,2)
    values=cell2mat(RightFB_type.input_superclass2{All_idx(i)}.syn_count);
    values=values/sum(values);
    b= bar(i,values,'stacked','EdgeColor','flat');


    for j = 1:length(b)
        b(j).FaceColor = 'flat';      % 'flat'으로 설정
        b(j).CData = colors(j, :); % 각 스택의 색 지정

    end
end
hold off;
ylim([0 1])
set(gca,'TickDir','out','XTick',1:1:size(RightFB_type,1),'XTickLabel',RightFB_type.type(All_idx))
title('center neuron replaced input')
%%
Me_FB_input_superclass2=zeros(size(FAFB_Superclass));
for i=1:1:size(Me_FB_idx,2)
    temp=cell2mat(RightFB_type.input_superclass2{Me_FB_idx(i)}.syn_count);
    temp=temp/sum(temp);
    Me_FB_input_superclass2=Me_FB_input_superclass2 +temp;
end
Me_FB_input_superclass2=Me_FB_input_superclass2./i;
Me_FB_input_superclass2=[FAFB_Superclass num2cell(Me_FB_input_superclass2)];

Lo_FB_input_superclass2=zeros(size(FAFB_Superclass));
for i=1:1:size(Lo_FB_idx,2)
    temp=cell2mat(RightFB_type.input_superclass2{Lo_FB_idx(i)}.syn_count);
    temp=temp/sum(temp);
    Lo_FB_input_superclass2=Lo_FB_input_superclass2 +temp;
end
Lo_FB_input_superclass2=Lo_FB_input_superclass2./i;
Lo_FB_input_superclass2=[FAFB_Superclass num2cell(Lo_FB_input_superclass2)];

Lop_FB_input_superclass2=zeros(size(FAFB_Superclass));
for i=1:1:size(Lop_FB_idx,2)
    temp=cell2mat(RightFB_type.input_superclass2{Lop_FB_idx(i)}.syn_count);
    temp=temp/sum(temp);
    Lop_FB_input_superclass2=Lop_FB_input_superclass2 +temp;
end
Lop_FB_input_superclass2=Lop_FB_input_superclass2./i;
Lop_FB_input_superclass2=[FAFB_Superclass num2cell(Lop_FB_input_superclass2)];

Multi_FB_input_superclass2=zeros(size(FAFB_Superclass));
for i=1:1:size(Multi_FB_idx,2)
    temp=cell2mat(RightFB_type.input_superclass2{Multi_FB_idx(i)}.syn_count);
    temp=temp/sum(temp);
    Multi_FB_input_superclass2=Multi_FB_input_superclass2 +temp;
end
Multi_FB_input_superclass2=Multi_FB_input_superclass2./i;
Multi_FB_input_superclass2=[FAFB_Superclass num2cell(Multi_FB_input_superclass2)];

figure(6);set(gcf,'Color','w'); hold on;

b= bar(1,cell2mat(Me_FB_input_superclass2(:,2)),'stacked','EdgeColor','flat');
for j = 1:length(b)
    b(j).FaceColor = 'flat';      % 'flat'으로 설정
    b(j).CData = colors(j, :); % 각 스택의 색 지정
end

b= bar(2,cell2mat(Lo_FB_input_superclass2(:,2)),'stacked','EdgeColor','flat');
for j = 1:length(b)
    b(j).FaceColor = 'flat';      % 'flat'으로 설정
    b(j).CData = colors(j, :); % 각 스택의 색 지정
end

b= bar(3,cell2mat(Lop_FB_input_superclass2(:,2)),'stacked','EdgeColor','flat');
for j = 1:length(b)
    b(j).FaceColor = 'flat';      % 'flat'으로 설정
    b(j).CData = colors(j, :); % 각 스택의 색 지정
end

b= bar(4,cell2mat(Multi_FB_input_superclass2(:,2)),'stacked','EdgeColor','flat');
for j = 1:length(b)
    b(j).FaceColor = 'flat';      % 'flat'으로 설정
    b(j).CData = colors(j, :); % 각 스택의 색 지정
end
hold off;

set(gca,'TickDir','out','XTick',1:1:4,'XTickLabel',{'Me FB','Lo FB','Lop FB','Multi FB'})
ylim([0 1])
title('Input')
title('center neuron replaced input')


%%
data1 = containers.Map(FAFB_Superclass, cell2mat(Me_FB_input_superclass(:,2)));
data2 = containers.Map(FAFB_Superclass, cell2mat(Lo_FB_input_superclass(:,2)));
data3 = containers.Map(FAFB_Superclass, cell2mat(Lop_FB_input_superclass(:,2)));
data4 = containers.Map(FAFB_Superclass, cell2mat(Multi_FB_input_superclass(:,2)));
dataList = {data1, data2, data3, data4};
selected = {'ascending' 'central' 'descending' 'visual_centrifugal' 'visual_projection'};

colorCell = mat2cell(colors, ones(1, size(colors, 1)), 3);
colorMap = containers.Map(FAFB_Superclass, colorCell);

% 선택된 label에 해당하는 색만 추출
SelectedColors = zeros(numel(selected), 3);
for i = 1:numel(selected)
    SelectedColors(i, :) = colorMap(selected{i});
end

% 그래프 그리기
plotStackedSelectedBar(dataList, selected, SelectedColors);
% print(gcf,'-depsc2','-vector','Neuropil_FB_input_superclass.eps')

%%
data1 = containers.Map(FAFB_Superclass, cell2mat(Me_FB_input_superclass2(:,2)));
data2 = containers.Map(FAFB_Superclass, cell2mat(Lo_FB_input_superclass2(:,2)));
data3 = containers.Map(FAFB_Superclass, cell2mat(Lop_FB_input_superclass2(:,2)));
data4 = containers.Map(FAFB_Superclass, cell2mat(Multi_FB_input_superclass2(:,2)));
dataList = {data1, data2, data3, data4};
selected = {'ascending' 'central' 'descending' 'visual_centrifugal' 'visual_projection'};

colorCell = mat2cell(colors, ones(1, size(colors, 1)), 3);
colorMap = containers.Map(FAFB_Superclass, colorCell);

% 선택된 label에 해당하는 색만 추출
SelectedColors = zeros(numel(selected), 3);
for i = 1:numel(selected)
    SelectedColors(i, :) = colorMap(selected{i});
end

% 그래프 그리기
plotStackedSelectedBar(dataList, selected, SelectedColors);
print(gcf,'-depsc2','-vector','Neuropil_FB_input_superclass2.eps')

%%
data1 = containers.Map(FAFB_Superclass, cell2mat(Me_FB_output_superclass(:,2)));
data2 = containers.Map(FAFB_Superclass, cell2mat(Lo_FB_output_superclass(:,2)));
data3 = containers.Map(FAFB_Superclass, cell2mat(Lop_FB_output_superclass(:,2)));
data4 = containers.Map(FAFB_Superclass, cell2mat(Multi_FB_output_superclass(:,2)));
dataList = {data1, data2, data3, data4};
selected = {'optic' 'visual_centrifugal' 'visual_projection'};

colorCell = mat2cell(colors, ones(1, size(colors, 1)), 3);
colorMap = containers.Map(FAFB_Superclass, colorCell);

% 선택된 label에 해당하는 색만 추출
SelectedColors = zeros(numel(selected), 3);
for i = 1:numel(selected)
    SelectedColors(i, :) = colorMap(selected{i});
end

% 그래프 그리기
plotStackedSelectedBar(dataList, selected, SelectedColors);
print(gcf,'-depsc2','-vector','Neuropil_FB_output_superclass.eps')
% print(gcf, 'Lo_Plane.pdf', '-dpdf', '-r600')


% view(-58.316292922587124,10.953666510684501)
% axis equal
%% 평균 + 표준편차 계산
n_superclass = length(FAFB_Superclass);
n_neurons = height(RightFB_type);
all_ratios = zeros(n_neurons, n_superclass); % 각 neuron의 normalized 비율 저장

for i = 1:n_neurons
    temp = cell2mat(RightFB_type.input_superclass2{i}.syn_count);
    temp = temp / sum(temp); % 정규화
    all_ratios(i, :) = temp;
end

% 평균과 표준편차
mean_input_superclass = mean(all_ratios, 1);
std_input_superclass = std(all_ratios, 0, 1); % 0: N-1로 나누는 표준편차

% 테이블로 저장 (선택사항)
Total_FB_input_superclass = [FAFB_Superclass num2cell(mean_input_superclass') num2cell(std_input_superclass')];
Total_FB_input_superclass.Properties.VariableNames = {'superclass', 'mean_ratio', 'std_ratio'};
