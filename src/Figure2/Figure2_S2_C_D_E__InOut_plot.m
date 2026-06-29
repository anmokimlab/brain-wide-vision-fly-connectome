clear all; close all ;clc
load('FF_Optic_Central_FanInOut_new.mat')
% % load('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\Allgraph_Thr0.mat')
% % % % load('RightFF_InOut.mat')
% % % % load('Optic_InOut.mat')
% % % % load('Central_InOut.mat')
% % % %%
% [UniqueRightFF_type,~,ic]=unique(RightFF_InOut.type);
% UniqueRightFF=table(UniqueRightFF_type,'VariableNames',{'type'});
% 
% 
% isEmptyIn = cellfun(@(c) isempty(c) || all(cellfun(@isempty, c(:))), RightFF_InOut.In);
% isEmptyOut = cellfun(@(c) isempty(c) || all(cellfun(@isempty, c(:))), RightFF_InOut.Out);
% 
% for i=1:1:size(UniqueRightFF,1)
%     idx=ic==i;
%     EIn_N=sum(idx&~isEmptyIn);
%     EOut_N=sum(idx&~isEmptyOut);
% 
%     UniqueRightFF.NofType(i)=sum(idx);
%     if ~EIn_N==0
% 
%         tempInStats=RightFF_InOut.In(idx);
%         tempInStats=vertcat(tempInStats{:});
%         tempInStats=sortrows(tempInStats,1,'ascend');
%         for j=size(tempInStats,1):-1:2
%             if strcmp(tempInStats{j,1},tempInStats{j-1,1})
%                 tempInStats{j-1,2}=[tempInStats{j,2};tempInStats{j-1,2}];
%                 tempInStats{j-1,3}=tempInStats{j-1,3}+tempInStats{j,3};
%                 tempInStats{j-1,4}=tempInStats{j-1,4}+tempInStats{j,4};
%                 tempInStats(j,:)=[];
%             end
%         end
%         tempInStats(:,3)=num2cell(cell2mat(tempInStats(:,3))/EIn_N);
%         UniqueRightFF.In{i}=tempInStats;
%     end
%     if ~EOut_N==0
% 
%         tempOutStats=RightFF_InOut.Out(idx);
%         tempOutStats=vertcat(tempOutStats{:});
%         tempOutStats=sortrows(tempOutStats,1,'ascend');
%         for j=size(tempOutStats,1):-1:2
%             if strcmp(tempOutStats{j,1},tempOutStats{j-1,1})
%                 tempOutStats{j-1,2}=[tempOutStats{j,2};tempOutStats{j-1,2}];
%                 tempOutStats{j-1,3}=tempOutStats{j-1,3}+tempOutStats{j,3};
%                 tempOutStats{j-1,4}=tempOutStats{j-1,4}+tempOutStats{j,4};
%                 tempOutStats(j,:)=[];
%             end
%         end
%         tempOutStats(:,3)=num2cell(cell2mat(tempOutStats(:,3))/EOut_N);
%         UniqueRightFF.Out{i}=tempOutStats;
%     end
% end
% 
% 
% %%
% [UniqueCentral_type,~,ic]=unique(Central_InOut.type);
% UniqueCentral=table(UniqueCentral_type,'VariableNames',{'type'});
% 
% 
% isEmptyIn = cellfun(@(c) isempty(c) || all(cellfun(@isempty, c(:))), Central_InOut.In);
% isEmptyOut = cellfun(@(c) isempty(c) || all(cellfun(@isempty, c(:))), Central_InOut.Out);
% 
% for i=1:1:size(UniqueCentral,1)
%     idx=ic==i;
%     EIn_N=sum(idx&~isEmptyIn);
%     EOut_N=sum(idx&~isEmptyOut);
%     UniqueCentral.NofType(i)=sum(idx);
%     if ~EIn_N==0
%         tempInStats=Central_InOut.In(idx);
%         tempInStats=vertcat(tempInStats{:});
%         tempInStats=sortrows(tempInStats,1,'ascend');
%         for j=size(tempInStats,1):-1:2
%             if strcmp(tempInStats{j,1},tempInStats{j-1,1})
%                 tempInStats{j-1,2}=[tempInStats{j,2};tempInStats{j-1,2}];
%                 tempInStats{j-1,3}=tempInStats{j-1,3}+tempInStats{j,3};
%                 tempInStats{j-1,4}=tempInStats{j-1,4}+tempInStats{j,4};
%                 tempInStats(j,:)=[];
%             end
%         end
%         tempInStats(:,3)=num2cell(cell2mat(tempInStats(:,3))/EIn_N);
%         UniqueCentral.In{i}=tempInStats;
%     end
%     if ~EOut_N==0
%         tempOutStats=Central_InOut.Out(idx);
%         tempOutStats=vertcat(tempOutStats{:});
%         tempOutStats=sortrows(tempOutStats,1,'ascend');
%         for j=size(tempOutStats,1):-1:2
%             if strcmp(tempOutStats{j,1},tempOutStats{j-1,1})
%                 tempOutStats{j-1,2}=[tempOutStats{j,2};tempOutStats{j-1,2}];
%                 tempOutStats{j-1,3}=tempOutStats{j-1,3}+tempOutStats{j,3};
%                 tempOutStats{j-1,4}=tempOutStats{j-1,4}+tempOutStats{j,4};
%                 tempOutStats(j,:)=[];
%             end
%         end
%         tempOutStats(:,3)=num2cell(cell2mat(tempOutStats(:,3))/EOut_N);
%         UniqueCentral.Out{i}=tempOutStats;
%     end
% end
% 
% %%
% [UniqueOptic_type,~,ic]=unique(Optic_InOut.type);
% UniqueOptic=table(UniqueOptic_type,'VariableNames',{'type'});
% 
% 
% isEmptyIn = cellfun(@(c) isempty(c) || all(cellfun(@isempty, c(:))), Optic_InOut.In);
% isEmptyOut = cellfun(@(c) isempty(c) || all(cellfun(@isempty, c(:))), Optic_InOut.Out);
% 
% for i=1:1:size(UniqueOptic,1)
%     idx=ic==i;
%     EIn_N=sum(idx&~isEmptyIn);
%     EOut_N=sum(idx&~isEmptyOut);
% 
%     UniqueOptic.NofType(i)=sum(idx);
%     if ~EIn_N==0
%         tempInStats=Optic_InOut.In(idx);
%         tempInStats=vertcat(tempInStats{:});
%         tempInStats=sortrows(tempInStats,1,'ascend');
%         for j=size(tempInStats,1):-1:2
%             if strcmp(tempInStats{j,1},tempInStats{j-1,1})
%                 tempInStats{j-1,2}=[tempInStats{j,2};tempInStats{j-1,2}];
%                 tempInStats{j-1,3}=tempInStats{j-1,3}+tempInStats{j,3};
%                 tempInStats{j-1,4}=tempInStats{j-1,4}+tempInStats{j,4};
%                 tempInStats(j,:)=[];
%             end
%         end
%         tempInStats(:,3)=num2cell(cell2mat(tempInStats(:,3))/EIn_N);
%         UniqueOptic.In{i}=tempInStats;
%     end
%     if ~EOut_N==0
%         tempOutStats=Optic_InOut.Out(idx);
%         tempOutStats=vertcat(tempOutStats{:});
%         tempOutStats=sortrows(tempOutStats,1,'ascend');
%         for j=size(tempOutStats,1):-1:2
%             if strcmp(tempOutStats{j,1},tempOutStats{j-1,1})
%                 tempOutStats{j-1,2}=[tempOutStats{j,2};tempOutStats{j-1,2}];
%                 tempOutStats{j-1,3}=tempOutStats{j-1,3}+tempOutStats{j,3};
%                 tempOutStats{j-1,4}=tempOutStats{j-1,4}+tempOutStats{j,4};
%                 tempOutStats(j,:)=[];
%             end
%         end
%         tempOutStats(:,3)=num2cell(cell2mat(tempOutStats(:,3))/EOut_N);
%         UniqueOptic.Out{i}=tempOutStats;
%     end
% end
% % %% load data
% % % load('FF_Optic_Central_FanInOut.mat')
% % idx = find(rootIds == target_root_id);
% % 
% % if ~isempty(idx)
% %     fprintf('Root ID: %d\n', target_root_id);
% %     fprintf('  PageRank: %.6f\n', pagerank_vals(idx));
% %     fprintf('  Betweenness: %.6f\n', betweenness_vals(idx));
% % else
% %     warning('해당 root_id를 찾을 수 없습니다.');
% % end
% %%
% for i=1:1:size(RightFF_InOut,1)
%         inStats=RightFF_InOut.In{i};
%         outStats=RightFF_InOut.Out{i};
%         if ~isempty(inStats)
%             RightFF_InOut.InNeuronNumber{i}=sum(cell2mat(inStats(:,3)));
%             RightFF_InOut.InNeuronTypeNumber{i}=size(inStats,1);
%             RightFF_InOut.InNeuronRatio{i}=RightFF_InOut.InNeuronNumber{i}/RightFF_InOut.InNeuronTypeNumber{i};
% 
%         end
%         if ~isempty(outStats)
%             RightFF_InOut.OutNeuronNumber{i}=sum(cell2mat(outStats(:,3)));
%             RightFF_InOut.OutNeuronTypeNumber{i}=size(outStats,1);
%             RightFF_InOut.OutNeuronRatio{i}=RightFF_InOut.OutNeuronNumber{i}/RightFF_InOut.OutNeuronTypeNumber{i};
% 
%         end
%     target_root_id=RightFF_InOut.root_id(i);
%     idx = find(rootIds == target_root_id);
%     RightFF_InOut.betweenness_unweighted(i)=betweenness_vals_unweighted(idx);
%     RightFF_InOut.betweenness_weighted(i)=betweenness_vals_weighted(idx);
%     RightFF_InOut.pagerank(i)=pagerank_vals(idx);
% 
% 
% end
% 
% for i=1:1:size(Central_InOut,1)
%         inStats=Central_InOut.In{i};
%         outStats=Central_InOut.Out{i};
%         if ~isempty(inStats)
%             Central_InOut.InNeuronNumber{i}=sum(cell2mat(inStats(:,3)));
%             Central_InOut.InNeuronTypeNumber{i}=size(inStats,1);
%             Central_InOut.InNeuronRatio{i}=Central_InOut.InNeuronNumber{i}/Central_InOut.InNeuronTypeNumber{i};
%         end
%         if ~isempty(outStats)
%             Central_InOut.OutNeuronNumber{i}=sum(cell2mat(outStats(:,3)));
%             Central_InOut.OutNeuronTypeNumber{i}=size(outStats,1);
%             Central_InOut.OutNeuronRatio{i}=Central_InOut.OutNeuronNumber{i}/Central_InOut.OutNeuronTypeNumber{i};
%         end
% 
%     target_root_id=Central_InOut.root_id(i);
%     idx = find(rootIds == target_root_id);
%     if ~isempty(idx)
%         Central_InOut.betweenness_unweighted(i)=betweenness_vals_unweighted(idx);
%         Central_InOut.betweenness_weighted(i)=betweenness_vals_weighted(idx);
%         Central_InOut.pagerank(i)=pagerank_vals(idx);
%     end
% end
% 
% for i=1:1:size(Optic_InOut,1)
%         inStats=Optic_InOut.In{i};
%         outStats=Optic_InOut.Out{i};
%         if ~isempty(inStats)
%             Optic_InOut.InNeuronNumber{i}=sum(cell2mat(inStats(:,3)));
%             Optic_InOut.InNeuronTypeNumber{i}=size(inStats,1);
%             Optic_InOut.InNeuronRatio{i}=Optic_InOut.InNeuronNumber{i}/Optic_InOut.InNeuronTypeNumber{i};
% 
%         end
%         if ~isempty(outStats)
%             Optic_InOut.OutNeuronNumber{i}=sum(cell2mat(outStats(:,3)));
%             Optic_InOut.OutNeuronTypeNumber{i}=size(outStats,1);
%             Optic_InOut.OutNeuronRatio{i}=Optic_InOut.OutNeuronNumber{i}/Optic_InOut.OutNeuronTypeNumber{i};
% 
%         end
%     target_root_id=Optic_InOut.root_id(i);
%     idx = find(rootIds == target_root_id);
%     if ~isempty(idx)
%         Optic_InOut.betweenness_unweighted(i)=betweenness_vals_unweighted(idx);
%         Optic_InOut.betweenness_weighted(i)=betweenness_vals_weighted(idx);
%         Optic_InOut.pagerank(i)=pagerank_vals(idx);
%     end
% end
% %
% [type_RightFF_InOut,~,ic]=unique(RightFF_InOut.type);
% type_RightFF_InOut=table(type_RightFF_InOut,'VariableNames',{'type'});
% for i=1:1:size(type_RightFF_InOut,1)
%     idx=ic==i;
%     type_RightFF_InOut.InNeuronNumber{i}=mean(cell2mat(RightFF_InOut.InNeuronNumber(idx)),'omitmissing');
%     type_RightFF_InOut.InNeuronTypeNumber{i}=mean(cell2mat(RightFF_InOut.InNeuronTypeNumber(idx)),'omitmissing');
%     type_RightFF_InOut.InNeuronRatio{i}=mean(cell2mat(RightFF_InOut.InNeuronRatio(idx)),'omitmissing');
% 
%     type_RightFF_InOut.OutNeuronNumber{i}=mean(cell2mat(RightFF_InOut.OutNeuronNumber(idx)),'omitmissing');
%     type_RightFF_InOut.OutNeuronTypeNumber{i}=mean(cell2mat(RightFF_InOut.OutNeuronTypeNumber(idx)),'omitmissing');
%     type_RightFF_InOut.OutNeuronRatio{i}=mean(cell2mat(RightFF_InOut.OutNeuronRatio(idx)),'omitmissing');
% 
%     type_RightFF_InOut.betweenness_unweighted{i}=mean((RightFF_InOut.betweenness_unweighted(idx)),'omitmissing');
%     type_RightFF_InOut.betweenness_weighted{i}=mean((RightFF_InOut.betweenness_weighted(idx)),'omitmissing');
%     type_RightFF_InOut.pagerank{i}=mean((RightFF_InOut.pagerank(idx)),'omitmissing');
% 
% end
% 
% [type_Central_InOut,~,ic]=unique(Central_InOut.type);
% type_Central_InOut=table(type_Central_InOut,'VariableNames',{'type'});
% for i=1:1:size(type_Central_InOut,1)
%     idx=ic==i;
%     type_Central_InOut.InNeuronNumber{i}=mean(cell2mat(Central_InOut.InNeuronNumber(idx)),'omitmissing');
%     type_Central_InOut.InNeuronTypeNumber{i}=mean(cell2mat(Central_InOut.InNeuronTypeNumber(idx)),'omitmissing');
%     type_Central_InOut.InNeuronRatio{i}=mean(cell2mat(Central_InOut.InNeuronRatio(idx)),'omitmissing');
% 
%     type_Central_InOut.OutNeuronNumber{i}=mean(cell2mat(Central_InOut.OutNeuronNumber(idx)),'omitmissing');
%     type_Central_InOut.OutNeuronTypeNumber{i}=mean(cell2mat(Central_InOut.OutNeuronTypeNumber(idx)),'omitmissing');
%     type_Central_InOut.OutNeuronRatio{i}=mean(cell2mat(Central_InOut.OutNeuronRatio(idx)),'omitmissing');
% 
%     type_Central_InOut.betweenness_unweighted{i}=mean((Central_InOut.betweenness_unweighted(idx)),'omitmissing');
%     type_Central_InOut.betweenness_weighted{i}=mean((Central_InOut.betweenness_weighted(idx)),'omitmissing');
%     type_Central_InOut.pagerank{i}=mean((Central_InOut.pagerank(idx)),'omitmissing');
% 
% end
% 
% [type_Optic_InOut,~,ic]=unique(Optic_InOut.type);
% type_Optic_InOut=table(type_Optic_InOut,'VariableNames',{'type'});
% for i=1:1:size(type_Optic_InOut,1)
%     idx=ic==i;
%     type_Optic_InOut.InNeuronNumber{i}=mean(cell2mat(Optic_InOut.InNeuronNumber(idx)),'omitmissing');
%     type_Optic_InOut.InNeuronTypeNumber{i}=mean(cell2mat(Optic_InOut.InNeuronTypeNumber(idx)),'omitmissing');
%     type_Optic_InOut.InNeuronRatio{i}=mean(cell2mat(Optic_InOut.InNeuronRatio(idx)),'omitmissing');
% 
%     type_Optic_InOut.OutNeuronNumber{i}=mean(cell2mat(Optic_InOut.OutNeuronNumber(idx)),'omitmissing');
%     type_Optic_InOut.OutNeuronTypeNumber{i}=mean(cell2mat(Optic_InOut.OutNeuronTypeNumber(idx)),'omitmissing');
%     type_Optic_InOut.OutNeuronRatio{i}=mean(cell2mat(Optic_InOut.OutNeuronRatio(idx)),'omitmissing');
% 
%     type_Optic_InOut.betweenness_unweighted{i}=mean((Optic_InOut.betweenness_unweighted(idx)),'omitmissing');
%     type_Optic_InOut.betweenness_weighted{i}=mean((Optic_InOut.betweenness_weighted(idx)),'omitmissing');
%     type_Optic_InOut.pagerank{i}=mean((Optic_InOut.pagerank(idx)),'omitmissing');
% 
% end

%% 1. in neuron number...

% 샘플 데이터 (길이가 다른 데이터)
data1 = cell2mat(type_Optic_InOut.InNeuronNumber);
data2 = cell2mat(type_RightFF_InOut.InNeuronNumber);
data3 = cell2mat(type_Central_InOut.InNeuronNumber);

% 데이터와 그룹 정의
data = {data1, data2, data3}; % 셀 배열로 저장
group = repelem(1:length(data), cellfun(@length, data)); % 그룹 레이블 생성
combinedData = vertcat(data{:}); % 데이터를 하나의 배열로 결합
figure(1);set(gcf,'Color','w')
% 박스플롯 생성
boxplot(combinedData, group, 'Labels', {'Optic', 'Feedforward', 'Central'},'Notch','on','Symbol','');
% boxplot(combinedData, group, 'Labels', {'Optic', 'Feedforward', 'Central'},'Notch','on');

ylabel('InNeuronNumber');
% xlim([0, 4]);
ylim([0, 2000]);
% 박스 색 변경
set(gca,'Box','off','TickDir','out')

ax = gca;
boxes = findobj(ax, 'Tag', 'Box');         % 박스 핸들
boxes(3).Color = [0.9290 0.6940 0.1250]; % 박스 색 (예: 청록색)
boxes(2).Color = [0 0.4470 0.7410]; % 박스 색 (예: 청록색)
boxes(1).Color = [0.8500 0.3250 0.0980]; % 박스 색 (예: 청록색)
medians = findobj(ax, 'Tag', 'Median');   % 중앙값 핸들
medians(3).Color = [0.9290 0.6940 0.1250]; % 박스 색 (예: 청록색)
medians(2).Color = [0 0.4470 0.7410]; % 박스 색 (예: 청록색)
medians(1).Color = [0.8500 0.3250 0.0980]; % 박스 색 (예: 청록색)

%%%% 비모수 검정
% 데이터 준비

% 그룹 이름 설정
groupNames = {'Optic', 'Feedforward', 'Central'};

% 비교쌍 정보
dataPairs = {
    data1, data2;
    data1, data3;
    data2, data3
    };

groupPairs = {
    [1, 2];
    [1, 3];
    [2, 3]
    };

% 보정된 유의수준 (Bonferroni)
alpha = 0.05;
numComparisons = size(dataPairs, 1);
correctedAlpha = alpha / numComparisons;

% 결과 출력
fprintf('=== Number of In Neuron Pairwise Comparison (Bonferroni corrected α = %.4f) ===\n\n', correctedAlpha);

for i = 1:numComparisons
    d1 = dataPairs{i, 1};
    d2 = dataPairs{i, 2};
    g1 = groupNames{groupPairs{i}(1)};
    g2 = groupNames{groupPairs{i}(2)};

    % Wilcoxon rank-sum test (비모수, 두 그룹 간 중앙값 비교)
    p = ranksum(d1, d2);

    % Probabilistic dominance (벡터화)
    p_dominance = mean(d1 > d2', 'all');
    p_dominance2 = mean(d1 < d2', 'all');

    % 결과 출력
    fprintf('%s vs %s:\n', g1, g2);
    fprintf('  p-value = %.4f --> %s\n', p, ...
        ternary(p < correctedAlpha, 'Significant ✅', 'Not significant ❌'));
    fprintf('  P(%s > %s) = %.3f\n', g1, g2, p_dominance);
    fprintf('  P(%s < %s) = %.3f\n\n', g1, g2, p_dominance2);

end
print(gcf,'-depsc2','-vector','figure2_InNeuronNumber.eps')

%% 2. in neuron types...

% 샘플 데이터 (길이가 다른 데이터)
data1 = cell2mat(type_Optic_InOut.InNeuronTypeNumber);
data2 = cell2mat(type_RightFF_InOut.InNeuronTypeNumber);
data3 = cell2mat(type_Central_InOut.InNeuronTypeNumber);

% 데이터와 그룹 정의
data = {data1, data2, data3}; % 셀 배열로 저장
group = repelem(1:length(data), cellfun(@length, data)); % 그룹 레이블 생성
combinedData = vertcat(data{:}); % 데이터를 하나의 배열로 결합
figure(2);set(gcf,'Color','w')
% 박스플롯 생성
boxplot(combinedData, group, 'Labels', {'Optic', 'Feedforward', 'Central'},'Notch','on','Symbol','');
% boxplot(combinedData, group, 'Labels', {'Optic', 'Feedforward', 'Central'},'Notch','on');

ylabel('InNeuronTypeNumber');
% xlim([0, 4]);
ylim([0, 500]);
% 박스 색 변경
set(gca,'Box','off','TickDir','out')

ax = gca;
boxes = findobj(ax, 'Tag', 'Box');         % 박스 핸들
boxes(3).Color = [0.9290 0.6940 0.1250]; % 박스 색 (예: 청록색)
boxes(2).Color = [0 0.4470 0.7410]; % 박스 색 (예: 청록색)
boxes(1).Color = [0.8500 0.3250 0.0980]; % 박스 색 (예: 청록색)
medians = findobj(ax, 'Tag', 'Median');   % 중앙값 핸들
medians(3).Color = [0.9290 0.6940 0.1250]; % 박스 색 (예: 청록색)
medians(2).Color = [0 0.4470 0.7410]; % 박스 색 (예: 청록색)
medians(1).Color = [0.8500 0.3250 0.0980]; % 박스 색 (예: 청록색)

%%%% 비모수 검정
% 데이터 준비

% 그룹 이름 설정
groupNames = {'Optic', 'Feedforward', 'Central'};

% 비교쌍 정보
dataPairs = {
    data1, data2;
    data1, data3;
    data2, data3
    };

groupPairs = {
    [1, 2];
    [1, 3];
    [2, 3]
    };

% 보정된 유의수준 (Bonferroni)
alpha = 0.05;
numComparisons = size(dataPairs, 1);
correctedAlpha = alpha / numComparisons;

% 결과 출력
fprintf('=== Number of In Neuron Type Pairwise Comparison (Bonferroni corrected α = %.4f) ===\n\n', correctedAlpha);

for i = 1:numComparisons
    d1 = dataPairs{i, 1};
    d2 = dataPairs{i, 2};
    g1 = groupNames{groupPairs{i}(1)};
    g2 = groupNames{groupPairs{i}(2)};

    % Wilcoxon rank-sum test (비모수, 두 그룹 간 중앙값 비교)
    p = ranksum(d1, d2);

    % Probabilistic dominance (벡터화)
    p_dominance = mean(d1 > d2', 'all');
    p_dominance2 = mean(d1 < d2', 'all');

    % 결과 출력
    fprintf('%s vs %s:\n', g1, g2);
    fprintf('  p-value = %.4f --> %s\n', p, ...
        ternary(p < correctedAlpha, 'Significant ✅', 'Not significant ❌'));
    fprintf('  P(%s > %s) = %.3f\n', g1, g2, p_dominance);
    fprintf('  P(%s < %s) = %.3f\n\n', g1, g2, p_dominance2);

end

print(gcf,'-depsc2','-vector','figure2_InNeuronType.eps')


%% 3. in neuron Ratio...

% 샘플 데이터 (길이가 다른 데이터)
data1 = cell2mat(type_Optic_InOut.InNeuronRatio);
data2 = cell2mat(type_RightFF_InOut.InNeuronRatio);
data3 = cell2mat(type_Central_InOut.InNeuronRatio);

% 데이터와 그룹 정의
data = {data1, data2, data3}; % 셀 배열로 저장
group = repelem(1:length(data), cellfun(@length, data)); % 그룹 레이블 생성
combinedData = vertcat(data{:}); % 데이터를 하나의 배열로 결합
figure(3);set(gcf,'Color','w')
% 박스플롯 생성
boxplot(combinedData, group, 'Labels', {'Optic', 'Feedforward', 'Central'},'Notch','on','Symbol','');
% boxplot(combinedData, group, 'Labels', {'Optic', 'Feedforward', 'Central'},'Notch','on');

ylabel('InNeuronRatio');
% xlim([0, 4]);
ylim([0,15]);
% 박스 색 변경
set(gca,'Box','off','TickDir','out')

ax = gca;
boxes = findobj(ax, 'Tag', 'Box');         % 박스 핸들
boxes(3).Color = [0.9290 0.6940 0.1250]; % 박스 색 (예: 청록색)
boxes(2).Color = [0 0.4470 0.7410]; % 박스 색 (예: 청록색)
boxes(1).Color = [0.8500 0.3250 0.0980]; % 박스 색 (예: 청록색)
medians = findobj(ax, 'Tag', 'Median');   % 중앙값 핸들
medians(3).Color = [0.9290 0.6940 0.1250]; % 박스 색 (예: 청록색)
medians(2).Color = [0 0.4470 0.7410]; % 박스 색 (예: 청록색)
medians(1).Color = [0.8500 0.3250 0.0980]; % 박스 색 (예: 청록색)

%%%% 비모수 검정
% 데이터 준비

% 그룹 이름 설정
groupNames = {'Optic', 'Feedforward', 'Central'};

% 비교쌍 정보
dataPairs = {
    data1, data2;
    data1, data3;
    data2, data3
    };

groupPairs = {
    [1, 2];
    [1, 3];
    [2, 3]
    };

% 보정된 유의수준 (Bonferroni)
alpha = 0.05;
numComparisons = size(dataPairs, 1);
correctedAlpha = alpha / numComparisons;

% 결과 출력
fprintf('=== InNeuron Ratio Pairwise Comparison (Bonferroni corrected α = %.4f) ===\n\n', correctedAlpha);

for i = 1:numComparisons
    d1 = dataPairs{i, 1};
    d2 = dataPairs{i, 2};
    g1 = groupNames{groupPairs{i}(1)};
    g2 = groupNames{groupPairs{i}(2)};

    % Wilcoxon rank-sum test (비모수, 두 그룹 간 중앙값 비교)
    p = ranksum(d1, d2);

    % Probabilistic dominance (벡터화)
    p_dominance = mean(d1 > d2', 'all');
    p_dominance2 = mean(d1 < d2', 'all');

    % 결과 출력
    fprintf('%s vs %s:\n', g1, g2);
    fprintf('  p-value = %.4f --> %s\n', p, ...
        ternary(p < correctedAlpha, 'Significant ✅', 'Not significant ❌'));
    fprintf('  P(%s > %s) = %.3f\n', g1, g2, p_dominance);
    fprintf('  P(%s < %s) = %.3f\n\n', g1, g2, p_dominance2);

end
print(gcf,'-depsc2','-vector','figure2_InNeuronRatio.eps')

%% OUT
%% 4. out neuron number...

% 샘플 데이터 (길이가 다른 데이터)
data1 = cell2mat(type_Optic_InOut.OutNeuronNumber);
data2 = cell2mat(type_RightFF_InOut.OutNeuronNumber);
data3 = cell2mat(type_Central_InOut.OutNeuronNumber);

% 데이터와 그룹 정의
data = {data1, data2, data3}; % 셀 배열로 저장
group = repelem(1:length(data), cellfun(@length, data)); % 그룹 레이블 생성
combinedData = vertcat(data{:}); % 데이터를 하나의 배열로 결합
figure(4);set(gcf,'Color','w')
% 박스플롯 생성
boxplot(combinedData, group, 'Labels', {'Optic', 'Feedforward', 'Central'},'Notch','on','Symbol','');
% boxplot(combinedData, group, 'Labels', {'Optic', 'Feedforward', 'Central'},'Notch','on');

ylabel('OutNeuronNumber');
% xlim([0, 4]);
ylim([0, 1200]);
% 박스 색 변경
set(gca,'Box','off','TickDir','out')

ax = gca;
boxes = findobj(ax, 'Tag', 'Box');         % 박스 핸들
boxes(3).Color = [0.9290 0.6940 0.1250]; % 박스 색 (예: 청록색)
boxes(2).Color = [0 0.4470 0.7410]; % 박스 색 (예: 청록색)
boxes(1).Color = [0.8500 0.3250 0.0980]; % 박스 색 (예: 청록색)
medians = findobj(ax, 'Tag', 'Median');   % 중앙값 핸들
medians(3).Color = [0.9290 0.6940 0.1250]; % 박스 색 (예: 청록색)
medians(2).Color = [0 0.4470 0.7410]; % 박스 색 (예: 청록색)
medians(1).Color = [0.8500 0.3250 0.0980]; % 박스 색 (예: 청록색)

%%%% 비모수 검정
% 데이터 준비

% 그룹 이름 설정
groupNames = {'Optic', 'Feedforward', 'Central'};

% 비교쌍 정보
dataPairs = {
    data1, data2;
    data1, data3;
    data2, data3
    };

groupPairs = {
    [1, 2];
    [1, 3];
    [2, 3]
    };

% 보정된 유의수준 (Bonferroni)
alpha = 0.05;
numComparisons = size(dataPairs, 1);
correctedAlpha = alpha / numComparisons;

% 결과 출력
fprintf('=== OutNeuronNumber Pairwise Comparison (Bonferroni corrected α = %.4f) ===\n\n', correctedAlpha);

for i = 1:numComparisons
    d1 = dataPairs{i, 1};
    d2 = dataPairs{i, 2};
    g1 = groupNames{groupPairs{i}(1)};
    g2 = groupNames{groupPairs{i}(2)};

    % Wilcoxon rank-sum test (비모수, 두 그룹 간 중앙값 비교)
    p = ranksum(d1, d2);

    % Probabilistic dominance (벡터화)
    p_dominance = mean(d1 > d2', 'all');
    p_dominance2 = mean(d1 < d2', 'all');

    % 결과 출력
    fprintf('%s vs %s:\n', g1, g2);
    fprintf('  p-value = %.4f --> %s\n', p, ...
        ternary(p < correctedAlpha, 'Significant ✅', 'Not significant ❌'));
    fprintf('  P(%s > %s) = %.3f\n', g1, g2, p_dominance);
    fprintf('  P(%s < %s) = %.3f\n\n', g1, g2, p_dominance2);

end
print(gcf,'-depsc2','-vector','figure2_OutNeuronNumber.eps')


%% 5. out neuron type...

% 샘플 데이터 (길이가 다른 데이터)
data1 = cell2mat(type_Optic_InOut.OutNeuronTypeNumber);
data2 = cell2mat(type_RightFF_InOut.OutNeuronTypeNumber);
data3 = cell2mat(type_Central_InOut.OutNeuronTypeNumber);

% 데이터와 그룹 정의
data = {data1, data2, data3}; % 셀 배열로 저장
group = repelem(1:length(data), cellfun(@length, data)); % 그룹 레이블 생성
combinedData = vertcat(data{:}); % 데이터를 하나의 배열로 결합
figure(5);set(gcf,'Color','w')
% 박스플롯 생성
boxplot(combinedData, group, 'Labels', {'Optic', 'Feedforward', 'Central'},'Notch','on','Symbol','');
% boxplot(combinedData, group, 'Labels', {'Optic', 'Feedforward', 'Central'},'Notch','on');

ylabel('OutNeuronType');
% xlim([0, 4]);
ylim([0, 600]);
% 박스 색 변경
set(gca,'Box','off','TickDir','out')

ax = gca;
boxes = findobj(ax, 'Tag', 'Box');         % 박스 핸들
boxes(3).Color = [0.9290 0.6940 0.1250]; % 박스 색 (예: 청록색)
boxes(2).Color = [0 0.4470 0.7410]; % 박스 색 (예: 청록색)
boxes(1).Color = [0.8500 0.3250 0.0980]; % 박스 색 (예: 청록색)
medians = findobj(ax, 'Tag', 'Median');   % 중앙값 핸들
medians(3).Color = [0.9290 0.6940 0.1250]; % 박스 색 (예: 청록색)
medians(2).Color = [0 0.4470 0.7410]; % 박스 색 (예: 청록색)
medians(1).Color = [0.8500 0.3250 0.0980]; % 박스 색 (예: 청록색)

%%%% 비모수 검정
% 데이터 준비

% 그룹 이름 설정
groupNames = {'Optic', 'Feedforward', 'Central'};

% 비교쌍 정보
dataPairs = {
    data1, data2;
    data1, data3;
    data2, data3
    };

groupPairs = {
    [1, 2];
    [1, 3];
    [2, 3]
    };

% 보정된 유의수준 (Bonferroni)
alpha = 0.05;
numComparisons = size(dataPairs, 1);
correctedAlpha = alpha / numComparisons;

% 결과 출력
fprintf('=== OutNeuronType Pairwise Comparison (Bonferroni corrected α = %.4f) ===\n\n', correctedAlpha);

for i = 1:numComparisons
    d1 = dataPairs{i, 1};
    d2 = dataPairs{i, 2};
    g1 = groupNames{groupPairs{i}(1)};
    g2 = groupNames{groupPairs{i}(2)};

    % Wilcoxon rank-sum test (비모수, 두 그룹 간 중앙값 비교)
    p = ranksum(d1, d2);

    % Probabilistic dominance (벡터화)
    p_dominance = mean(d1 > d2', 'all');
    p_dominance2 = mean(d1 < d2', 'all');

    % 결과 출력
    fprintf('%s vs %s:\n', g1, g2);
    fprintf('  p-value = %.4f --> %s\n', p, ...
        ternary(p < correctedAlpha, 'Significant ✅', 'Not significant ❌'));
    fprintf('  P(%s > %s) = %.3f\n', g1, g2, p_dominance);
    fprintf('  P(%s < %s) = %.3f\n\n', g1, g2, p_dominance2);

end
print(gcf,'-depsc2','-vector','figure2_OutNeuronType.eps')


%% 5. out neuron Ratio...

% 샘플 데이터 (길이가 다른 데이터)
data1 = cell2mat(type_Optic_InOut.OutNeuronRatio);
data2 = cell2mat(type_RightFF_InOut.OutNeuronRatio);
data3 = cell2mat(type_Central_InOut.OutNeuronRatio);

% 데이터와 그룹 정의
data = {data1, data2, data3}; % 셀 배열로 저장
group = repelem(1:length(data), cellfun(@length, data)); % 그룹 레이블 생성
combinedData = vertcat(data{:}); % 데이터를 하나의 배열로 결합
figure(6);set(gcf,'Color','w')
% 박스플롯 생성
boxplot(combinedData, group, 'Labels', {'Optic', 'Feedforward', 'Central'},'Notch','on','Symbol','');
% boxplot(combinedData, group, 'Labels', {'Optic', 'Feedforward', 'Central'},'Notch','on');

ylabel('OutNeuronRatio');
% xlim([0, 4]);
ylim([0, 15]);
% 박스 색 변경
set(gca,'Box','off','TickDir','out')

ax = gca;
boxes = findobj(ax, 'Tag', 'Box');         % 박스 핸들
boxes(3).Color = [0.9290 0.6940 0.1250]; % 박스 색 (예: 청록색)
boxes(2).Color = [0 0.4470 0.7410]; % 박스 색 (예: 청록색)
boxes(1).Color = [0.8500 0.3250 0.0980]; % 박스 색 (예: 청록색)
medians = findobj(ax, 'Tag', 'Median');   % 중앙값 핸들
medians(3).Color = [0.9290 0.6940 0.1250]; % 박스 색 (예: 청록색)
medians(2).Color = [0 0.4470 0.7410]; % 박스 색 (예: 청록색)
medians(1).Color = [0.8500 0.3250 0.0980]; % 박스 색 (예: 청록색)

%%%% 비모수 검정
% 데이터 준비

% 그룹 이름 설정
groupNames = {'Optic', 'Feedforward', 'Central'};

% 비교쌍 정보
dataPairs = {
    data1, data2;
    data1, data3;
    data2, data3
    };

groupPairs = {
    [1, 2];
    [1, 3];
    [2, 3]
    };

% 보정된 유의수준 (Bonferroni)
alpha = 0.05;
numComparisons = size(dataPairs, 1);
correctedAlpha = alpha / numComparisons;

% 결과 출력
fprintf('=== OutNeuronRatio Pairwise Comparison (Bonferroni corrected α = %.4f) ===\n\n', correctedAlpha);

for i = 1:numComparisons
    d1 = dataPairs{i, 1};
    d2 = dataPairs{i, 2};
    g1 = groupNames{groupPairs{i}(1)};
    g2 = groupNames{groupPairs{i}(2)};

    % Wilcoxon rank-sum test (비모수, 두 그룹 간 중앙값 비교)
    p = ranksum(d1, d2);

    % Probabilistic dominance (벡터화)
    p_dominance = mean(d1 > d2', 'all');
    p_dominance2 = mean(d1 < d2', 'all');

    % 결과 출력
    fprintf('%s vs %s:\n', g1, g2);
    fprintf('  p-value = %.4f --> %s\n', p, ...
        ternary(p < correctedAlpha, 'Significant ✅', 'Not significant ❌'));
    fprintf('  P(%s > %s) = %.3f\n', g1, g2, p_dominance);
    fprintf('  P(%s < %s) = %.3f\n\n', g1, g2, p_dominance2);

end
print(gcf,'-depsc2','-vector','figure2_OutNeuronRatio.eps')
%

%% 5. betweeness unweighted...

% 샘플 데이터 (길이가 다른 데이터)
data1 = cell2mat(type_Optic_InOut.betweenness_unweighted);
data2 = cell2mat(type_RightFF_InOut.betweenness_unweighted);
data3 = cell2mat(type_Central_InOut.betweenness_unweighted);

% 데이터와 그룹 정의
data = {data1, data2, data3}; % 셀 배열로 저장
group = repelem(1:length(data), cellfun(@length, data)); % 그룹 레이블 생성
combinedData = vertcat(data{:}); % 데이터를 하나의 배열로 결합
figure(7);set(gcf,'Color','w')
% 박스플롯 생성
boxplot(combinedData, group, 'Labels', {'Optic', 'Feedforward', 'Central'},'Notch','on','Symbol','');
% boxplot(combinedData, group, 'Labels', {'Optic', 'Feedforward', 'Central'},'Notch','on');

ylabel('betweenness_unweighted');
% xlim([0, 4]);
ylim([0, 1.2e7]);
% 박스 색 변경
set(gca,'Box','off','TickDir','out')

ax = gca;
boxes = findobj(ax, 'Tag', 'Box');         % 박스 핸들
boxes(3).Color = [0.9290 0.6940 0.1250]; % 박스 색 (예: 청록색)
boxes(2).Color = [0 0.4470 0.7410]; % 박스 색 (예: 청록색)
boxes(1).Color = [0.8500 0.3250 0.0980]; % 박스 색 (예: 청록색)
medians = findobj(ax, 'Tag', 'Median');   % 중앙값 핸들
medians(3).Color = [0.9290 0.6940 0.1250]; % 박스 색 (예: 청록색)
medians(2).Color = [0 0.4470 0.7410]; % 박스 색 (예: 청록색)
medians(1).Color = [0.8500 0.3250 0.0980]; % 박스 색 (예: 청록색)

%%%% 비모수 검정
% 데이터 준비

% 그룹 이름 설정
groupNames = {'Optic', 'Feedforward', 'Central'};

% 비교쌍 정보
dataPairs = {
    data1, data2;
    data1, data3;
    data2, data3
    };

groupPairs = {
    [1, 2];
    [1, 3];
    [2, 3]
    };

% 보정된 유의수준 (Bonferroni)
alpha = 0.05;
numComparisons = size(dataPairs, 1);
correctedAlpha = alpha / numComparisons;

% 결과 출력
fprintf('=== Betweeness unweighted Comparison (Bonferroni corrected α = %.4f) ===\n\n', correctedAlpha);

for i = 1:numComparisons
    d1 = dataPairs{i, 1};
    d2 = dataPairs{i, 2};
    g1 = groupNames{groupPairs{i}(1)};
    g2 = groupNames{groupPairs{i}(2)};

    % Wilcoxon rank-sum test (비모수, 두 그룹 간 중앙값 비교)
    p = ranksum(d1, d2);

    % Probabilistic dominance (벡터화)
    p_dominance = mean(d1 > d2', 'all');
    p_dominance2 = mean(d1 < d2', 'all');

    % 결과 출력
    fprintf('%s vs %s:\n', g1, g2);
    fprintf('  p-value = %.4f --> %s\n', p, ...
        ternary(p < correctedAlpha, 'Significant ✅', 'Not significant ❌'));
    fprintf('  P(%s > %s) = %.3f\n', g1, g2, p_dominance);
    fprintf('  P(%s < %s) = %.3f\n\n', g1, g2, p_dominance2);

end
print(gcf,'-depsc2','-vector','betweeness unweighted.eps')
%
%% 5. betweeness weighted...

% 샘플 데이터 (길이가 다른 데이터)
data1 = cell2mat(type_Optic_InOut.betweenness_weighted);
data2 = cell2mat(type_RightFF_InOut.betweenness_weighted);
data3 = cell2mat(type_Central_InOut.betweenness_weighted);

% 데이터와 그룹 정의
data = {data1, data2, data3}; % 셀 배열로 저장
group = repelem(1:length(data), cellfun(@length, data)); % 그룹 레이블 생성
combinedData = vertcat(data{:}); % 데이터를 하나의 배열로 결합
figure(8);set(gcf,'Color','w')
% 박스플롯 생성
boxplot(combinedData, group, 'Labels', {'Optic', 'Feedforward', 'Central'},'Notch','on','Symbol','');
% boxplot(combinedData, group, 'Labels', {'Optic', 'Feedforward', 'Central'},'Notch','on');

ylabel('betweenness_weighted');
% xlim([0, 4]);
ylim([0, 1.2e7]);
% 박스 색 변경
set(gca,'Box','off','TickDir','out')

ax = gca;
boxes = findobj(ax, 'Tag', 'Box');         % 박스 핸들
boxes(3).Color = [0.9290 0.6940 0.1250]; % 박스 색 (예: 청록색)
boxes(2).Color = [0 0.4470 0.7410]; % 박스 색 (예: 청록색)
boxes(1).Color = [0.8500 0.3250 0.0980]; % 박스 색 (예: 청록색)
medians = findobj(ax, 'Tag', 'Median');   % 중앙값 핸들
medians(3).Color = [0.9290 0.6940 0.1250]; % 박스 색 (예: 청록색)
medians(2).Color = [0 0.4470 0.7410]; % 박스 색 (예: 청록색)
medians(1).Color = [0.8500 0.3250 0.0980]; % 박스 색 (예: 청록색)

%%%% 비모수 검정
% 데이터 준비

% 그룹 이름 설정
groupNames = {'Optic', 'Feedforward', 'Central'};

% 비교쌍 정보
dataPairs = {
    data1, data2;
    data1, data3;
    data2, data3
    };

groupPairs = {
    [1, 2];
    [1, 3];
    [2, 3]
    };

% 보정된 유의수준 (Bonferroni)
alpha = 0.05;
numComparisons = size(dataPairs, 1);
correctedAlpha = alpha / numComparisons;

% 결과 출력
fprintf('=== Betweeness weighted Comparison (Bonferroni corrected α = %.4f) ===\n\n', correctedAlpha);

for i = 1:numComparisons
    d1 = dataPairs{i, 1};
    d2 = dataPairs{i, 2};
    g1 = groupNames{groupPairs{i}(1)};
    g2 = groupNames{groupPairs{i}(2)};

    % Wilcoxon rank-sum test (비모수, 두 그룹 간 중앙값 비교)
    p = ranksum(d1, d2);

    % Probabilistic dominance (벡터화)
    p_dominance = mean(d1 > d2', 'all');
    p_dominance2 = mean(d1 < d2', 'all');

    % 결과 출력
    fprintf('%s vs %s:\n', g1, g2);
    fprintf('  p-value = %.4f --> %s\n', p, ...
        ternary(p < correctedAlpha, 'Significant ✅', 'Not significant ❌'));
    fprintf('  P(%s > %s) = %.3f\n', g1, g2, p_dominance);
    fprintf('  P(%s < %s) = %.3f\n\n', g1, g2, p_dominance2);

end
print(gcf,'-depsc2','-vector','betweeness weighted.eps')

%% 5. page rank...

% 샘플 데이터 (길이가 다른 데이터)
data1 = cell2mat(type_Optic_InOut.pagerank);
data2 = cell2mat(type_RightFF_InOut.pagerank);
data3 = cell2mat(type_Central_InOut.pagerank);

% 데이터와 그룹 정의
data = {data1, data2, data3}; % 셀 배열로 저장
group = repelem(1:length(data), cellfun(@length, data)); % 그룹 레이블 생성
combinedData = vertcat(data{:}); % 데이터를 하나의 배열로 결합
figure(9);set(gcf,'Color','w')
% 박스플롯 생성
boxplot(combinedData, group, 'Labels', {'Optic', 'Feedforward', 'Central'},'Notch','on','Symbol','');
% boxplot(combinedData, group, 'Labels', {'Optic', 'Feedforward', 'Central'},'Notch','on');

ylabel('pagerank');
% xlim([0, 4]);
ylim([0, 10e-5]);
% 박스 색 변경
set(gca,'Box','off','TickDir','out')

ax = gca;
boxes = findobj(ax, 'Tag', 'Box');         % 박스 핸들
boxes(3).Color = [0.9290 0.6940 0.1250]; % 박스 색 (예: 청록색)
boxes(2).Color = [0 0.4470 0.7410]; % 박스 색 (예: 청록색)
boxes(1).Color = [0.8500 0.3250 0.0980]; % 박스 색 (예: 청록색)
medians = findobj(ax, 'Tag', 'Median');   % 중앙값 핸들
medians(3).Color = [0.9290 0.6940 0.1250]; % 박스 색 (예: 청록색)
medians(2).Color = [0 0.4470 0.7410]; % 박스 색 (예: 청록색)
medians(1).Color = [0.8500 0.3250 0.0980]; % 박스 색 (예: 청록색)

%%%% 비모수 검정
% 데이터 준비

% 그룹 이름 설정
groupNames = {'Optic', 'Feedforward', 'Central'};

% 비교쌍 정보
dataPairs = {
    data1, data2;
    data1, data3;
    data2, data3
    };

groupPairs = {
    [1, 2];
    [1, 3];
    [2, 3]
    };

% 보정된 유의수준 (Bonferroni)
alpha = 0.05;
numComparisons = size(dataPairs, 1);
correctedAlpha = alpha / numComparisons;

% 결과 출력
fprintf('=== Pagerank Comparison (Bonferroni corrected α = %.4f) ===\n\n', correctedAlpha);

for i = 1:numComparisons
    d1 = dataPairs{i, 1};
    d2 = dataPairs{i, 2};
    g1 = groupNames{groupPairs{i}(1)};
    g2 = groupNames{groupPairs{i}(2)};

    % Wilcoxon rank-sum test (비모수, 두 그룹 간 중앙값 비교)
    p = ranksum(d1, d2);

    % Probabilistic dominance (벡터화)
    p_dominance = mean(d1 > d2', 'all');
    p_dominance2 = mean(d1 < d2', 'all');

    % 결과 출력
    fprintf('%s vs %s:\n', g1, g2);
    fprintf('  p-value = %.4f --> %s\n', p, ...
        ternary(p < correctedAlpha, 'Significant ✅', 'Not significant ❌'));
    fprintf('  P(%s > %s) = %.3f\n', g1, g2, p_dominance);
    fprintf('  P(%s < %s) = %.3f\n\n', g1, g2, p_dominance2);

end
print(gcf,'-depsc2','-vector','pagerank.eps')

%% Figure s2. ff 하나씩 그리기
N_bar=10;
figure(10);set(gcf,'Color','w','InnerPosition',[551,202,621,291],'OuterPosition',[543,194,637,385])
temp=cell2mat(type_RightFF_InOut.InNeuronNumber);
[sorted_InNeuronNumber,sorted_idx]=sort(temp,'descend');

bar(type_RightFF_InOut.type(sorted_idx(1:1:N_bar)),sorted_InNeuronNumber(1:1:N_bar),'FaceColor','#0072BD','EdgeColor','flat');
set(gca,'Box','off','TickDir','out')
title('In Neuron Number FF')
print(gcf,'-depsc2','-vector','FigureS2_InNeuronNumber FF.eps')

figure(11);set(gcf,'Color','w','InnerPosition',[551,202,621,291],'OuterPosition',[543,194,637,385])
temp=cell2mat(type_Optic_InOut.InNeuronNumber);
[sorted_InNeuronNumber,sorted_idx]=sort(temp,'descend','MissingPlacement','last');

bar(type_Optic_InOut.type(sorted_idx(1:1:N_bar)),sorted_InNeuronNumber(1:1:N_bar),'FaceColor','#EDB120','EdgeColor','flat');
set(gca,'Box','off','TickDir','out')
title('In Neuron Number Optic')
print(gcf,'-depsc2','-vector','FigureS2_InNeuronNumber Optic.eps')

figure(12);set(gcf,'Color','w','InnerPosition',[551,202,621,291],'OuterPosition',[543,194,637,385])
temp=cell2mat(type_Central_InOut.InNeuronNumber);
[sorted_InNeuronNumber,sorted_idx]=sort(temp,'descend','MissingPlacement','last');

bar(type_Central_InOut.type(sorted_idx(1:1:N_bar)),sorted_InNeuronNumber(1:1:N_bar),'FaceColor','#7E2F8E','EdgeColor','flat');
set(gca,'Box','off','TickDir','out')
title('In Neuron Number Central')
print(gcf,'-depsc2','-vector','FigureS2_InNeuronNumber Central.eps')

%%
figure(13);set(gcf,'Color','w','InnerPosition',[551,202,621,291],'OuterPosition',[543,194,637,385])
temp=cell2mat(type_RightFF_InOut.InNeuronTypeNumber);
[sorted_InNeuronTypeNumber,sorted_idx]=sort(temp,'descend');

bar(type_RightFF_InOut.type(sorted_idx(1:1:N_bar)),sorted_InNeuronTypeNumber(1:1:N_bar),'FaceColor','#0072BD','EdgeColor','flat');
set(gca,'Box','off','TickDir','out')
title('In Neuron Type Number FF')
print(gcf,'-depsc2','-vector','FigureS2_InNeuronTypeNumber FF.eps')
% 
figure(14);set(gcf,'Color','w','InnerPosition',[551,202,621,291],'OuterPosition',[543,194,637,385])
temp=cell2mat(type_Optic_InOut.InNeuronTypeNumber);
[sorted_InNeuronTypeNumber,sorted_idx]=sort(temp,'descend','MissingPlacement','last');

bar(type_Optic_InOut.type(sorted_idx(1:1:N_bar)),sorted_InNeuronTypeNumber(1:1:N_bar),'FaceColor','#EDB120','EdgeColor','flat');
set(gca,'Box','off','TickDir','out')
title('In Neuron Type Number Optic')
print(gcf,'-depsc2','-vector','FigureS2_InNeuronTypeNumber Optic.eps')

figure(15);set(gcf,'Color','w','InnerPosition',[551,202,621,291],'OuterPosition',[543,194,637,385])
temp=cell2mat(type_Central_InOut.InNeuronTypeNumber);
[sorted_InNeuronTypeNumber,sorted_idx]=sort(temp,'descend','MissingPlacement','last');

bar(type_Central_InOut.type(sorted_idx(1:1:N_bar)),sorted_InNeuronTypeNumber(1:1:N_bar),'FaceColor','#7E2F8E','EdgeColor','flat');
set(gca,'Box','off','TickDir','out')
title('In Neuron Type Number Central')
print(gcf,'-depsc2','-vector','FigureS2_InNeuronTypeNumber Central.eps')
% 
%%
figure(16);set(gcf,'Color','w','InnerPosition',[551,202,621,291],'OuterPosition',[543,194,637,385])
temp=cell2mat(type_RightFF_InOut.OutNeuronNumber);
[sorted_OutNeuronNumber,sorted_idx]=sort(temp,'descend');

bar(type_RightFF_InOut.type(sorted_idx(1:1:N_bar)),sorted_OutNeuronNumber(1:1:N_bar),'FaceColor','#0072BD','EdgeColor','flat');
set(gca,'Box','off','TickDir','out')
title('Out Neuron Number FF')
print(gcf,'-depsc2','-vector','FigureS2_OutNeuronNumber FF.eps')

figure(17);set(gcf,'Color','w','InnerPosition',[551,202,621,291],'OuterPosition',[543,194,637,385])
temp=cell2mat(type_Optic_InOut.OutNeuronNumber);
[sorted_OutNeuronNumber,sorted_idx]=sort(temp,'descend','MissingPlacement','last');

bar(type_Optic_InOut.type(sorted_idx(1:1:N_bar)),sorted_OutNeuronNumber(1:1:N_bar),'FaceColor','#EDB120','EdgeColor','flat');
set(gca,'Box','off','TickDir','out')
title('Out Neuron Number Optic')
print(gcf,'-depsc2','-vector','FigureS2_OutNeuronNumber Optic.eps')

figure(18);set(gcf,'Color','w','InnerPosition',[551,202,621,291],'OuterPosition',[543,194,637,385])
temp=cell2mat(type_Central_InOut.OutNeuronNumber);
[sorted_OutNeuronNumber,sorted_idx]=sort(temp,'descend','MissingPlacement','last');

bar(type_Central_InOut.type(sorted_idx(1:1:N_bar)),sorted_OutNeuronNumber(1:1:N_bar),'FaceColor','#7E2F8E','EdgeColor','flat');
set(gca,'Box','off','TickDir','out')
title('Out Neuron Number Central')
print(gcf,'-depsc2','-vector','FigureS2_OutNeuronNumber Central.eps')
%%
figure(19);set(gcf,'Color','w','InnerPosition',[551,202,621,291],'OuterPosition',[543,194,637,385])
temp=cell2mat(type_RightFF_InOut.OutNeuronTypeNumber);
[sorted_OutNeuronTypeNumber,sorted_idx]=sort(temp,'descend');

bar(type_RightFF_InOut.type(sorted_idx(1:1:N_bar)),sorted_OutNeuronTypeNumber(1:1:N_bar),'FaceColor','#0072BD','EdgeColor','flat');
set(gca,'Box','off','TickDir','out')
title('Out Neuron Type Number FF')
print(gcf,'-depsc2','-vector','FigureS2_OutNeuronTypeNumber FF.eps')

figure(20);set(gcf,'Color','w','InnerPosition',[551,202,621,291],'OuterPosition',[543,194,637,385])
temp=cell2mat(type_Optic_InOut.OutNeuronTypeNumber);
[sorted_OutNeuronTypeNumber,sorted_idx]=sort(temp,'descend','MissingPlacement','last');

bar(type_Optic_InOut.type(sorted_idx(1:1:N_bar)),sorted_OutNeuronTypeNumber(1:1:N_bar),'FaceColor','#EDB120','EdgeColor','flat');
set(gca,'Box','off','TickDir','out')
title('Out Neuron Type Number Optic')
print(gcf,'-depsc2','-vector','FigureS2_OutNeuronTypeNumber Optic.eps')

figure(21);set(gcf,'Color','w','InnerPosition',[551,202,621,291],'OuterPosition',[543,194,637,385])
temp=cell2mat(type_Central_InOut.OutNeuronTypeNumber);
[sorted_OutNeuronTypeNumber,sorted_idx]=sort(temp,'descend','MissingPlacement','last');

bar(type_Central_InOut.type(sorted_idx(1:1:N_bar)),sorted_OutNeuronTypeNumber(1:1:N_bar),'FaceColor','#7E2F8E','EdgeColor','flat');
set(gca,'Box','off','TickDir','out')
title('Out Neuron Type Number Central')
print(gcf,'-depsc2','-vector','FigureS2_OutNeuronTypeNumber Central.eps')
%%
figure(22);set(gcf,'Color','w','InnerPosition',[551,202,621,291],'OuterPosition',[543,194,637,385])
temp=cell2mat(type_RightFF_InOut.betweenness_unweighted);
[sorted_betweenness_unweighted,sorted_idx]=sort(temp,'descend');

bar(type_RightFF_InOut.type(sorted_idx(1:1:N_bar)),sorted_betweenness_unweighted(1:1:N_bar),'FaceColor','#0072BD','EdgeColor','flat');
set(gca,'Box','off','TickDir','out')
title('betweenness_unweighted FF')
print(gcf,'-depsc2','-vector','FigureS2_betweenness_unweighted FF.eps')

figure(23);set(gcf,'Color','w','InnerPosition',[551,202,621,291],'OuterPosition',[543,194,637,385])
temp=cell2mat(type_Optic_InOut.betweenness_unweighted);
[sorted_betweenness_unweighted,sorted_idx]=sort(temp,'descend','MissingPlacement','last');

bar(type_Optic_InOut.type(sorted_idx(1:1:N_bar)),sorted_betweenness_unweighted(1:1:N_bar),'FaceColor','#EDB120','EdgeColor','flat');
set(gca,'Box','off','TickDir','out')
title('betweenness_unweighted Optic')
print(gcf,'-depsc2','-vector','FigureS2_betweenness_unweighted Optic.eps')

figure(24);set(gcf,'Color','w','InnerPosition',[551,202,621,291],'OuterPosition',[543,194,637,385])
temp=cell2mat(type_Central_InOut.betweenness_unweighted);
[sorted_betweenness_unweighted,sorted_idx]=sort(temp,'descend','MissingPlacement','last');

bar(type_Central_InOut.type(sorted_idx(1:1:N_bar)),sorted_betweenness_unweighted(1:1:N_bar),'FaceColor','#7E2F8E','EdgeColor','flat');
set(gca,'Box','off','TickDir','out')
title('betweenness_unweighted Central')
print(gcf,'-depsc2','-vector','FigureS2_betweenness_unweighted Central.eps')


%%
figure(25);set(gcf,'Color','w','InnerPosition',[551,202,621,291],'OuterPosition',[543,194,637,385])
temp=cell2mat(type_RightFF_InOut.pagerank);
[sorted_pagerank,sorted_idx]=sort(temp,'descend');

bar(type_RightFF_InOut.type(sorted_idx(1:1:N_bar)),sorted_pagerank(1:1:N_bar),'FaceColor','#0072BD','EdgeColor','flat');
set(gca,'Box','off','TickDir','out')
title('pagerank FF')
print(gcf,'-depsc2','-vector','FigureS2_pagerank FF.eps')

figure(26);set(gcf,'Color','w','InnerPosition',[551,202,621,291],'OuterPosition',[543,194,637,385])
temp=cell2mat(type_Optic_InOut.pagerank);
[sorted_pagerank,sorted_idx]=sort(temp,'descend','MissingPlacement','last');

bar(type_Optic_InOut.type(sorted_idx(1:1:N_bar)),sorted_pagerank(1:1:N_bar),'FaceColor','#EDB120','EdgeColor','flat');
set(gca,'Box','off','TickDir','out')
title('pagerank Optic')
print(gcf,'-depsc2','-vector','FigureS2_pagerank Optic.eps')

figure(27);set(gcf,'Color','w','InnerPosition',[551,202,621,291],'OuterPosition',[543,194,637,385])
temp=cell2mat(type_Central_InOut.pagerank);
[sorted_pagerank,sorted_idx]=sort(temp,'descend','MissingPlacement','last');

bar(type_Central_InOut.type(sorted_idx(1:1:N_bar)),sorted_pagerank(1:1:N_bar),'FaceColor','#7E2F8E','EdgeColor','flat');
set(gca,'Box','off','TickDir','out')
title('pagerank Central')
print(gcf,'-depsc2','-vector','FigureS2_pagerank Central.eps')

%%
%% ===== Figure S2 — FF only, 8 metrics × 10 bars (fixed order) =====
% 요구 테이블: type_RightFF_InOut (FF만)
% 필수 변수:
%   - type (cellstr or string)
%   - InNeuronNumber, InNeuronTypeNumber, InNeuronRatio
%   - OutNueronNumber(오탈자) 또는 OutNeuronNumber
%   - OutNueronTypeNumber(오탈자) 또는 OutNeuronTypeNumber
%   - OutNeuronRatio
%   - betweenness_unweighted, pagerank
% 선택/안전: NaN/빈값은 뒤로 보냄

N_bar = 10;   % 각 지표에서 뽑을 Top-K
T = type_RightFF_InOut;

% -- 도우미: 변수 존재/선택 & cell->double 캐스팅
toNum = @(x) double(cell2mat(x));
hasVar = @(nm) any(strcmpi(T.Properties.VariableNames, nm));
pickVar = @(cands) cands{find(cellfun(hasVar,cands),1,'first')}; % 후보 중 첫 존재 변수명

% getNumeric = @(nm) ...
%     ( iscell(T.(nm)) ? toNum(T.(nm)) : double(T.(nm)) ); %#ok<NASGU>

% MATLAB에는 삼항연산자(?)가 없으므로 함수로 대체
function v = getNumVec(T, nm)
    if iscell(T.(nm)), v = double(cell2mat(T.(nm)));
    else,               v = double(T.(nm));
    end
end

% -- 각 지표 벡터 (오탈자 호환)
InNum   = getNumVec(T, pickVar({"InNeuronNumber"}));
InTyp   = getNumVec(T, pickVar({"InNeuronTypeNumber"}));
InRatio = getNumVec(T, pickVar({"InNeuronRatio"}));

OutNumVar = pickVar({"OutNeuronNumber"});
OutNum    = getNumVec(T, OutNumVar);

OutTypVar = pickVar({"OutNeuronTypeNumber"});
OutTyp    = getNumVec(T, OutTypVar);

OutRatio  = getNumVec(T, pickVar({"OutNeuronRatio"}));

Between   = getNumVec(T, pickVar({"betweenness_unweighted"}));
PageRank  = getNumVec(T, pickVar({"pagerank"}));

% -- 상위 10개 index씩 추출 (NaN/결측은 뒤로)
getTopIdx = @(v) sortrows([(1:numel(v))', v(:)], -2, 'MissingPlacement','last');

blk1 = getTopIdx(InNum );     blk1 = blk1(1:min(N_bar,end),1);
blk2 = getTopIdx(InTyp );     blk2 = blk2(1:min(N_bar,end),1);
blk3 = getTopIdx(InRatio);    blk3 = blk3(1:min(N_bar,end),1);
blk4 = getTopIdx(OutNum);     blk4 = blk4(1:min(N_bar,end),1);
blk5 = getTopIdx(OutTyp);     blk5 = blk5(1:min(N_bar,end),1);
blk6 = getTopIdx(OutRatio);   blk6 = blk6(1:min(N_bar,end),1);
blk7 = getTopIdx(Between);    blk7 = blk7(1:min(N_bar,end),1);
blk8 = getTopIdx(PageRank);   blk8 = blk8(1:min(N_bar,end),1);

% -- 고정 x-순서(총 80개; 중복 허용)
orderIdx = [blk1; blk2; blk3; blk4; blk5; blk6; blk7; blk8];

% -- 라벨(타입명만; 필요시 블록/순위 접두어 복원 가능)
types = string(T.type(orderIdx));
labels = types;

% -- x 구분선 위치 (8블록 → 7개 구분선)
nBlocks = 8;
cutXs = (N_bar+0.5) + (0:(nBlocks-2))*N_bar;  % 10.5,20.5,...,70.5

% -- 그릴 지표들(각 지표마다 y만 바꿔서 동일 x-순서로 그림)
METRICS = {InNum, InTyp, InRatio, OutNum, OutTyp, OutRatio, Between, PageRank};
TITLES  = { ...
    'FF — Presynaptic neuron count (Top10×8 order)', ...
    'FF — Presynaptic type count (Top10×8 order)', ...
    'FF — Presynaptic neurons per type (Top10×8 order)', ...
    'FF — Postsynaptic neuron count (Top10×8 order)', ...
    'FF — Postsynaptic type count (Top10×8 order)', ...
    'FF — Postsynaptic neurons per type (Top10×8 order)', ...
    'FF — Betweenness (unweighted) (Top10×8 order)', ...
    'FF — PageRank (Top10×8 order)'};
FILES   = { ...
    'FigureS2_FF_only_PresynNeuronCount.eps', ...
    'FigureS2_FF_only_PresynTypeCount.eps', ...
    'FigureS2_FF_only_PresynPerType.eps', ...
    'FigureS2_FF_only_PostsynNeuronCount.eps', ...
    'FigureS2_FF_only_PostsynTypeCount.eps', ...
    'FigureS2_FF_only_PostsynPerType.eps', ...
    'FigureS2_FF_only_Betweenness.eps', ...
    'FigureS2_FF_only_PageRank.eps'};

% -- 블록명(상단 텍스트용)
blockNames = ["Pre n","Pre type","Pre perType","Post n","Post type","Post perType","Betweenness","PageRank"];

% -- 공통 Figure 스타일
figPosInner = [551,202,1300,380];  % 80 bars라 가로 더 넉넉히
barColor    = '#0072BD';           % FF 컬러 고정

for m = 1:numel(METRICS)
    y = METRICS{m}(orderIdx);

    figure(100+m); clf; set(gcf,'Color','w','InnerPosition',figPosInner);
    bar(y, 'FaceColor', barColor, 'EdgeColor','flat');
    xlim([0.5, numel(orderIdx)+0.5])
    set(gca,'Box','off','TickDir','out')
    title(TITLES{m}, 'Interpreter','none')
    xticks(1:numel(orderIdx));
    xticklabels(labels);
    xtickangle(90);
    ylabel('Value')

    % 블록 구분선
    hold on;
    for cx = cutXs
        xline(cx, '-', 'Color', [0.5 0.5 0.5], 'Alpha', 0.3);
    end
    % 블록명 상단 표시
    yTop = max(y(~isnan(y))) * 1.05;
    for b = 1:nBlocks
        xc = ( (b-1)*N_bar + 0.5 + (b*N_bar) ) / 2;  % 블록 중앙
        text(xc, yTop, blockNames(b), 'HorizontalAlignment','center');
    end

    % 저장 (벡터)
    print(gcf,'-depsc2','-painters','-vector', FILES{m});
end

%% 간단한 ternary function 정의
function out = ternary(cond, valTrue, valFalse)
if cond
    out = valTrue;
else
    out = valFalse;
end
end

