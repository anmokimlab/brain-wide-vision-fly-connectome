%% 기존 코드 동일
clear all; close all; clc;

opt = detectImportOptions('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\connections_no_threshold.csv');
opt = setvartype(opt,'pre_root_id','int64');
opt = setvartype(opt,'post_root_id','int64');
FAFBConnections = readtable('C:\Users\user\Documents\MATLAB\FAFB_v2\V738_V5_20250730\CodexData\connections_no_threshold.csv',opt);

FAFBConnections = sortrows(FAFBConnections, {'pre_root_id','post_root_id','syn_count'}, {'ascend','ascend','descend'});

% Pre = FAFBConnections.pre_root_id;
% Pre_shifted = Pre([2:end end]);
% Post = FAFBConnections.post_root_id;
% Post_shifted = Post([2:end end]);
% Pre_Post_Same = find((Pre == Pre_shifted) & (Post == Post_shifted));
% Pre_Post_Same(end) = [];
% 
% for i = size(Pre_Post_Same, 1):-1:1
%     FAFBConnections.syn_count(Pre_Post_Same(i)) = ...
%         FAFBConnections.syn_count(Pre_Post_Same(i)) + ...
%         FAFBConnections.syn_count(Pre_Post_Same(i) + 1);
%     disp(i / size(Pre_Post_Same, 1) * 100)
% end
% 
% FAFBConnections(Pre_Post_Same+1, :) = [];


% 2. 이전 row와 같은 (pre, post) 쌍인지 검사
same_as_previous = [false; ...
    FAFBConnections.pre_root_id(2:end) == FAFBConnections.pre_root_id(1:end-1) & ...
    FAFBConnections.post_root_id(2:end) == FAFBConnections.post_root_id(1:end-1)];

% 3. 같은 쌍마다 고유한 group ID 부여
group_id = cumsum(~same_as_previous);

% 4. group마다 syn_count 합산
[G, pre_group] = findgroups(group_id);
syn_sum = splitapply(@sum, FAFBConnections.syn_count, G);

% 5. group 내 첫 번째 row만 남기고 나머지 제거
first_in_group = [true; diff(group_id) ~= 0];
FAFBConnections = FAFBConnections(first_in_group, :);

% 6. 합쳐진 syn_count 대입
FAFBConnections.syn_count = syn_sum;
%% Make Graph
rootIds = unique([FAFBConnections.pre_root_id; FAFBConnections.post_root_id]);
[~, preIdx] = ismember(FAFBConnections.pre_root_id, rootIds);
[~, postIdx] = ismember(FAFBConnections.post_root_id, rootIds);
weight = FAFBConnections.syn_count;

G = digraph(preIdx, postIdx, weight);

%% 🧠 Centrality 계산

% (1) PageRank: 가중치 반영
pagerank_vals = centrality(G, 'pagerank');

% (2) Betweenness (비가중치)
betweenness_vals_unweighted = centrality(G, 'betweenness');

% (3) Betweenness (가중치 고려: weight가 클수록 짧은 거리로 간주)
edge_costs = 1 ./ (G.Edges.Weight + eps);  % 역가중치, 0 나눔 방지
betweenness_vals_weighted = centrality(G, 'betweenness', 'Cost', edge_costs);

%% 결과 저장
save("Allgraph_Thr0.mat", "G", "postIdx", "preIdx", "rootIds", "weight", ...
    "pagerank_vals", "betweenness_vals_unweighted", "betweenness_vals_weighted");

%% 🔍 특정 root_id의 centrality 값 추출
% % 예시: 특정 root_id의 PageRank와 Betweenness 값을 보고 싶을 때
% target_root_id = 720575940612764252;  % 원하는 root_id로 바꿔
% 
% % 인덱스 찾기
% idx = find(rootIds == target_root_id);
% 
% if ~isempty(idx)
%     fprintf('Root ID: %d\n', target_root_id);
%     fprintf('  PageRank: %.6f\n', pagerank_vals(idx));
%     fprintf('  Betweenness: %.6f\n', betweenness_vals(idx));
% else
%     warning('해당 root_id를 찾을 수 없습니다.');
% end