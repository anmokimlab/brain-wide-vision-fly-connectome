%% =========================================================
%  SuperCluster x Neuropil 값을 직접 입력 → 임계값 이하 제거 → 유향 그래프
%  작성자: ChatGPT
% =========================================================
close all; clear; clc;

%% 1) 노드 이름 정의
SC  = {'SuperCluster1','SuperCluster2','SuperCluster3','SuperCluster4', ...
       'SuperCluster5','SuperCluster6','SuperCluster7'};
NP  = {'PVLP','AVLP','PLP','IPS','SPS','AOTU','SCL','SLP','ICL','VES'};

%% 2) 값 직접 입력 (행=SuperCluster, 열=Neuropil)
W = [
  86.08922206  10.94896344   2.552094156   0            0            0            0.070861905  0.320210561  0.015983888  0
  67.38585739  21.08198577   6.307809693   0.029043294  3.280415467  0.002461296  0.000984518  0            0.599571734  0.050210441
  40.57079211  56.58075281   1.879489091   0.002117496  0.086817319  0.000423499  0            0.001693996  0.020327958  0
   1.043310117  0.573844615 41.96665641  26.33239697  21.72288067   0            0.009620195  0            0.165467349  0.56470543
   1.841138285  1.993240753 21.15271941   0.309637166 43.75545648   0.102048339  0.871097039  0.516837976  2.723875819  7.923296899
   5.460599041  1.991077366 50.39482434   0.336910901  4.775326671  0           11.95395111   9.159132047  7.196240691  0
   0.015611828  0.028881881  0.242763918  0.003122366  0           97.3295969    0            0.601835951  0.001561183  0
];

%% 3) 임계값 이하 제거(엣지 미생성)
thr = 1.0;   % <= thr 인 값은 그래프에서 제거. 필요에 맞게 조정

src = {}; tgt = {}; wt = [];
for i = 1:size(W,1)        % SuperCluster
    for j = 1:size(W,2)    % Neuropil
        if W(i,j) > thr
            src{end+1} = SC{i}; %#ok<SAGROW>
            tgt{end+1} = NP{j}; %#ok<SAGROW>
            wt(end+1,1) = W(i,j); %#ok<SAGROW>
        end
    end
end

% 엣지가 하나도 없으면 종료
if isempty(wt)
    error('임계값 %.3g 초과 엣지가 없습니다. thr를 낮춰주세요.', thr);
end

%% 4) 유향 그래프 구성
G = digraph(src, tgt, wt);

%% 5) 노드 배치 (이분 그래프: 왼쪽=SC, 오른쪽=NP)
nSC = numel(SC); nNP = numel(NP);
nodeNames = G.Nodes.Name;
isSC = ismember(nodeNames, SC);
idxSC = find(isSC);
idxNP = find(~isSC);

% 좌우로 두 레이어 배치
X = zeros(nSC+nNP,1);
Y = zeros(nSC+nNP,1);
X(idxSC) = 0;                              % SuperCluster 레이어
Y(idxSC) = linspace(nSC,1,numel(idxSC));   % 보기 좋게 위→아래로 정렬
X(idxNP) = 1;                              % Neuropil 레이어
Y(idxNP) = linspace(nNP,1,numel(idxNP));

%% 6) 플롯
figure('Color','w');
p = plot(G, 'XData',X, 'YData',Y, 'ArrowSize',10, ...
            'NodeLabel', nodeNames);
axis off; title(sprintf('Directed bipartite graph (values > %g)', thr));

% 노드 스타일(왼쪽=파랑 원, 오른쪽=주황 네모)
highlight(p, idxSC, 'NodeColor',[0.12 0.47 0.71], 'Marker','o', 'MarkerSize',7);
highlight(p, idxNP, 'NodeColor',[0.85 0.33 0.10], 'Marker','s', 'MarkerSize',7);

% 엣지 두께를 가중치에 비례하게
w = G.Edges.Weight;
if numel(w) == 1
    lw = 2.5;
else
    wN = (w - min(w)) / (max(w) - min(w) + eps); % [0,1]
    lw = 0.8 + 4.2*wN;
end
p.LineWidth = lw;
p.EdgeColor = [0.2 0.2 0.2];

% (선택) 엣지 라벨: 값 표시 (자리수 조정)
showEdgeLabels = true;
if showEdgeLabels
    p.EdgeLabel = arrayfun(@(x) sprintf('%.2f', x), w, 'UniformOutput', false);
end

% (선택) 범례 비슷한 표시
text(-0.05, max(Y(idxSC))+0.5, 'SuperCluster →', 'FontWeight','bold');
text( 1.05, max(Y(idxNP))+0.5, 'Neuropil', 'FontWeight','bold');

%% 7) (선택) 임계값 이후의 이분 인접행렬도 보관
% 생성된 그래프의 노드 순서에 맞춘 희소행렬(옵션)
% 행: SuperCluster, 열: Neuropil만 추출
% (필요하면 사용하세요)
A = adjacency(G, 'weighted');
% 현재 G에는 SC와 NP가 섞여 있으므로, 원래 순서를 기준으로 다시 만들려면:
BiAdj = zeros(nSC, nNP);
for e = 1:numedges(G)
    i = find(strcmp(G.Edges.EndNodes{e,1}, SC));
    j = find(strcmp(G.Edges.EndNodes{e,2}, NP));
    if ~isempty(i) && ~isempty(j)
        BiAdj(i,j) = G.Edges.Weight(e);
    end
end
% disp(BiAdj);  % 확인용

%% 8) (선택) 그림 저장
% exportgraphics(gcf, sprintf('bipartite_digraph_thr_%g.pdf', thr), 'ContentType','vector');
