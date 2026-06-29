function plotStackedSelectedBar(dataList, selectedLabels, colors)
    % dataList: containers.Map 객체들의 cell array 또는 구조체 배열
    % selectedLabels: 관심 있는 label들의 cell array
    % colors: 각 selected label에 대한 RGB 색상 (n x 3 matrix)

    numSamples = numel(dataList);
    numSelected = numel(selectedLabels);
    stackedData = zeros(numSamples, numSelected + 1);  % +1 for "Others"

    allKeys = dataList{1}.keys;  % assume all maps have same keys

    % 데이터 채우기
    for i = 1:numSamples
        if isstruct(dataList)
            dataMap = dataList(i).data;
        else
            dataMap = dataList{i};
        end
        otherSum = 0;
        for k = 1:numel(allKeys)
            label = allKeys{k};
            value = dataMap(label);
            selIdx = find(strcmp(selectedLabels, label));
            if ~isempty(selIdx)
                stackedData(i, selIdx) = value;
            else
                otherSum = otherSum + value;
            end
        end
        stackedData(i, end) = otherSum;  % 마지막 열에 Others 저장
    end

    % 색상 정의 (선택된 색 + Others 색)
    othersColor = [137 137 137]/255;
    allColors = [colors; othersColor];

    % 그래프 그리기
    figure; set(gcf,'Color','w');
    b = bar(stackedData, 'stacked', 'FaceColor', 'flat','EdgeColor','flat');
    for j = 1:(numSelected + 1)
        b(j).FaceColor = allColors(j, :);
    end
    set(gca,'TickDir','out','Box','off');
    xlabel('Sample Index');
    ylabel('Value');
    legend([selectedLabels, {'Others'}], 'Location', 'northeastoutside');
    % title('Stacked Bar Graph with Others');
    ylim([0 1])
end
