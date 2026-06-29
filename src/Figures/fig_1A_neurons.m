clear all; clc; close all

% =========================================================================
% [1. Settings] Edit the paths and styles here
% =========================================================================
baseDir = fileparts(fileparts(mfilename('fullpath')));  % resolved from this script's location (portable; no absolute paths)

% Data load settings (force root_id to int64)
typeOpts = detectImportOptions(fullfile(baseDir, 'Codex_Data', 'consolidated_cell_types.csv'));
typeOpts = setvartype(typeOpts, 'root_id', 'int64');
FAFBConsolidated_type = readtable(fullfile(baseDir, 'Codex_Data', 'consolidated_cell_types.csv'), typeOpts);

% Background (Neuropil) settings
neuropilDir   = fullfile(baseDir, 'Processed_Data', 'fig1a_brain_mesh');
neuropilColor = [0.7 0.7 0.7];
neuropilAlpha = 0.05;

% Foreground (Neuron) settings
neuronDir      = fullfile(baseDir, 'Processed_Data', 'fig1a_neuron_mesh');
neuronAlpha    = 1.0;
neuronEdge     = 'none';
useSingleColor = false;
singleColor    = [0.8 0.8 0.8];

% =========================================================================
% [2. Create figure]
% =========================================================================
figure('Color','w','WindowState','maximized', 'InvertHardcopy', 'off');
hold on;

% =========================================================================
% [Part 1] Draw the neuropil (background)
% =========================================================================
vFilesNP = dir(fullfile(neuropilDir, "*_vertices.csv"));

if ~isempty(vFilesNP)
    fprintf("Plotting Neuropil...\n");

    for i = 1:numel(vFilesNP)
        vPath = fullfile(vFilesNP(i).folder, vFilesNP(i).name);
        fPath = fullfile(vFilesNP(i).folder, ...
            strrep(vFilesNP(i).name, "_vertices.csv", "_faces.csv"));

        if ~isfile(fPath)
            continue;
        end

        V_np = readmatrix(vPath);
        F_np = readmatrix(fPath);

        if isempty(V_np) || isempty(F_np)
            continue;
        end

        F_np = round(F_np);

        if min(F_np(:)) == 0
            F_np = F_np + 1;
        end

        trisurf(F_np, V_np(:,1), V_np(:,2), V_np(:,3), ...
            'EdgeColor', 'none', ...
            'FaceColor', neuropilColor, ...
            'FaceAlpha', neuropilAlpha, ...
            'FaceLighting', 'none');
    end
end

% =========================================================================
% [Part 2] Draw the neurons (foreground) - precise int64 conversion logic
% =========================================================================
vFilesNeu = dir(fullfile(neuronDir, "*_vertices.csv"));

if isempty(vFilesNeu)

    fprintf("[WARN] No neuron files found in %s\n", neuronDir);

else

    fprintf("Mapping Types using Custom int64 Logic...\n");

    numFiles = numel(vFilesNeu);
    neuronTypes = cell(numFiles, 1);

    % ---------------------------------------------------------------------
    % Step 1: extract the ID from the file name and convert it precisely
    % ---------------------------------------------------------------------
    for i = 1:numFiles

        idStrList = regexp(vFilesNeu(i).name, '\d+', 'match');

        if ~isempty(idStrList)

            % Precise string -> int64 conversion
            s = char(idStrList{1});
            neg = (s(1) == '-');

            if neg
                s = s(2:end);
            end

            acc = int64(0);

            for k = 1:numel(s)
                d = int64(s(k) - '0');

                % overflow check
                if acc > idivide(intmax('int64') - d, int64(10))
                    error('int64 overflow detected for ID: %s', s);
                end

                acc = acc * 10 + d;
            end

            if neg
                acc = -acc;
            end

            currID = acc;  % precisely converted int64 ID

            % Table matching
            rowIdx = (FAFBConsolidated_type.root_id == currID);

            if any(rowIdx)
                t = FAFBConsolidated_type.primary_type(rowIdx);

                if iscell(t)
                    neuronTypes{i} = char(t{1});
                elseif iscategorical(t)
                    neuronTypes{i} = char(t);
                else
                    neuronTypes{i} = char(string(t));
                end
            else
                neuronTypes{i} = 'Unknown';
            end

        else
            neuronTypes{i} = 'Unknown';
        end
    end

    % ---------------------------------------------------------------------
    % Step 2: assign a color per type
    % ---------------------------------------------------------------------
    uniqueTypes = unique(neuronTypes);
    numUnique = numel(uniqueTypes);

    % Fixed per-type colors (blue / green / red)
    %   LT51 = blue, VCH = green, LC14 = red
    blue  = [0      0.4470 0.7410];
    green = [0.4660 0.6740 0.1880];
    red   = [0.8500 0.3250 0.0980];
    fixedColors = containers.Map( ...
        {'LT51', 'VCH', 'LC14'}, ...
        {blue,   green,  red});

    % Fallback color order for any other types (gray, then the rest, repeating)
    customColors = [
        114/255 113/255 113/255;
        0.8500  0.3250  0.0980;
        0       0.4470  0.7410;
        0.4660 0.6740 0.1880
    ];

    typeColorMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
    fallbackIdx = 0;

    for i = 1:numUnique
        tname = uniqueTypes{i};

        if isKey(fixedColors, tname)
            typeColorMap(tname) = fixedColors(tname);
        else
            fallbackIdx = fallbackIdx + 1;
            row = mod(fallbackIdx - 1, size(customColors, 1)) + 1;
            typeColorMap(tname) = customColors(row, :);
        end
    end

    % ---------------------------------------------------------------------
    % Step 3: actually draw
    % ---------------------------------------------------------------------
    for i = 1:numFiles

        vPath = fullfile(vFilesNeu(i).folder, vFilesNeu(i).name);
        fPath = fullfile(vFilesNeu(i).folder, ...
            strrep(vFilesNeu(i).name, "_vertices.csv", "_faces.csv"));

        if ~isfile(fPath)
            continue;
        end

        V = readmatrix(vPath);
        F = readmatrix(fPath);

        if isempty(V) || isempty(F)
            continue;
        end

        V = double(V(:,1:3));
        F = round(double(F(:,1:3)));

        if min(F(:)) == 0
            F = F + 1;
        end

        % Index validity check
        maxIdx = size(V,1);
        valid = all(F >= 1 & F <= maxIdx, 2);

        if ~all(valid)
            F = F(valid,:);
        end

        % Color decision
        fc = singleColor;

        if ~useSingleColor
            fc = typeColorMap(neuronTypes{i});
        end

        trisurf(F, V(:,1), V(:,2), V(:,3), ...
            'FaceColor', fc, ...
            'EdgeColor', neuronEdge, ...
            'FaceAlpha', neuronAlpha, ...
            'DisplayName', neuronTypes{i});
    end
end

% =========================================================================
% [Part 3] Finalize and save
% =========================================================================
view(0, -90);
axis equal;
axis tight;
axis off;

delete(findall(gcf, 'Type', 'light'));
camlight('headlight');
lighting gouraud;
material shiny;

set(gca, 'Units', 'normalized', 'Position', [0 0 1 1]);
