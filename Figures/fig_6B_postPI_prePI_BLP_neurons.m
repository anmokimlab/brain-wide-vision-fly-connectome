%% 1. Load data and initialize
clear all; close all; clc

baseDir = fileparts(fileparts(mfilename('fullpath')));  % resolved from this script's location (portable; no absolute paths)

% Per-neuron synapse / NPI table from Data_Processing/s02_compute_postPI_prePI_all_neurons.m
load(fullfile(baseDir, 'Processed_Data', 'FAFB_NPI_thr0.mat'))

%% 2. BLP (bilateral) PostPI / PrePI (right vs left optic-lobe lateralization)
FAFBNPIs.In_Synapse_Total  = FAFBNPIs.In_Synapse_Optic_R  + FAFBNPIs.In_Synapse_Central  + FAFBNPIs.In_Synapse_Optic_L;
FAFBNPIs.Out_Synapse_Total = FAFBNPIs.Out_Synapse_Optic_R + FAFBNPIs.Out_Synapse_Central + FAFBNPIs.Out_Synapse_Optic_L;

FAFBNPIs.BLP_PostPI = (FAFBNPIs.In_Synapse_Optic_R  - FAFBNPIs.In_Synapse_Optic_L)  ./ (FAFBNPIs.In_Synapse_Optic_R  + FAFBNPIs.In_Synapse_Optic_L);
FAFBNPIs.BLP_PrePI  = (FAFBNPIs.Out_Synapse_Optic_R - FAFBNPIs.Out_Synapse_Optic_L) ./ (FAFBNPIs.Out_Synapse_Optic_R + FAFBNPIs.Out_Synapse_Optic_L);

FAFBNPIs_All = FAFBNPIs;   % keep the unfiltered table for per-type neuron counts

%% 3. Select BLP (bilateral) neurons
% Condition 1: at least 5 synapses bridging the right and left optic lobes
idx = (FAFBNPIs.In_Synapse_Optic_R >= 5  & FAFBNPIs.Out_Synapse_Optic_L >= 5) | ...
      (FAFBNPIs.Out_Synapse_Optic_R >= 5 & FAFBNPIs.In_Synapse_Optic_L  >= 5);
FAFBNPIs = FAFBNPIs(idx, :);

% Condition 2: drop neurons whose input or output is mostly central (> Thr_L)
Thr_L = 0.7;
idx = (FAFBNPIs.In_Synapse_Central  > (FAFBNPIs.In_Synapse_Optic_R  + FAFBNPIs.In_Synapse_Optic_L  + FAFBNPIs.In_Synapse_Central)  * Thr_L) | ...
      (FAFBNPIs.Out_Synapse_Central > (FAFBNPIs.Out_Synapse_Optic_R + FAFBNPIs.Out_Synapse_Optic_L + FAFBNPIs.Out_Synapse_Central) * Thr_L);
FAFBNPIs(idx, :) = [];

% If a cell type has a 'right' neuron, drop its 'left' neurons
rightSideTypes = FAFBNPIs.type(strcmp(FAFBNPIs.side, 'right'));
idxToRemove = strcmp(FAFBNPIs.side, 'left') & ismember(FAFBNPIs.type, rightSideTypes);
FAFBNPIs(idxToRemove, :) = [];

%% 4. Aggregate per cell type
[UniqueTypes, ~, ~] = unique(FAFBNPIs.type);
FAFBNPIs_by_type = table(UniqueTypes, 'VariableNames', {'type'});

for i = 1:size(FAFBNPIs_by_type, 1)
    idx = strcmp(FAFBNPIs.type, FAFBNPIs_by_type.type{i});

    currentSuperclass = FAFBNPIs.superclass(idx);
    [uniqueSC, ~, ic] = unique(currentSuperclass);
    [~, maxidx] = max(accumarray(ic, 1));
    FAFBNPIs_by_type.SuperClass{i} = uniqueSC{maxidx};

    FAFBNPIs_by_type.Mean_BLP_PostPI(i) = mean(FAFBNPIs.BLP_PostPI(idx));
    FAFBNPIs_by_type.Std_BLP_PostPI(i)  = std(FAFBNPIs.BLP_PostPI(idx));
    FAFBNPIs_by_type.Mean_BLP_PrePI(i)  = mean(FAFBNPIs.BLP_PrePI(idx));
    FAFBNPIs_by_type.Std_BLP_PrePI(i)   = std(FAFBNPIs.BLP_PrePI(idx));
    FAFBNPIs_by_type.Number_of_neurons(i)       = sum(idx);
    FAFBNPIs_by_type.Total_number_of_neurons(i) = sum(strcmp(FAFBNPIs_All.type, FAFBNPIs_by_type.type{i}));
end

% Condition 3: drop types kept in fewer than 20% of their neurons
idx = (FAFBNPIs_by_type.Number_of_neurons ./ FAFBNPIs_by_type.Total_number_of_neurons * 100) < 20;
FAFBNPIs_by_type(idx, :) = [];
FAFBNPIs(~ismember(FAFBNPIs.type, FAFBNPIs_by_type.type), :) = [];

%% 5. Figure 1 (panel 6B): per-type BLP PostPI vs PrePI, colored by superclass
uniqueSuperclass = {'endocrine', 'motor', 'sensory', 'visual_projection', ...
    'visual_centrifugal', 'optic', 'central', 'descending', 'ascending'};
colors = [ 0.3010, 0.7450, 0.9330;   % endocrine
           0.2780, 0.6000, 0.8000;   % motor
           0.8500, 0.3250, 0.0980;   % sensory
           0.0000, 0.4470, 0.7410;   % visual_projection
           0.4660, 0.6740, 0.1880;   % visual_centrifugal
           0.9290, 0.6940, 0.1250;   % optic
           0.4940, 0.1840, 0.5560;   % central
           0.8500, 0.1500, 0.2000;   % descending
           0.6350, 0.5090, 0.2540 ]; % ascending

figure(1); set(gcf,'Color','w'); hold on;
plot(linspace(-1,0.8,10), linspace(-0.8,1,10), 'Color','#EAEBEB', 'LineWidth',1);
plot(linspace(-0.8,1,10), linspace(-1,0.8,10), 'Color','#EAEBEB', 'LineWidth',1);

for i = 1:numel(uniqueSuperclass)
    idx = strcmp(FAFBNPIs_by_type.SuperClass, uniqueSuperclass{i});
    scatter(FAFBNPIs_by_type.Mean_BLP_PostPI(idx), FAFBNPIs_by_type.Mean_BLP_PrePI(idx), ...
        27, colors(i,:), 'filled', 'MarkerFaceAlpha', 0.65);
end
hold off;
axis square; grid on;
xlabel('Input CenterOfMass'); ylabel('Output CenterOfMass');
set(gca,'TickDir','Out','Box','off','XTick',-1:0.2:1,'YTick',-1:0.2:1)
xlim([-1.05 1.05]); ylim([-1.05 1.05]);

%% 6. Classify BLP types into right (R) / left (L) / bidirectional (BD) and save
BLP_R_by_type = FAFBNPIs_by_type((FAFBNPIs_by_type.Mean_BLP_PostPI - FAFBNPIs_by_type.Mean_BLP_PrePI) >= 0.2, :);
BLP_L_by_type = FAFBNPIs_by_type((FAFBNPIs_by_type.Mean_BLP_PostPI - FAFBNPIs_by_type.Mean_BLP_PrePI) <= -0.2, :);
BLP_BD_type   = FAFBNPIs_by_type(((FAFBNPIs_by_type.Mean_BLP_PostPI - FAFBNPIs_by_type.Mean_BLP_PrePI) < 0.2) & ...
                                       ((FAFBNPIs_by_type.Mean_BLP_PostPI - FAFBNPIs_by_type.Mean_BLP_PrePI) > -0.2), :);

BLP_R_NPIs  = FAFBNPIs(ismember(FAFBNPIs.type, BLP_R_by_type.type), :);
BLP_L_NPIs  = FAFBNPIs(ismember(FAFBNPIs.type, BLP_L_by_type.type), :);
BLP_BD_NPIs = FAFBNPIs(ismember(FAFBNPIs.type, BLP_BD_type.type), :);

save(fullfile(baseDir, 'Processed_Data', 'BLP_neurons_thr0.mat'), ...
    'BLP_R_by_type', 'BLP_L_by_type', 'BLP_BD_type', ...
    'BLP_R_NPIs', 'BLP_L_NPIs', 'BLP_BD_NPIs');

%% Save the root_ids of BLP_R neurons as CSV (no header, root_id only)
writematrix(BLP_R_NPIs.root_id, fullfile(baseDir, 'Processed_Data', 'BLP_R_root_ids.csv'));
