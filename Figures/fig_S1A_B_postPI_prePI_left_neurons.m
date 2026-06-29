%% 1. Load data and initialize
clear all; close all; clc

baseDir = fileparts(fileparts(mfilename('fullpath')));  % resolved from this script's location (portable; no absolute paths)

% Per-neuron synapse / NPI table from Data_Processing/s02_compute_postPI_prePI_all_neurons.m
load(fullfile(baseDir, 'Processed_Data', 'FAFB_NPI_thr0.mat'))

%% 2. Left PostPI / PrePI (left optic lobe vs central)
FAFBNPIs.In_Synapse_Total  = FAFBNPIs.In_Synapse_Optic_R  + FAFBNPIs.In_Synapse_Central  + FAFBNPIs.In_Synapse_Optic_L;
FAFBNPIs.Out_Synapse_Total = FAFBNPIs.Out_Synapse_Optic_R + FAFBNPIs.Out_Synapse_Central + FAFBNPIs.Out_Synapse_Optic_L;

FAFBNPIs.Left_PostPI = (FAFBNPIs.In_Synapse_Optic_L  - FAFBNPIs.In_Synapse_Central)  ./ (FAFBNPIs.In_Synapse_Optic_L  + FAFBNPIs.In_Synapse_Central);
FAFBNPIs.Left_PrePI  = (FAFBNPIs.Out_Synapse_Optic_L - FAFBNPIs.Out_Synapse_Central) ./ (FAFBNPIs.Out_Synapse_Optic_L + FAFBNPIs.Out_Synapse_Central);

FAFBNPIs_All = FAFBNPIs;   % keep the unfiltered table for per-type neuron counts

%% 3. Select left optic-lobe neurons
% Condition 1: at least 5 synapses bridging the left optic lobe and the central brain
idx = (FAFBNPIs.In_Synapse_Optic_L >= 5  & FAFBNPIs.Out_Synapse_Central >= 5) | ...
      (FAFBNPIs.Out_Synapse_Optic_L >= 5 & FAFBNPIs.In_Synapse_Central  >= 5);
FAFBNPIs = FAFBNPIs(idx, :);

% Condition 2: drop neurons whose input or output is mostly in the right optic lobe (> Thr_R)
Thr_R = 0.7;
idx = (FAFBNPIs.In_Synapse_Optic_R  > (FAFBNPIs.In_Synapse_Optic_R  + FAFBNPIs.In_Synapse_Optic_L  + FAFBNPIs.In_Synapse_Central)  * Thr_R) | ...
      (FAFBNPIs.Out_Synapse_Optic_R > (FAFBNPIs.Out_Synapse_Optic_R + FAFBNPIs.Out_Synapse_Optic_L + FAFBNPIs.Out_Synapse_Central) * Thr_R);
FAFBNPIs(idx, :) = [];

%% 4. Aggregate per cell type
[UniqueTypes, ~, ~] = unique(FAFBNPIs.type);
FAFBNPIs_by_type = table(UniqueTypes, 'VariableNames', {'type'});
for i = 1:size(FAFBNPIs_by_type, 1)
    idx = strcmp(FAFBNPIs.type, FAFBNPIs_by_type.type{i});

    currentSuperclass = FAFBNPIs.superclass(idx);
    [uniqueSC, ~, ic] = unique(currentSuperclass);
    [~, maxidx] = max(accumarray(ic, 1));
    FAFBNPIs_by_type.SuperClass{i} = uniqueSC{maxidx};

    FAFBNPIs_by_type.Mean_Left_PostPI(i) = mean(FAFBNPIs.Left_PostPI(idx));
    FAFBNPIs_by_type.Std_Left_PostPI(i)  = std(FAFBNPIs.Left_PostPI(idx));
    FAFBNPIs_by_type.Mean_Left_PrePI(i)  = mean(FAFBNPIs.Left_PrePI(idx));
    FAFBNPIs_by_type.Std_Left_PrePI(i)   = std(FAFBNPIs.Left_PrePI(idx));
    FAFBNPIs_by_type.Number_of_neurons(i)       = sum(idx);
    FAFBNPIs_by_type.Total_number_of_neurons(i) = sum(strcmp(FAFBNPIs_All.type, FAFBNPIs_by_type.type{i}));
end

% Condition 3: drop types kept in fewer than 20% of their neurons
idx = (FAFBNPIs_by_type.Number_of_neurons ./ FAFBNPIs_by_type.Total_number_of_neurons * 100) < 20;
FAFBNPIs_by_type(idx, :) = [];
FAFBNPIs(~ismember(FAFBNPIs.type, FAFBNPIs_by_type.type), :) = [];

%% Superclass color map
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

% Histograms (panel S1B) use only these four superclasses (rows 4-7 of colors)
histSuperclass = {'visual_projection', 'visual_centrifugal', 'optic', 'central'};

%% 5. Figure 1 (panel S1A): per-type Left PostPI vs PrePI, colored by superclass
figure(1); set(gcf,'Color','w'); hold on;
plot(linspace(-1,0.8,10), linspace(-0.8,1,10), 'Color','#EAEBEB', 'LineWidth',1);
plot(linspace(-0.8,1,10), linspace(-1,0.8,10), 'Color','#EAEBEB', 'LineWidth',1);
for i = 1:numel(uniqueSuperclass)
    idx = strcmp(FAFBNPIs_by_type.SuperClass, uniqueSuperclass{i});
    scatter(FAFBNPIs_by_type.Mean_Left_PostPI(idx), FAFBNPIs_by_type.Mean_Left_PrePI(idx), ...
        27, colors(i,:), 'filled', 'MarkerFaceAlpha', 0.65);
end
hold off;
axis square; grid on;
xlabel('Input CenterOfMass'); ylabel('Output CenterOfMass');
set(gca,'TickDir','Out','Box','off','XTick',-1:0.2:1,'YTick',-1:0.2:1)
xlim([-1.05 1.05]); ylim([-1.05 1.05]);

%% 6. Figure 2 (panel S1B): Left PostPI + PrePI histogram
figure(2); set(gcf,'Color','w'); hold on;
for i = 1:numel(histSuperclass)
    idx = strcmp(FAFBNPIs.superclass, histSuperclass{i});
    if any(idx)
        h = histogram(FAFBNPIs.Left_PostPI(idx) + FAFBNPIs.Left_PrePI(idx), -2:0.4:2, 'Normalization','probability');
        h.FaceColor = colors(i+3,:); h.EdgeColor = colors(i+3,:); h.FaceAlpha = 0.8;
    end
end
hold off;
axis square
xlabel('PostPI+PrePI')
set(gca,'TickDir','Out','Box','off','XTick',-2:0.4:2,'YTick',0:0.2:1)
xlim([-2.05 2.05]); ylim([0 1]);

%% 7. Figure 3 (panel S1B): Left PostPI - PrePI histogram
figure(3); set(gcf,'Color','w'); hold on;
for i = 1:numel(histSuperclass)
    idx = strcmp(FAFBNPIs.superclass, histSuperclass{i});
    if any(idx)
        h = histogram(FAFBNPIs.Left_PostPI(idx) - FAFBNPIs.Left_PrePI(idx), -2:0.4:2, 'Normalization','probability');
        h.FaceColor = colors(i+3,:); h.EdgeColor = colors(i+3,:); h.FaceAlpha = 0.8;
    end
end
hold off;
axis square
xlabel('PostPI-PrePI')
set(gca,'TickDir','Out','Box','off','XTick',-2:0.4:2,'YTick',0:0.2:1)
xlim([-2.05 2.05]); ylim([0 1]);
