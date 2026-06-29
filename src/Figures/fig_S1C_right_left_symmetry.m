%% 1. Load data and initialize
clear all; close all; clc

baseDir = fileparts(fileparts(mfilename('fullpath')));  % resolved from this script's location (portable; no absolute paths)

% Per-neuron synapse / NPI table from Data_Processing/s02_compute_postPI_prePI_all_neurons.m
load(fullfile(baseDir, 'Processed_Data', 'FAFB_NPI_thr0.mat'))

% Right / Left PostPI / PrePI (optic lobe vs central, per hemisphere)
FAFBNPIs.In_Synapse_Total  = FAFBNPIs.In_Synapse_Optic_R  + FAFBNPIs.In_Synapse_Central  + FAFBNPIs.In_Synapse_Optic_L;
FAFBNPIs.Out_Synapse_Total = FAFBNPIs.Out_Synapse_Optic_R + FAFBNPIs.Out_Synapse_Central + FAFBNPIs.Out_Synapse_Optic_L;

FAFBNPIs.Right_PostPI = (FAFBNPIs.In_Synapse_Optic_R  - FAFBNPIs.In_Synapse_Central)  ./ (FAFBNPIs.In_Synapse_Optic_R  + FAFBNPIs.In_Synapse_Central);
FAFBNPIs.Right_PrePI  = (FAFBNPIs.Out_Synapse_Optic_R - FAFBNPIs.Out_Synapse_Central) ./ (FAFBNPIs.Out_Synapse_Optic_R + FAFBNPIs.Out_Synapse_Central);
FAFBNPIs.Left_PostPI  = (FAFBNPIs.In_Synapse_Optic_L  - FAFBNPIs.In_Synapse_Central)  ./ (FAFBNPIs.In_Synapse_Optic_L  + FAFBNPIs.In_Synapse_Central);
FAFBNPIs.Left_PrePI   = (FAFBNPIs.Out_Synapse_Optic_L - FAFBNPIs.Out_Synapse_Central) ./ (FAFBNPIs.Out_Synapse_Optic_L + FAFBNPIs.Out_Synapse_Central);

FAFBNPIs_All = FAFBNPIs;   % unfiltered table (for per-type total neuron counts)

Thr = 0.7;

% Right optic-lobe neurons (conditions 1-3, as in fig_1D_E_postPI_prePI_right_neurons.m)
Right_Neurons_All = FAFBNPIs;
idx = (Right_Neurons_All.In_Synapse_Optic_R >= 5  & Right_Neurons_All.Out_Synapse_Central >= 5) | ...
      (Right_Neurons_All.Out_Synapse_Optic_R >= 5 & Right_Neurons_All.In_Synapse_Central  >= 5);
Right_Neurons_All = Right_Neurons_All(idx, :);
idx = (Right_Neurons_All.In_Synapse_Optic_L  > (Right_Neurons_All.In_Synapse_Optic_R  + Right_Neurons_All.In_Synapse_Optic_L  + Right_Neurons_All.In_Synapse_Central)  * Thr) | ...
      (Right_Neurons_All.Out_Synapse_Optic_L > (Right_Neurons_All.Out_Synapse_Optic_R + Right_Neurons_All.Out_Synapse_Optic_L + Right_Neurons_All.Out_Synapse_Central) * Thr);
Right_Neurons_All(idx, :) = [];
Right_Neurons_All = drop_small_types(Right_Neurons_All, FAFBNPIs_All, 20);   % drop types kept in <20% of their neurons

% Left optic-lobe neurons (conditions 1-3, as in fig_S1A_B_postPI_prePI_left_neurons.m)
Left_Neurons_All = FAFBNPIs;
idx = (Left_Neurons_All.In_Synapse_Optic_L >= 5  & Left_Neurons_All.Out_Synapse_Central >= 5) | ...
      (Left_Neurons_All.Out_Synapse_Optic_L >= 5 & Left_Neurons_All.In_Synapse_Central  >= 5);
Left_Neurons_All = Left_Neurons_All(idx, :);
idx = (Left_Neurons_All.In_Synapse_Optic_R  > (Left_Neurons_All.In_Synapse_Optic_R  + Left_Neurons_All.In_Synapse_Optic_L  + Left_Neurons_All.In_Synapse_Central)  * Thr) | ...
      (Left_Neurons_All.Out_Synapse_Optic_R > (Left_Neurons_All.Out_Synapse_Optic_R + Left_Neurons_All.Out_Synapse_Optic_L + Left_Neurons_All.Out_Synapse_Central) * Thr);
Left_Neurons_All(idx, :) = [];
Left_Neurons_All = drop_small_types(Left_Neurons_All, FAFBNPIs_All, 20);     % drop types kept in <20% of their neurons

%% 2. Per-type means (right)
[UniqueTypesRight,~,~] = unique(Right_Neurons_All.type);
Right_by_type = table(UniqueTypesRight, 'VariableNames', {'type'});
for i = 1:size(Right_by_type,1)
    idx = strcmp(Right_Neurons_All.type, Right_by_type.type{i});
    [uniqueSC,~,ic] = unique(Right_Neurons_All.superclass(idx));
    [~,maxidx] = max(accumarray(ic,1));
    Right_by_type.SuperClass{i}        = uniqueSC{maxidx};
    Right_by_type.Mean_Right_PostPI{i} = mean(Right_Neurons_All.Right_PostPI(idx), 'omitmissing');
    Right_by_type.Std_Right_PostPI{i}  = std(Right_Neurons_All.Right_PostPI(idx),  'omitmissing');
    Right_by_type.Mean_Right_PrePI{i}  = mean(Right_Neurons_All.Right_PrePI(idx),  'omitmissing');
    Right_by_type.Std_Right_PrePI{i}   = std(Right_Neurons_All.Right_PrePI(idx),   'omitmissing');
    Right_by_type.Number_of_neurons(i)       = sum(idx);
    Right_by_type.Total_number_of_neurons(i) = sum(strcmp(FAFBNPIs_All.type, Right_by_type.type{i}));
end

%% 3. Per-type means (left)
[UniqueTypesLeft,~,~] = unique(Left_Neurons_All.type);
Left_by_type = table(UniqueTypesLeft, 'VariableNames', {'type'});
for i = 1:size(Left_by_type,1)
    idx = strcmp(Left_Neurons_All.type, Left_by_type.type{i});
    [uniqueSC,~,ic] = unique(Left_Neurons_All.superclass(idx));
    [~,maxidx] = max(accumarray(ic,1));
    Left_by_type.SuperClass{i}       = uniqueSC{maxidx};
    Left_by_type.Mean_Left_PostPI{i} = mean(Left_Neurons_All.Left_PostPI(idx), 'omitmissing');
    Left_by_type.Std_Left_PostPI{i}  = std(Left_Neurons_All.Left_PostPI(idx),  'omitmissing');
    Left_by_type.Mean_Left_PrePI{i}  = mean(Left_Neurons_All.Left_PrePI(idx),  'omitmissing');
    Left_by_type.Std_Left_PrePI{i}   = std(Left_Neurons_All.Left_PrePI(idx),   'omitmissing');
    Left_by_type.Number_of_neurons(i)       = sum(idx);
    Left_by_type.Total_number_of_neurons(i) = sum(strcmp(FAFBNPIs_All.type, Left_by_type.type{i}));
end

%% 4. Right-vs-left distance per type (types present on both sides)
idx = ismember(Left_by_type.type, Right_by_type.type);
uniqueTypes = table(Left_by_type.type(idx), 'VariableNames', {'type'});
uniqueTypes.SuperClass = Left_by_type.SuperClass(idx);
for i = 1:size(uniqueTypes,1)
    idx_L = strcmp(Left_by_type.type,  uniqueTypes.type{i});
    idx_R = strcmp(Right_by_type.type, uniqueTypes.type{i});
    uniqueTypes.Distance(i) = sqrt( (Left_by_type.Mean_Left_PostPI{idx_L} - Right_by_type.Mean_Right_PostPI{idx_R})^2 + ...
                                    (Left_by_type.Mean_Left_PrePI{idx_L}  - Right_by_type.Mean_Right_PrePI{idx_R})^2 );
end

% Group distances by superclass
VPN_Table     = uniqueTypes(strcmp(uniqueTypes.SuperClass,'visual_projection'), :);
VCN_Table     = uniqueTypes(strcmp(uniqueTypes.SuperClass,'visual_centrifugal'), :);
Optic_Table   = uniqueTypes(strcmp(uniqueTypes.SuperClass,'optic'), :);
Central_Table = uniqueTypes(strcmp(uniqueTypes.SuperClass,'central'), :);

%% 5. Figure 1 (panel S1C): right-vs-left distance per superclass
colors = [ 0.0000, 0.4470, 0.7410;   % VPN
           0.4660, 0.6740, 0.1880;   % VCN
           0.9290, 0.6940, 0.1250;   % Optic
           0.4940, 0.1840, 0.5560 ]; % Central
scatter_colors = colors * 0.8;

figure(1); clf; set(gcf,'Color','w'); hold on;

group_data = {VPN_Table.Distance, VCN_Table.Distance, Optic_Table.Distance, Central_Table.Distance};
data_means = cellfun(@mean, group_data);
data_stds  = cellfun(@std,  group_data);
x_labels   = ["VPN", "VCN", "Optic", "Central"];

% Bar chart
bar_handle = bar(x_labels, data_means, 'FaceColor', 'flat', 'EdgeColor', 'none');
for i = 1:numel(data_means)
    bar_handle.CData(i, :) = colors(i, :);
end
xtick_positions = bar_handle.XEndPoints;

% Error bars
errorbar(xtick_positions, data_means, data_stds, 'k.', 'LineWidth', 1.5);

% Individual data points (jittered)
jitter_width = 0.15;
for i = 1:numel(group_data)
    n_points = numel(group_data{i});
    jittered_x = xtick_positions(i) + (rand(n_points,1) - 0.5) * jitter_width;
    scatter(jittered_x, group_data{i}, 15, scatter_colors(i,:), 'filled', 'MarkerFaceAlpha', 0.6);
end

set(gca, 'Box','off', 'TickDir','out', 'FontSize', 12);
ylim([0 2*sqrt(2)]);
ylabel('Distance (L-R pair)');
title('Mean Distance with Individual Data Points');
grid on;

%% ===================== Local functions =====================
function T = drop_small_types(T, T_All, minPct)
% Keep only cell types whose retained fraction (count in T / count in T_All) is >= minPct (%).
types = unique(T.type);
keep = false(height(T),1);
for k = 1:numel(types)
    nKept  = sum(strcmp(T.type, types{k}));
    nTotal = sum(strcmp(T_All.type, types{k}));
    if nTotal > 0 && (nKept/nTotal*100) >= minPct
        keep = keep | strcmp(T.type, types{k});
    end
end
T = T(keep, :);
end
