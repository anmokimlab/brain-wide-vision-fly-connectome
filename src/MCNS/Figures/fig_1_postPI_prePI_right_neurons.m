%% fig_1_postPI_prePI_right_neurons (MCNS)
% MCNS analogue of the FAFB Figures/fig_1D_E_postPI_prePI_right_neurons.m.
clear all; close all; clc

%% 1. Load data and initialize

baseDir = fileparts(fileparts(mfilename('fullpath')));  % resolved from this script's location (portable; no absolute paths)

% Per-neuron synapse / NPI table from MCNS/Data_Processing/s05_compute_postPI_prePI_all_neurons.m
load(fullfile(baseDir, 'Processed_Data', 'MCNS_NPI_thr0.mat'))   % MCNSNPIs

% Label neurons with no superclass as 'undecided'
MCNSNPIs.superclass(cellfun(@isempty, MCNSNPIs.superclass)) = {'undecided'};

%% 2. Right / Left PostPI / PrePI (optic lobe vs central)
MCNSNPIs.In_Synapse_Total  = MCNSNPIs.In_Synapse_Optic_R  + MCNSNPIs.In_Synapse_Central  + MCNSNPIs.In_Synapse_Optic_L;
MCNSNPIs.Out_Synapse_Total = MCNSNPIs.Out_Synapse_Optic_R + MCNSNPIs.Out_Synapse_Central + MCNSNPIs.Out_Synapse_Optic_L;

MCNSNPIs.Right_PostPI = (MCNSNPIs.In_Synapse_Optic_R  - MCNSNPIs.In_Synapse_Central)  ./ (MCNSNPIs.In_Synapse_Optic_R  + MCNSNPIs.In_Synapse_Central);
MCNSNPIs.Right_PrePI  = (MCNSNPIs.Out_Synapse_Optic_R - MCNSNPIs.Out_Synapse_Central) ./ (MCNSNPIs.Out_Synapse_Optic_R + MCNSNPIs.Out_Synapse_Central);
MCNSNPIs.Left_PostPI  = (MCNSNPIs.In_Synapse_Optic_L  - MCNSNPIs.In_Synapse_Central)  ./ (MCNSNPIs.In_Synapse_Optic_L  + MCNSNPIs.In_Synapse_Central);
MCNSNPIs.Left_PrePI   = (MCNSNPIs.Out_Synapse_Optic_L - MCNSNPIs.Out_Synapse_Central) ./ (MCNSNPIs.Out_Synapse_Optic_L + MCNSNPIs.Out_Synapse_Central);

MCNSNPIs_All = MCNSNPIs;   % keep the unfiltered table for per-type neuron counts

%% 3. Select right optic-lobe neurons
% Condition 1: at least 5 synapses bridging the right optic lobe and central brain
idx = (MCNSNPIs.In_Synapse_Optic_R >= 5  & MCNSNPIs.Out_Synapse_Central >= 5) | ...
      (MCNSNPIs.Out_Synapse_Optic_R >= 5 & MCNSNPIs.In_Synapse_Central  >= 5);
MCNSNPIs = MCNSNPIs(idx, :);

% Condition 2: drop neurons whose input or output is mostly in the left optic lobe (> Thr_L)
Thr_L = 0.7;
idx = (MCNSNPIs.In_Synapse_Optic_L  > (MCNSNPIs.In_Synapse_Optic_R  + MCNSNPIs.In_Synapse_Optic_L  + MCNSNPIs.In_Synapse_Central)  * Thr_L) | ...
      (MCNSNPIs.Out_Synapse_Optic_L > (MCNSNPIs.Out_Synapse_Optic_R + MCNSNPIs.Out_Synapse_Optic_L + MCNSNPIs.Out_Synapse_Central) * Thr_L);
MCNSNPIs(idx, :) = [];

%% 4. Aggregate per cell type
[UniqueTypes,~,~] = unique(MCNSNPIs.type);
MCNSNPIs_by_type = table(UniqueTypes, 'VariableNames', {'type'});
for i = 1:size(MCNSNPIs_by_type,1)
    idx = strcmp(MCNSNPIs.type, MCNSNPIs_by_type.type{i});
    [uniqueSC,~,ic] = unique(MCNSNPIs.superclass(idx));
    [~,maxidx] = max(accumarray(ic,1));
    MCNSNPIs_by_type.SuperClass{i}        = uniqueSC{maxidx};
    MCNSNPIs_by_type.Mean_Right_PostPI(i) = mean(MCNSNPIs.Right_PostPI(idx));
    MCNSNPIs_by_type.Std_Right_PostPI(i)  = std(MCNSNPIs.Right_PostPI(idx));
    MCNSNPIs_by_type.Mean_Right_PrePI(i)  = mean(MCNSNPIs.Right_PrePI(idx));
    MCNSNPIs_by_type.Std_Right_PrePI(i)   = std(MCNSNPIs.Right_PrePI(idx));
    MCNSNPIs_by_type.Number_of_neurons(i)       = sum(idx);
    MCNSNPIs_by_type.Total_number_of_neurons(i) = sum(strcmp(MCNSNPIs_All.type, MCNSNPIs_by_type.type{i}));
end

% Condition 3: drop types kept in fewer than 20% of their neurons
idx = (MCNSNPIs_by_type.Number_of_neurons ./ MCNSNPIs_by_type.Total_number_of_neurons * 100) < 20;
MCNSNPIs_by_type(idx, :) = [];
MCNSNPIs(~ismember(MCNSNPIs.type, MCNSNPIs_by_type.type), :) = [];

%% Superclass color map
uniqueSuperclass = {'endocrine', 'motor', 'sensory', 'visual_projection', ...
    'visual_centrifugal', 'optic', 'central', 'descending', 'ascending', 'undecided'};
colors = [ 0.3010, 0.7450, 0.9330;   % endocrine
           0.2780, 0.6000, 0.8000;   % motor
           0.8500, 0.3250, 0.0980;   % sensory
           0.0000, 0.4470, 0.7410;   % visual_projection
           0.4660, 0.6740, 0.1880;   % visual_centrifugal
           0.9290, 0.6940, 0.1250;   % optic
           0.4940, 0.1840, 0.5560;   % central
           0.8500, 0.1500, 0.2000;   % descending
           0.6350, 0.5090, 0.2540;   % ascending
           0.5000, 0.5000, 0.5000 ]; % undecided

% Histograms use only these four superclasses (rows 4-7 of colors)
histSuperclass = {'visual_projection', 'visual_centrifugal', 'optic', 'central'};

%% 5. Figure 1: per-type Right PostPI vs PrePI, colored by superclass
figure(1); set(gcf,'Color','w'); hold on;
plot(linspace(-1,0.8,10), linspace(-0.8,1,10), 'Color','#EAEBEB', 'LineWidth',1);
plot(linspace(-0.8,1,10), linspace(-1,0.8,10), 'Color','#EAEBEB', 'LineWidth',1);
for i = 1:numel(uniqueSuperclass)
    idx = strcmp(MCNSNPIs_by_type.SuperClass, uniqueSuperclass{i});
    scatter(MCNSNPIs_by_type.Mean_Right_PostPI(idx), MCNSNPIs_by_type.Mean_Right_PrePI(idx), ...
        27, colors(i,:), 'filled', 'MarkerFaceAlpha', 0.65);
end
hold off;
axis square; grid on;
xlabel('postPI'); ylabel('prePI');
set(gca,'TickDir','Out','Box','off','XTick',-1:0.2:1,'YTick',-1:0.2:1)
xlim([-1.05 1.05]); ylim([-1.05 1.05]);

%% 6. Figure 2: Right PostPI + PrePI histogram
figure(2); set(gcf,'Color','w'); hold on;
for i = 1:numel(histSuperclass)
    idx = strcmp(MCNSNPIs.superclass, histSuperclass{i});
    if any(idx)
        h = histogram(MCNSNPIs.Right_PostPI(idx) + MCNSNPIs.Right_PrePI(idx), -2:0.4:2, 'Normalization','probability');
        h.FaceColor = colors(i+3,:); h.EdgeColor = colors(i+3,:); h.FaceAlpha = 0.8;
    end
end
hold off;
axis square
xlabel('PostPI+PrePI')
set(gca,'TickDir','Out','Box','off','XTick',-2:0.4:2,'YTick',0:0.2:1)
xlim([-2.05 2.05]); ylim([0 1]);

%% 7. Figure 3: Right PostPI - PrePI histogram
figure(3); set(gcf,'Color','w'); hold on;
for i = 1:numel(histSuperclass)
    idx = strcmp(MCNSNPIs.superclass, histSuperclass{i});
    if any(idx)
        h = histogram(MCNSNPIs.Right_PostPI(idx) - MCNSNPIs.Right_PrePI(idx), -2:0.4:2, 'Normalization','probability');
        h.FaceColor = colors(i+3,:); h.EdgeColor = colors(i+3,:); h.FaceAlpha = 0.8;
    end
end
hold off;
axis square
xlabel('PostPI-PrePI')
set(gca,'TickDir','Out','Box','off','XTick',-2:0.4:2,'YTick',0:0.2:1)
xlim([-2.05 2.05]); ylim([0 1]);

%% 8. Classify right-hemisphere types (FFP / FBP / BDP) by postPI - prePI
RightFFP_by_type = MCNSNPIs_by_type((MCNSNPIs_by_type.Mean_Right_PostPI - MCNSNPIs_by_type.Mean_Right_PrePI) >= 0.2, :);
RightFBP_by_type = MCNSNPIs_by_type((MCNSNPIs_by_type.Mean_Right_PostPI - MCNSNPIs_by_type.Mean_Right_PrePI) <= -0.2, :);
RightBDP_by_type = MCNSNPIs_by_type(((MCNSNPIs_by_type.Mean_Right_PostPI - MCNSNPIs_by_type.Mean_Right_PrePI) < 0.2) & ...
                                    ((MCNSNPIs_by_type.Mean_Right_PostPI - MCNSNPIs_by_type.Mean_Right_PrePI) > -0.2), :);

RightFFP_NPIs = MCNSNPIs(ismember(MCNSNPIs.type, RightFFP_by_type.type), :);
RightFBP_NPIs = MCNSNPIs(ismember(MCNSNPIs.type, RightFBP_by_type.type), :);
RightBDP_NPIs = MCNSNPIs(ismember(MCNSNPIs.type, RightBDP_by_type.type), :);

% Split the bidirectional (BDP) types into optic / central / real by postPI + prePI
RightBDP_optic_type   = RightBDP_by_type((RightBDP_by_type.Mean_Right_PostPI + RightBDP_by_type.Mean_Right_PrePI) >= 1.2, :);
RightBDP_central_type = RightBDP_by_type((RightBDP_by_type.Mean_Right_PostPI + RightBDP_by_type.Mean_Right_PrePI) <= -1.2, :);
RightBDP_real_type    = RightBDP_by_type(((RightBDP_by_type.Mean_Right_PostPI + RightBDP_by_type.Mean_Right_PrePI) < 1.2) & ...
                                         ((RightBDP_by_type.Mean_Right_PostPI + RightBDP_by_type.Mean_Right_PrePI) > -1.2), :);

RightBDP_optic_NPIs   = RightBDP_NPIs(ismember(RightBDP_NPIs.type, RightBDP_optic_type.type), :);
RightBDP_central_NPIs = RightBDP_NPIs(ismember(RightBDP_NPIs.type, RightBDP_central_type.type), :);
RightBDP_real_NPIs    = RightBDP_NPIs(ismember(RightBDP_NPIs.type, RightBDP_real_type.type), :);

%% Save right-hemisphere neuron classification
save(fullfile(baseDir, 'Processed_Data', 'right_neurons_thr0.mat'), ...
    'RightFFP_by_type','RightFFP_NPIs', ...
    'RightFBP_by_type','RightFBP_NPIs', ...
    'RightBDP_optic_type','RightBDP_optic_NPIs', ...
    'RightBDP_central_type','RightBDP_central_NPIs', ...
    'RightBDP_real_type','RightBDP_real_NPIs');

%% Save the root_ids of FFP / FBP / BDP / Others neurons as CSV (no header, root_id only)
% Others = other bidirectional neurons (optic + central BD, i.e. everything but real BD)
writematrix(RightFFP_NPIs.root_id, fullfile(baseDir, 'Processed_Data', 'right_FFP_root_ids.csv'));
writematrix(RightFBP_NPIs.root_id, fullfile(baseDir, 'Processed_Data', 'right_FBP_root_ids.csv'));
writematrix(RightBDP_real_NPIs.root_id, fullfile(baseDir, 'Processed_Data', 'right_BDP_root_ids.csv'));
writematrix([RightBDP_optic_NPIs.root_id; RightBDP_central_NPIs.root_id], ...
    fullfile(baseDir, 'Processed_Data', 'right_Others_root_ids.csv'));
