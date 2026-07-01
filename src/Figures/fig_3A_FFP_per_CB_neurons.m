%% 1. Load data and initial setup
clear all; close all; clc;

baseDir = fileparts(fileparts(mfilename('fullpath')));  % resolved from this script's location (portable; no absolute paths)

% Load Codex data
opt = detectImportOptions(fullfile(baseDir, 'Codex_Data', 'connections_princeton.csv'));
opt = setvartype(opt, {'pre_root_id', 'post_root_id'}, 'int64');
FAFBConnections = readtable(fullfile(baseDir, 'Codex_Data', 'connections_princeton.csv'), opt);

opt = detectImportOptions(fullfile(baseDir, 'Codex_Data', 'consolidated_cell_types.csv'));
opt = setvartype(opt, 'root_id', 'int64');
FAFBConsolidatedTypes = readtable(fullfile(baseDir, 'Codex_Data', 'consolidated_cell_types.csv'), opt);

% FFP / FBP neuron classification (provides RightFFP_NPIs), saved by
% Figures/fig_1D_E_postPI_prePI_right_neurons.m
load(fullfile(baseDir, 'Processed_Data', 'right_neurons_thr0.mat'));

% Central neuropils = all neuropils except the right/left optic lobe and UNASGD
FAFBNeuropils = unique(FAFBConnections.neuropil);
FAFBNeuropil_OpticLobeRight = {'LA_R';'AME_R';'ME_R';'LO_R';'LOP_R'};
FAFBNeuropil_OpticLobeLeft  = {'LA_L';'AME_L';'ME_L';'LO_L';'LOP_L'};
FAFBNeuropil_Central = FAFBNeuropils(~ismember(FAFBNeuropils, FAFBNeuropil_OpticLobeRight));
FAFBNeuropil_Central = FAFBNeuropil_Central(~ismember(FAFBNeuropil_Central, FAFBNeuropil_OpticLobeLeft));
FAFBNeuropil_Central = FAFBNeuropil_Central(~ismember(FAFBNeuropil_Central, 'UNASGD'));

%% 2. Data processing (FFP -> central brain connections)

% Keep connections whose presynaptic neuron is an FFP neuron, drop those
% whose postsynaptic neuron is also an FFP neuron, and restrict to central
% neuropils.
FFP_C_idx = ismember(FAFBConnections.pre_root_id, RightFFP_NPIs.root_id);
FFP_Central_sum = FAFBConnections(FFP_C_idx, :);
FFP_Central_sum(ismember(FFP_Central_sum.post_root_id, RightFFP_NPIs.root_id), :) = [];
FFP_Central_sum(~ismember(FFP_Central_sum.neuropil, FAFBNeuropil_Central), :) = [];
FFP_Central_sum = groupsummary(FFP_Central_sum, {'pre_root_id', 'post_root_id'}, 'sum', 'syn_count');

%% 3. Type mapping (presynaptic FFP cell type)

% Map presynaptic FFP neurons to their consolidated cell type.
[tf, loc] = ismember(FFP_Central_sum.pre_root_id, FAFBConsolidatedTypes.root_id);
FFP_Central_sum.pre_type = repmat("Unidentified", height(FFP_Central_sum), 1);
FFP_Central_sum.pre_type(tf) = string(FAFBConsolidatedTypes.primary_type(loc(tf)));

%% 4. Visualization

% Highlight color used for the FFP-type distributions (purple)
highlightColor = [0.49, 0.18, 0.56];

%% --- Figure 1: number of FFP cell types converging onto each CB neuron ---
% Count, for every central-brain (CB) postsynaptic neuron, how many distinct
% FFP presynaptic cell types provide input.

% 1. Group by central neuron (post_root_id) and presynaptic FFP type (pre_type)
Type_to_Neuron_Table = groupsummary(FFP_Central_sum, {'post_root_id', 'pre_type'});

% 2. Count the number of distinct pre_type per post_root_id (CB neuron)
Type_to_Neuron_Count = groupsummary(Type_to_Neuron_Table, 'post_root_id');

figure(1); set(gcf,'Color','w'); hold on;
histogram(Type_to_Neuron_Count.GroupCount, 'BinWidth', 1, 'Normalization', 'count', ...
          'DisplayStyle', 'bar', 'LineWidth', 2, 'EdgeColor', highlightColor);
xlim([0 40]);
set(gca, 'XScale', 'linear', 'TickDir', 'out'); grid on;
xlabel('Number of FFP Types');
ylabel('Number of CB neurons');
title('Count (Type-to-Neuron Level)');


%% --- Figure 2: convergence of FFP cell types onto CB neurons (3 groups) ---
counts = Type_to_Neuron_Count.GroupCount;

% Compute fractions (1, 2-10, >10)
n_total  = numel(counts);
n_single = sum(counts == 1);
n_mid    = sum((counts >= 2) & (counts <= 10));
n_high   = sum(counts > 10);

fractions = [n_single, n_mid, n_high] / n_total;

% Plot
figure(2); clf;
set(gcf,'Color','w');

b = bar(fractions, 'FaceColor', 'flat');

% Colors
b.CData(1,:) = [0.7 0.7 0.7];   % 1 FFP type (gray)
b.CData(2,:) = highlightColor;  % 2-10 (highlight)
b.CData(3,:) = [0.3 0.3 0.3];   % >10 (dark gray)

ylim([0 1]);
set(gca, 'TickDir', 'out', ...
         'XTick', 1:3, ...
         'Box','off',...
         'XTickLabel', {'1 FFP cell type','2-10 FFP cell types','>10 FFP cell types'});
ylabel('Fraction of CB postsynaptic neurons');

grid off;
title('Convergence of FFP cell types onto CB neurons');

