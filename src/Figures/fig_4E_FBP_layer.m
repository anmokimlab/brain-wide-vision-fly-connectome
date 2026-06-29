clear all; close all; clc;

baseDir = fileparts(fileparts(mfilename('fullpath')));  % resolved from this script's location (portable; no absolute paths)

% FBP output-synapse depth profiles (FBP_Synapse_*), saved by
% Data_Processing/s11_FBP_output_layer.m
load(fullfile(baseDir, 'Processed_Data', 'FBP_output_synapse_layer.mat'));

% FBP upstream-synapse depth profiles (FBP_Upstream_Synapse_*), saved by
% Data_Processing/s12_FBP_upstream_layer.m
load(fullfile(baseDir, 'Processed_Data', 'FBP_upstream_synapse_layer.mat'));

% FFP / FBP / BDP neuron classification (provides RightFBP_NPIs), saved by
% Figures/fig_1D_E_postPI_prePI_right_neurons.m
load(fullfile(baseDir, 'Processed_Data', 'right_neurons_thr0.mat'),'RightFBP_NPIs')

%% Group the FBP neurons by cell type
[RightFBP_type,~,ic]=unique(RightFBP_NPIs.type);
RightFBP_type=table(RightFBP_type,'VariableNames',{'type'});
for i=1:1:size(RightFBP_type,1)
    idx=ic==i;
    RightFBP_type.root_id{i}=RightFBP_NPIs.root_id(idx);
end

% Group FBP types by their target optic lobe (Me / Lo / Lop / Multi).
% (index lists into RightFBP_type, from the FBP output-neuropil classification, s10)
Me_FBP_idx=[4 22 25 26 29 30 31 32 33 68 69 70 71 72 73 74 75 76 78 81 82 83 84 85 86 87 93];
Lo_FBP_idx=[2 10 11 12 13 14 15 16 17 18 19 21 34 35 36 37 38 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 95];
Lop_FBP_idx=[3 7 8 9 28 62 63 64 65 66];
Multi_FBP_idx=[1 5 6 23 24 27 39 60 61 67 77 79 80 88 89 90 91 92 94 96];

FBP_Me=RightFBP_type(Me_FBP_idx,:);
FBP_Lo=RightFBP_type(Lo_FBP_idx,:);
FBP_Lop=RightFBP_type(Lop_FBP_idx,:);
FBP_Multi=RightFBP_type(Multi_FBP_idx,:);

%% Custom colormaps (each column of the triplet gets its own color ramp)
n = 256;
startColor = [0, 0, 0];

% Out (orange) — the FBP neuron's own output synapses
endColorOut = [0.8500, 0.3250, 0.0980];
customColormapOut = zeros(n,3);
for i = 1:3
    customColormapOut(:,i) = linspace(startColor(i), endColorOut(i), n);
end

% In - ipsi (green) — upstream input on the right (ipsilateral) side
endColorIn_ipsi = [0.4660, 0.6740, 0.1880];
customColormapIn_ipsi = zeros(n,3);
for i = 1:3
    customColormapIn_ipsi(:,i) = linspace(startColor(i), endColorIn_ipsi(i), n);
end

% In - contra (blue) — upstream input on the left (contralateral) side
endColorIn = [0, 0.4470, 0.7410];
customColormapIn = zeros(n,3);
for i = 1:3
    customColormapIn(:,i) = linspace(startColor(i), endColorIn(i), n);
end

edges = -7e-5:1e-6:5e-5;

%% ME
figure(1);set(gcf,'Color','w')
numCols = size(FBP_Synapse_ME_R, 2);
numNeurons = size(FBP_Synapse_ME_R, 1);
RGBimage = zeros(numNeurons, 3 * numCols, 3);

for i = 1:numCols
    idx = 3*(i-1);
    dataL = FBP_Upstream_Synapse_ME_L(:, i);   % upstream input, contralateral (L)
    dataM = FBP_Synapse_ME_R(:, i);            % FBP output (R)
    dataR = FBP_Upstream_Synapse_ME_R(:, i);   % upstream input, ipsilateral (R)

    RGBimage(:,idx+1,:) = customColormapIn_ipsi( max(1, round(dataL*(n-1))) , :);
    RGBimage(:,idx+2,:) = customColormapOut(     max(1, round(dataM*(n-1))) , :);
    RGBimage(:,idx+3,:) = customColormapIn(      max(1, round(dataR*(n-1))) , :);
end

image(RGBimage);
axis tight;
xlabel('Input/Output columns (triplets: ipsi-out-contra)');
ylabel('Neuron index');
title('Disynaptic inputs and outputs (ME)');
set(gca, 'TickDir', 'out','YDir','normal','Box','off');
set(gca,'YTick',0.5:10:size(RGBimage,1)+0.5,'YTickLabel',edges(1:10:end));
set(gca,'XTick',2:3:size(RGBimage,2),'XTickLabel',FBP_Me.type);
hold on;
for i = 1:numCols-1
    xline(i*3 + 0.5, 'w--', 'LineWidth', 0.05);
end
ylim([0.5 90.5])
%% LO
figure(2);set(gcf,'Color','w')
numCols = size(FBP_Synapse_LO_R, 2);
numNeurons = size(FBP_Synapse_LO_R, 1);
RGBimage = zeros(numNeurons, 3 * numCols, 3);

for i = 1:numCols
    idx = 3*(i-1);
    dataL = FBP_Upstream_Synapse_LO_L(:, i);   % upstream input, contralateral (L)
    dataM = FBP_Synapse_LO_R(:, i);            % FBP output (R)
    dataR = FBP_Upstream_Synapse_LO_R(:, i);   % upstream input, ipsilateral (R)

    RGBimage(:,idx+1,:) = customColormapIn_ipsi( max(1, round(dataL*(n-1))) , :);
    RGBimage(:,idx+2,:) = customColormapOut(     max(1, round(dataM*(n-1))) , :);
    RGBimage(:,idx+3,:) = customColormapIn(      max(1, round(dataR*(n-1))) , :);
end

image(RGBimage);
axis tight;
xlabel('Input/Output columns (triplets: ipsi-out-contra)');
ylabel('Neuron index');
title('Disynaptic inputs and outputs (LO)');
set(gca, 'TickDir', 'out','YDir','normal','Box','off');
set(gca,'YTick',0.5:10:size(RGBimage,1)+0.5,'YTickLabel',edges(1:10:end));
set(gca,'XTick',2:3:size(RGBimage,2),'XTickLabel',FBP_Lo.type);
hold on;
for i = 1:numCols-1
    xline(i*3 + 0.5, 'w--', 'LineWidth', 0.05);
end
ylim([25.5 90.5])

%% LOP
figure(3);set(gcf,'Color','w')
numCols = size(FBP_Synapse_LOP_R, 2);
numNeurons = size(FBP_Synapse_LOP_R, 1);
RGBimage = zeros(numNeurons, 3 * numCols, 3);

for i = 1:numCols
    idx = 3*(i-1);
    dataL = FBP_Upstream_Synapse_LOP_L(:, i);  % upstream input, contralateral (L)
    dataM = FBP_Synapse_LOP_R(:, i);           % FBP output (R)
    dataR = FBP_Upstream_Synapse_LOP_R(:, i);  % upstream input, ipsilateral (R)

    RGBimage(:,idx+1,:) = customColormapIn_ipsi( max(1, round(dataL*(n-1))) , :);
    RGBimage(:,idx+2,:) = customColormapOut(     max(1, round(dataM*(n-1))) , :);
    RGBimage(:,idx+3,:) = customColormapIn(      max(1, round(dataR*(n-1))) , :);
end

image(RGBimage);
axis tight;
xlabel('Input/Output columns (triplets: ipsi-out-contra)');
ylabel('Neuron index');
title('Disynaptic inputs and outputs (LOP)');
set(gca, 'TickDir', 'out','Box','off');
set(gca,'YTick',0.5:10:size(RGBimage,1)+0.5,'YTickLabel',edges(1:10:end));
set(gca,'XTick',2:3:size(RGBimage,2),'XTickLabel',FBP_Lop.type);
hold on;
for i = 1:numCols-1
    xline(i*3 + 0.5, 'w--', 'LineWidth', 0.05);
end
ylim([60.5 100.5])
