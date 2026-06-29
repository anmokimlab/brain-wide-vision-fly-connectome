%% fig_4H_FBP_loop
%  Loop-forming FBP neuron percentage (paper panel 4H).
%
%  For every right-hemisphere FBP (feedback) neuron `a`, look at its
%  "disynaptic-visual-input" partners -- the central-brain (CB) input neurons
%  `b` (b -> a in CB) that themselves receive optic-lobe (OL) input -- and ask
%  whether `a` closes a feedback loop back onto any of them through its OL output:
%     direct   :  a -> b           (one OL output hop)        => a -> b -> a
%     indirect :  a -> x -> b      (two OL hops, x ~= b)      => a -> x -> b -> a
%  A focal neuron "forms a loop" (direct / indirect) if it reaches >= 1 such
%  partner. The panel reports, per FB neuropil subgroup (Me / Lo / Lop / Multi),
%  the PERCENTAGE of that subgroup's focal neurons that form a direct / indirect
%  loop (count-based; each focal neuron counted once).
%
%  Reads connections_no_threshold.csv and consolidated_cell_types.csv from
%  Codex_Data\ and the FBP neuron set (RightFBP_NPIs) from
%  Processed_Data\right_neurons_thr0.mat (produced by
%  Figures\fig_1D_E_postPI_prePI_right_neurons.m). FBP types are grouped into
%  Me / Lo / Lop / Multi by their target optic lobe (hardcoded index lists from
%  the FBP output-neuropil classification, s10).
clear all; close all; clc

baseDir  = fileparts(fileparts(mfilename('fullpath')));  % resolved from this script's location (portable; no absolute paths)
codexDir = fullfile(baseDir, 'Codex_Data');

THRESHOLD = 5;   % drop every edge with syn_count < THRESHOLD

%% --- Focal FBP neuron set ---
load(fullfile(baseDir, 'Processed_Data', 'right_neurons_thr0.mat'), 'RightFBP_NPIs')

%% --- Load connections / cell types ---
fprintf('Loading connections ...\n');
opt = detectImportOptions(fullfile(codexDir, 'connections_no_threshold.csv'));
opt.SelectedVariableNames = {'pre_root_id','post_root_id','neuropil','syn_count'};
opt = setvartype(opt,'pre_root_id','int64');
opt = setvartype(opt,'post_root_id','int64');
opt = setvartype(opt,'syn_count','double');
FAFBConnections = readtable(fullfile(codexDir, 'connections_no_threshold.csv'), opt);
FAFBConnections(FAFBConnections.syn_count < THRESHOLD, :) = [];

opt2 = detectImportOptions(fullfile(codexDir, 'consolidated_cell_types.csv'));
opt2 = setvartype(opt2,'root_id','int64');
FAFBConsolidatedTypes = readtable(fullfile(codexDir, 'consolidated_cell_types.csv'), opt2);

opticLobeNPL = {'ME_R','AME_R','LO_R','LOP_R','LA_R', ...
                'ME_L','AME_L','LO_L','LOP_L','LA_L'};
excludeNPL   = {'UNASGD'};

%% --- Neuron-level adjacency (binary, OL vs CB) ---
fprintf('Building neuron-level adjacency ...\n');
ids = unique([FAFBConnections.pre_root_id; FAFBConnections.post_root_id]);
N   = numel(ids);
[~, preI]  = ismember(FAFBConnections.pre_root_id,  ids);
[~, postI] = ismember(FAFBConnections.post_root_id, ids);
isOL = ismember(FAFBConnections.neuropil, opticLobeNPL);
isCB = ~isOL & ~ismember(FAFBConnections.neuropil, excludeNPL);   % central brain (drop UNASGD)

A_OL = spones(sparse(preI(isOL), postI(isOL), 1, N, N));   % optic lobe (i -> j exists)
A_CB = spones(sparse(preI(isCB), postI(isCB), 1, N, N));   % central brain

% neurons receiving >= 1 incoming OL synapse (restricts the input base to
% "disynaptic visual input" neurons)
receivesOL = full(sum(A_OL, 1).' > 0);   % N x 1 logical

% input base = CB, both output hops = OL
A_first = A_OL;  A_second = A_OL;  A_in = A_CB;

% neuron index -> primary_type
ptype = repmat("UNTYPED", N, 1);
[tfC, locC] = ismember(ids, FAFBConsolidatedTypes.root_id);
pt = string(FAFBConsolidatedTypes.primary_type);
pt(strlength(pt)==0 | ismissing(pt)) = "UNTYPED";
ptype(tfC) = pt(locC(tfC));

%% --- Group the FBP neurons by cell type ---
[RightFBP_type,~,ic] = unique(RightFBP_NPIs.type);
RightFBP_type = table(RightFBP_type,'VariableNames',{'type'});
for i = 1:size(RightFBP_type,1)
    RightFBP_type.root_id{i} = RightFBP_NPIs.root_id(ic==i);
end

% Group FBP types by their target optic lobe (Me / Lo / Lop / Multi).
% (index lists into RightFBP_type, from the FBP output-neuropil classification, s10)
Me_FBP_idx   = [4 22 25 26 29 30 31 32 33 68 69 70 71 72 73 74 75 76 78 81 82 83 84 85 86 87 93];
Lo_FBP_idx   = [2 10 11 12 13 14 15 16 17 18 19 21 34 35 36 37 38 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 95];
Lop_FBP_idx  = [3 7 8 9 28 62 63 64 65 66];
Multi_FBP_idx= [1 5 6 23 24 27 39 60 61 67 77 79 80 88 89 90 91 92 94 96];

FBP_Me    = RightFBP_type(Me_FBP_idx,:);
FBP_Lo    = RightFBP_type(Lo_FBP_idx,:);
FBP_Lop   = RightFBP_type(Lop_FBP_idx,:);
FBP_Multi = RightFBP_type(Multi_FBP_idx,:);
FBP_subgroups = {FBP_Me, FBP_Lo, FBP_Lop, FBP_Multi};
group_labels  = {'FBP\_Me','FBP\_Lo','FBP\_Lop','FBP\_Multi'};

%% --- Per-neuron loop participation (count of reached input partners) ---
npiRid      = int64(RightFBP_NPIs.root_id);
FOCAL_TYPES = RightFBP_type.type;
perNeuron   = table();
for it = 1:numel(FOCAL_TYPES)
    tname  = FOCAL_TYPES{it};
    isFoc  = (ptype == tname) & ismember(ids, npiRid);
    focIdx = find(isFoc);
    n = numel(focIdx);
    if n == 0, continue; end

    [nDirect, nVia1] = deal(zeros(n,1));
    for k = 1:n
        ia = focIdx(k);
        out1 = find(A_first(ia,:));  out1 = out1(:);  out1(out1==ia) = [];   % direct OL output targets
        out1Vec = false(N,1);  out1Vec(out1) = true;
        out2 = find(any(A_second(out1,:), 1));  out2(out2==ia) = [];         % a -> x -> b (OL -> OL)
        out2Vec = false(N,1);  out2Vec(out2) = true;
        via1Vec = out2Vec & ~out1Vec;                                        % via 1 intermediary, excl. direct

        inIdx = find(A_in(:,ia));  inIdx(inIdx==ia) = [];                    % CB inputs b -> a
        inIdx = inIdx(receivesOL(inIdx));                                    % keep only b that receive OL input
        if isempty(inIdx), continue; end

        nDirect(k) = nnz(out1Vec(inIdx));
        nVia1(k)   = nnz(via1Vec(inIdx));
    end
    perNeuron = [perNeuron; table(repmat({tname},n,1), nDirect, nVia1, ...
        'VariableNames', {'type','n_direct','n_via1'})]; %#ok<AGROW>
end

%% --- Figure 1 (panel 4H): loop-forming FBP neuron percentage per subgroup ---
%  Per FB neuropil subgroup, % of focal neurons that take part in >= 1 feedback
%  cycle (direct = a->b->a ; indirect = a->x->b->a). Each focal neuron counted once.
loopLbls = {'direct','indirect'};
[cnt, tot] = subgroupLoopCounts(perNeuron, FBP_subgroups);
bar_data = 100 * cnt ./ max(tot, 1);
glab = erase(string(group_labels), '\');

fprintf('=== Loop-forming FBP neurons (%% of focal neurons) ===\n');
for g = 1:numel(FBP_subgroups)
    fprintf('   %-9s : direct %5.1f%%  indirect %5.1f%%   (n=%d)\n', ...
        glab(g), bar_data(g,1), bar_data(g,2), tot(g));
end

figure(1); set(gcf, 'Color','w', 'Position',[100 100 700 500]);
b = bar(bar_data, 'grouped');
grayTones = linspace(0.75, 0.25, 2)';
for m = 1:2, b(m).FaceColor = grayTones(m)*[1 1 1]; end
set(gca,'XTick',1:numel(FBP_subgroups),'XTickLabel',cellstr(glab),'TickLabelInterpreter','none');
ylabel('loop-forming focal neurons (%)');
title('Loop-forming FBP neurons (\geq1 feedback cycle) per FB neuropil subgroup');
legend(loopLbls,'Location','best');  grid on;
set(gca,'Box','off','TickDir','out');


%% ===== Helper functions =====
function [cnt, tot] = subgroupLoopCounts(perNeuron, subgroups)
%  Pool ALL focal neurons of each subgroup's types and count each neuron ONCE per
%  class if it takes part in >= 1 cycle:
%     direct   (col 1) = n_direct > 0   (a -> b -> a)
%     indirect (col 2) = n_via1  > 0    (a -> x -> b -> a)
%  Returns cnt (nG x 2 looped-neuron counts) and tot (nG x 1 total focal neurons).
    nG    = numel(subgroups);
    types = string(perNeuron.type);
    cnt = zeros(nG, 2);
    tot = zeros(nG, 1);
    for g = 1:nG
        sel = ismember(types, string(subgroups{g}.type));
        tot(g)   = nnz(sel);
        cnt(g,1) = nnz(sel & perNeuron.n_direct > 0);
        cnt(g,2) = nnz(sel & perNeuron.n_via1  > 0);
    end
end
