%% Per-neuron input/output composition pie charts (panels S3C / S3D / S3E)
% For three example bidirectional neuron types — LC9 (S3C), LT43 (S3D), LT52 (S3E) —
% draws pie charts of their input and output partner composition. Each neuron yields
% 8 pies: input and output optic-vs-central overviews, then partner-type breakdowns
% (above a per-neuron synapse threshold) for total / central-only / optic-only inputs
% and outputs. Slices are colored by the partner type's most common superclass.
clear all; close all; clc;

baseDir = fileparts(fileparts(mfilename('fullpath')));  % resolved from this script's location (portable; no absolute paths)

opt = detectImportOptions(fullfile(baseDir, 'Codex_Data', 'consolidated_cell_types.csv'));
opt = setvartype(opt,'root_id','int64');
FAFB_consolidated_cell_types = readtable(fullfile(baseDir, 'Codex_Data', 'consolidated_cell_types.csv'), opt);

opt = detectImportOptions(fullfile(baseDir, 'Codex_Data', 'classification.csv'));
opt = setvartype(opt,'root_id','int64');
FAFB_classification = readtable(fullfile(baseDir, 'Codex_Data', 'classification.csv'), opt);

opt = detectImportOptions(fullfile(baseDir, 'Codex_Data', 'connections_no_threshold.csv'));
opt = setvartype(opt,'pre_root_id','int64');
opt = setvartype(opt,'post_root_id','int64');
FAFBConnections = readtable(fullfile(baseDir, 'Codex_Data', 'connections_no_threshold.csv'), opt);

% Build each example neuron's input/output partner tables (struct `Wantsee`) directly,
% via the local function seeConnection_by_region (no precomputed .mat needed).
WantToSee = {'LC9';'LT43';'LT52'};   % index order: 1=LC9 (S3C), 2=LT43 (S3D), 3=LT52 (S3E)
Wantsee = seeConnection_by_region(WantToSee, FAFBConnections, FAFB_consolidated_cell_types);

% Superclass color map
uniqueSuperclass = {'endocrine','motor','sensory','visual_projection', ...
    'visual_centrifugal','optic','central','descending','ascending'};
colors = [ 0.3010, 0.7450, 0.9330;   % endocrine
           0.2780, 0.6000, 0.8000;   % motor
           0.8500, 0.3250, 0.0980;   % sensory
           0.0000, 0.4470, 0.7410;   % visual_projection
           0.4660, 0.6740, 0.1880;   % visual_centrifugal
           0.9290, 0.6940, 0.1250;   % optic
           0.4940, 0.1840, 0.5560;   % central
           0.8500, 0.1500, 0.2000;   % descending
           0.6350, 0.5090, 0.2540 ]; % ascending

% Per-neuron config: Wantsee index, panel letter, and the 6 partner-type thresholds
% [Inputs Total, Outputs Total, Inputs Central, Inputs Optic, Outputs Central, Outputs Optic]
neurons = struct( ...
    'name',  {'LC9', 'LT43', 'LT52'}, ...
    'idx',   {1, 2, 3}, ...
    'panel', {'S3C', 'S3D', 'S3E'}, ...
    'thr',   {[1000 800 800 1100 700 500], ...
              [40 80 20 50 45 60], ...
              [700 500 160 700 300 400]});

for n = 1:numel(neurons)
    nm  = neurons(n).name;
    idx = neurons(n).idx;
    thr = neurons(n).thr;

    % --- Extract this neuron's partner tables ---
    Total_in   = cell2table(Wantsee.InNeuronTypes{idx}(:,[1 2]),  'VariableNames', {'Type','Synapse'});
    Total_out  = cell2table(Wantsee.OutNeuronTypes{idx}(:,[1 2]), 'VariableNames', {'Type','Synapse'});
    in_optic    = cell2table(Wantsee.unique_InOptic_Types{idx},   'VariableNames', {'Type','Synapse','Number'});
    in_central  = cell2table(Wantsee.unique_InCentral_Types{idx}, 'VariableNames', {'Type','Synapse','Number'});
    out_optic   = cell2table(Wantsee.unique_OutOptic_Types{idx},  'VariableNames', {'Type','Synapse','Number'});
    out_central = cell2table(Wantsee.unique_OutCentral_Types{idx},'VariableNames', {'Type','Synapse','Number'});

    % --- Optic vs central overviews ---
    draw_overview_pie(sum(in_central.Synapse),  sum(in_optic.Synapse),  sprintf('%s Input',  nm));
    draw_overview_pie(sum(out_central.Synapse), sum(out_optic.Synapse), sprintf('%s Output', nm));

    % --- Partner-type breakdowns ---
    draw_type_pie(Total_in,    thr(1), sprintf('%s Inputs Total',               nm), uniqueSuperclass, colors, FAFB_consolidated_cell_types, FAFB_classification);
    draw_type_pie(Total_out,   thr(2), sprintf('%s Outputs Total',              nm), uniqueSuperclass, colors, FAFB_consolidated_cell_types, FAFB_classification);
    draw_type_pie(in_central,  thr(3), sprintf('%s Inputs from Central brains', nm), uniqueSuperclass, colors, FAFB_consolidated_cell_types, FAFB_classification);
    draw_type_pie(in_optic,    thr(4), sprintf('%s Inputs from Optic lobes',    nm), uniqueSuperclass, colors, FAFB_consolidated_cell_types, FAFB_classification);
    draw_type_pie(out_central, thr(5), sprintf('%s Outputs to Central brains',  nm), uniqueSuperclass, colors, FAFB_consolidated_cell_types, FAFB_classification);
    draw_type_pie(out_optic,   thr(6), sprintf('%s Outputs to Optic lobes',     nm), uniqueSuperclass, colors, FAFB_consolidated_cell_types, FAFB_classification);
end

%% ===================== Local functions =====================
function draw_overview_pie(centralSyn, opticSyn, titleStr)
% Two-slice pie: total synapses in the central brain vs the optic lobes.
figure; set(gcf,'Color','w');
piechart([centralSyn, opticSyn], Names={'Central','Optic'});
title(titleStr);
end

function draw_type_pie(T, thr, titleStr, uniqueSuperclass, colors, FAFBConsolidated_type, FAFBClassification)
% Pie of partner types with synapse count >= thr, colored by each type's most common
% superclass and sorted by superclass.
T(T.Synapse < thr, :) = [];
if isempty(T), return; end

for i = 1:height(T)
    root_ids = FAFBConsolidated_type.root_id(strcmp(FAFBConsolidated_type.primary_type, T.Type{i}));
    sc   = FAFBClassification.super_class(ismember(FAFBClassification.root_id, root_ids));
    if isempty(sc)
        T.superclass{i} = '';   % unmatched type (no color)
    else
        [u,~,ix] = unique(sc);
        [~,mi]   = max(histcounts(ix, 1:numel(u)+1));
        T.superclass{i} = sc{mi};
    end
end

T = sortrows(T,'superclass','ascend');
T = sortrows(T,'Type','ascend');
T = sortrows(T,'superclass','ascend');

pieColors = zeros(height(T),3);
for i = 1:height(T)
    ci = matches(uniqueSuperclass, T.superclass{i});
    if any(ci), pieColors(i,:) = colors(ci,:); end
end

figure; set(gcf,'Color','w');
p = piechart(T, "Synapse", "Type", 'LabelStyle','name');
title(sprintf('%s (synapse > %d)', titleStr, thr));
p.ColorOrder = pieColors;
end

function [WantSee] = seeConnection_by_region(WantSeetypes, FAFBConnections, FAFB_consolidated_cell_types)
% For each requested cell type, collect its input/output partners (by partner type and
% by neuropil), split into optic-lobe vs central-brain, and return per-type/per-neuron
% summaries. Output struct columns include InNeuronTypes / OutNeuronTypes and
% unique_{In,Out}{Optic,Central}_Types ({type, total synapse, # roots}).

FAFBNeuropils=unique(FAFBConnections.neuropil);
FAFBNeuropil_OpticLobeRight={'AME_R','ME_R','LO_R','LOP_R','LA_R'};
FAFBNeuropil_OpticLobeLeft={'AME_L','ME_L','LO_L','LOP_L','LA_L'};
FAFBNeuropil_OpticLobe=[FAFBNeuropil_OpticLobeLeft, FAFBNeuropil_OpticLobeRight];
FAFBNeuropil_Central=FAFBNeuropils(~ismember(FAFBNeuropils,FAFBNeuropil_OpticLobeRight));
FAFBNeuropil_Central=FAFBNeuropil_Central(~ismember(FAFBNeuropil_Central,FAFBNeuropil_OpticLobeLeft));
FAFBNeuropil_Central=FAFBNeuropil_Central(~ismember(FAFBNeuropil_Central,'UNASGD'));
WantSee=table(WantSeetypes,'VariableNames',{'type'});

for i=1:1:size(WantSee,1)
    WantrootidsConsoltype=FAFB_consolidated_cell_types.root_id(strcmpi(FAFB_consolidated_cell_types.primary_type,WantSee.type{i}));
    Wantrootids=unique(WantrootidsConsoltype);
    WantSee.root_ids{i}=Wantrootids;
    if isempty(WantrootidsConsoltype)
        WantSee.root_ids{i}=sscanf(WantSee.type{i},'%ld');
    end

end

for i=1:1:size(WantSee,1)
    Want_root_ids=WantSee.root_ids{i};
    In_Connections=FAFBConnections(ismember(FAFBConnections.post_root_id,Want_root_ids),:);
    Out_Connections=FAFBConnections(ismember(FAFBConnections.pre_root_id,Want_root_ids),:);
    for j=1:1:size(In_Connections,1)
        In_root_ids=In_Connections.pre_root_id(j);
        idx_Consol=find(FAFB_consolidated_cell_types.root_id==In_root_ids);
        if ~isempty(idx_Consol)
            consolType=FAFB_consolidated_cell_types.primary_type{idx_Consol};

        else
            consolType='';

        end
        In_Connections.consoltype{j}=consolType;
        if ~isempty(consolType)
            In_Connections.type{j}=consolType;
        else
            In_Connections.type{j}=num2str(In_root_ids);
        end

    end

    for j=1:1:size(Out_Connections,1)
        Out_root_ids=Out_Connections.post_root_id(j);

        idx_Consol=find(FAFB_consolidated_cell_types.root_id==Out_root_ids);
        if ~isempty(idx_Consol)
            consolType=FAFB_consolidated_cell_types.primary_type{idx_Consol};

        else
            consolType='';

        end

        if ~isempty(consolType)
            Out_Connections.type{j}=consolType;
        else
            Out_Connections.type{j}=num2str(Out_root_ids);
        end
    end

    if isempty(In_Connections)
        WantSee.InNeuronTypes{i}={};
    else
        [uniqueInTypes,~,ic_in]=unique(In_Connections.type);

        for j=1:1:size(uniqueInTypes,1)
            uniqueInTypes{j,2}=sum(In_Connections.syn_count(ic_in==j));
            uniqueInTypes{j,3}=unique(In_Connections.pre_root_id(ic_in==j));
            temp=In_Connections(ic_in==j,:);
            temp=sortrows(temp,"syn_count","descend");
            temp=sortrows(temp,"post_root_id","ascend");
            temp=sortrows(temp,"pre_root_id","ascend");
            uniqueInTypes{j,4}=temp;
        end
        for j=1:1:size(uniqueInTypes,1)
            uniqueInTypes{j,5}=uniqueInTypes{j,2}/sum(cell2mat(uniqueInTypes(:,2)))*100;
        end
        WantSee.InNeuronTypes{i}=sortrows(uniqueInTypes,2,'descend');
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%% OUT
    if isempty(Out_Connections)
        WantSee.OutNeuronTypes{i}={};
    else
        [uniqueOutTypes,~,ic_out]=unique(Out_Connections.type);

        for j=1:1:size(uniqueOutTypes,1)
            uniqueOutTypes{j,2}=sum(Out_Connections.syn_count(ic_out==j));
            uniqueOutTypes{j,3}=unique(Out_Connections.post_root_id(ic_out==j));
            temp=Out_Connections(ic_out==j,:);
            temp=sortrows(temp,"syn_count","descend");
            temp=sortrows(temp,"post_root_id","ascend");
            temp=sortrows(temp,"pre_root_id","ascend");
            uniqueOutTypes{j,4}=temp;
        end
        for j=1:1:size(uniqueOutTypes,1)
            uniqueOutTypes{j,5}=uniqueOutTypes{j,2}/sum(cell2mat(uniqueOutTypes(:,2)))*100;
        end
        WantSee.OutNeuronTypes{i}=sortrows(uniqueOutTypes,2,'descend');
    end

    [unique_in_neuropil,~,ic_in]=unique(In_Connections.neuropil);

    for j=1:1:size(unique_in_neuropil,1)
        idx=ic_in==j;
        unique_in_neuropil{j,2}=sum(In_Connections.syn_count(idx));

    end
    WantSee.InNeuropils{i}=unique_in_neuropil;

    [unique_out_neuropil,~,ic_out]=unique(Out_Connections.neuropil);

    for j=1:1:size(unique_out_neuropil,1)
        idx=ic_out==j;
        unique_out_neuropil{j,2}=sum(Out_Connections.syn_count(idx));

    end
    WantSee.OutNeuropils{i}=unique_out_neuropil;

    if ~isempty(In_Connections)
        In_Connections = sortrows(In_Connections,"syn_count","descend");
        In_Connections = sortrows(In_Connections,"post_root_id","ascend");
        In_Connections = sortrows(In_Connections,"pre_root_id","ascend");
        In_Connections = sortrows(In_Connections,"type","ascend");

        WantSee.InConnections{i}=In_Connections;

    end

    if ~isempty(Out_Connections)
        Out_Connections = sortrows(Out_Connections,"syn_count","descend");
        Out_Connections = sortrows(Out_Connections,"post_root_id","ascend");
        Out_Connections = sortrows(Out_Connections,"pre_root_id","ascend");
        Out_Connections = sortrows(Out_Connections,"type","ascend");

        WantSee.OutConnections{i}=Out_Connections;

    end

end
for i = 1:1:size(WantSee,1)
    % ===== IN: connections into the current neuron =====
    % (NOTE) field name 'InConnections' follows the existing structure
    current_InConnections = WantSee.InConnections{i};

    InConnections_Optic   = current_InConnections(ismember(current_InConnections.neuropil, FAFBNeuropil_OpticLobe), :);
    InConnections_Central = current_InConnections(ismember(current_InConnections.neuropil, FAFBNeuropil_Central), :);

    % --- Optic IN: sum per pre_root_id + aggregate by type ---
    [unique_optic_in_neurons, ~, ic] = unique(InConnections_Optic.pre_root_id);
    unique_optic_in_neurons = [unique_optic_in_neurons, zeros(numel(unique_optic_in_neurons),1)]; % [root_id, sum_syn]
    Optic_in_neurons_Types = cell(0,3); % {type, total_syn, n_roots}

    for j = 1:1:size(unique_optic_in_neurons,1)
        idx = (ic == j);
        unique_optic_in_neurons(j,2) = sum(InConnections_Optic.syn_count(idx));

        % map to type
        rid = unique_optic_in_neurons(j,1);
        idx_type = (FAFB_consolidated_cell_types.root_id == rid);
        if any(idx_type)
            Optic_in_neurons_Types(j,1) = FAFB_consolidated_cell_types.primary_type(idx_type);
        else
            Optic_in_neurons_Types{j,1} = 'UnKnown';
        end
    end

    if ~isempty(Optic_in_neurons_Types)
        [unique_optic_in_neurons_Types, ~, icType] = unique(Optic_in_neurons_Types(:,1));
        unique_optic_in_neurons_Types(:,2) = {0}; % total syn
        unique_optic_in_neurons_Types(:,3) = {0}; % n roots

        for j = 1:1:size(unique_optic_in_neurons_Types,1)
            idx = (icType == j);
            % total syn of the roots of this type
            unique_optic_in_neurons_Types{j,2} = sum(unique_optic_in_neurons(idx,2));
            % number of unique roots of this type
            unique_optic_in_neurons_Types{j,3} = sum(idx);
        end

        unique_optic_in_neurons_Types = sortrows(unique_optic_in_neurons_Types, 2, 'descend');
    else
        unique_optic_in_neurons_Types = cell(0,3);
    end
    idx_unknown=strcmp(unique_optic_in_neurons_Types(:,1),'UnKnown');
    unique_optic_in_neurons_Types(idx_unknown,:)=[];
    WantSee.unique_InOptic_neurons{i} = unique_optic_in_neurons;
    WantSee.unique_InOptic_Types{i}   = unique_optic_in_neurons_Types;

    % --- Central IN: sum per pre_root_id + aggregate by type ---
    [unique_central_in_neurons, ~, ic] = unique(InConnections_Central.pre_root_id);
    unique_central_in_neurons = [unique_central_in_neurons, zeros(numel(unique_central_in_neurons),1)];
    Central_in_neurons_Types = cell(0,3);

    for j = 1:1:size(unique_central_in_neurons,1)
        idx = (ic == j);
        unique_central_in_neurons(j,2) = sum(InConnections_Central.syn_count(idx));

        rid = unique_central_in_neurons(j,1);
        idx_type = (FAFB_consolidated_cell_types.root_id == rid);
        if any(idx_type)
            Central_in_neurons_Types(j,1) = FAFB_consolidated_cell_types.primary_type(idx_type);
        else
            Central_in_neurons_Types{j,1} = 'UnKnown';
        end
    end

    if ~isempty(Central_in_neurons_Types)
        [unique_central_in_neurons_Types, ~, icType] = unique(Central_in_neurons_Types(:,1));
        unique_central_in_neurons_Types(:,2) = {0};
        unique_central_in_neurons_Types(:,3) = {0};

        for j = 1:1:size(unique_central_in_neurons_Types,1)
            idx = (icType == j);
            unique_central_in_neurons_Types{j,2} = sum(unique_central_in_neurons(idx,2));
            unique_central_in_neurons_Types{j,3} = sum(idx);
        end

        unique_central_in_neurons_Types = sortrows(unique_central_in_neurons_Types, 2, 'descend');
    else
        unique_central_in_neurons_Types = cell(0,3);
    end
    idx_unknown=strcmp(unique_central_in_neurons_Types(:,1),'UnKnown');
    unique_central_in_neurons_Types(idx_unknown,:)=[];

    WantSee.unique_InCentral_neurons{i} = unique_central_in_neurons;
    WantSee.unique_InCentral_Types{i}   = unique_central_in_neurons_Types;

    % ===== OUT: connections out of the current neuron =====
    current_OutConnections = WantSee.OutConnections{i};
    OutConnections_Optic   = current_OutConnections(ismember(current_OutConnections.neuropil, FAFBNeuropil_OpticLobe), :);
    OutConnections_Central = current_OutConnections(ismember(current_OutConnections.neuropil, FAFBNeuropil_Central), :);

    % --- Optic OUT: sum per post_root_id + aggregate by type ---
    [unique_optic_out_neurons, ~, ic] = unique(OutConnections_Optic.post_root_id); % post_root_id is correct here
    unique_optic_out_neurons = [unique_optic_out_neurons, zeros(numel(unique_optic_out_neurons),1)];
    Optic_out_neurons_Types = cell(0,3);

    for j = 1:1:size(unique_optic_out_neurons,1)
        idx = (ic == j);
        unique_optic_out_neurons(j,2) = sum(OutConnections_Optic.syn_count(idx));

        rid = unique_optic_out_neurons(j,1);
        idx_type = (FAFB_consolidated_cell_types.root_id == rid);
        if any(idx_type)
            Optic_out_neurons_Types(j,1) = FAFB_consolidated_cell_types.primary_type(idx_type);
        else
            Optic_out_neurons_Types{j,1} = 'UnKnown';
        end
    end

    if ~isempty(Optic_out_neurons_Types)
        [unique_optic_out_neurons_Types, ~, icType] = unique(Optic_out_neurons_Types(:,1));
        unique_optic_out_neurons_Types(:,2) = {0};
        unique_optic_out_neurons_Types(:,3) = {0};

        for j = 1:1:size(unique_optic_out_neurons_Types,1)
            idx = (icType == j);
            unique_optic_out_neurons_Types{j,2} = sum(unique_optic_out_neurons(idx,2));
            unique_optic_out_neurons_Types{j,3} = sum(idx);
        end

        unique_optic_out_neurons_Types = sortrows(unique_optic_out_neurons_Types, 2, 'descend');
    else
        unique_optic_out_neurons_Types = cell(0,3);
    end
    idx_unknown=strcmp(unique_optic_out_neurons_Types(:,1),'UnKnown');
    unique_optic_out_neurons_Types(idx_unknown,:)=[];

    WantSee.unique_OutOptic_neurons{i} = unique_optic_out_neurons;
    WantSee.unique_OutOptic_Types{i}   = unique_optic_out_neurons_Types;

    % --- Central OUT: sum per post_root_id + aggregate by type ---
    [unique_central_out_neurons, ~, ic] = unique(OutConnections_Central.post_root_id); % post_root_id
    unique_central_out_neurons = [unique_central_out_neurons, zeros(numel(unique_central_out_neurons),1)];
    Central_out_neurons_Types = cell(0,3);

    for j = 1:1:size(unique_central_out_neurons,1)
        idx = (ic == j);
        unique_central_out_neurons(j,2) = sum(OutConnections_Central.syn_count(idx));

        rid = unique_central_out_neurons(j,1);
        idx_type = (FAFB_consolidated_cell_types.root_id == rid);
        if any(idx_type)
            Central_out_neurons_Types(j,1) = FAFB_consolidated_cell_types.primary_type(idx_type);
        else
            Central_out_neurons_Types{j,1} = 'UnKnown';
        end
    end

    if ~isempty(Central_out_neurons_Types)
        [unique_central_out_neurons_Types, ~, icType] = unique(Central_out_neurons_Types(:,1));
        unique_central_out_neurons_Types(:,2) = {0};
        unique_central_out_neurons_Types(:,3) = {0};

        for j = 1:1:size(unique_central_out_neurons_Types,1)
            idx = (icType == j);
            unique_central_out_neurons_Types{j,2} = sum(unique_central_out_neurons(idx,2));
            unique_central_out_neurons_Types{j,3} = sum(idx);
        end

        unique_central_out_neurons_Types = sortrows(unique_central_out_neurons_Types, 2, 'descend');
    else
        unique_central_out_neurons_Types = cell(0,3);
    end
    idx_unknown=strcmp(unique_central_out_neurons_Types(:,1),'UnKnown');
    unique_central_out_neurons_Types(idx_unknown,:)=[];
    WantSee.unique_OutCentral_neurons{i} = unique_central_out_neurons;
    WantSee.unique_OutCentral_Types{i}   = unique_central_out_neurons_Types;
end
end
