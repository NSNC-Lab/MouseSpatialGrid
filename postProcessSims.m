
% trialStartTimes and trialEndTimes need to be the same cumsum as the trial
% lengths
trialStartTimes = zeros(1,length(subz)); %ms
for i = 1:length(subz) 
    trialStartTimes(i) = padToTime; %3500 ms
end

% start slicing data by trial lengths
options.trial_length = padToTime; % [ms]
PPtrialStartTimes = [1 cumsum(trialStartTimes)/dt+1]; % in [samples]
PPtrialEndTimes = PPtrialStartTimes(2:end)-(padToTime/dt-padToTime/dt+1);

% ICfiles contains names for spatial grid configs
load('ICfiles.mat','ICfiles');
configName = cellfun(@(x) x(1:end-4),{ICfiles.name}','UniformOutput',false);

% find varied parameter, excluding trials (1st entry in varies)
% if varied_param is empty, settle on 2nd entry in varies
varied_param = find( (cellfun(@length,{varies(2:end).range}) > 1 & ~cellfun(@iscolumn,{varies(2:end).range}))) + 1;
if isempty(varied_param), varied_param = 2; end
expVar = [varies(varied_param(1)).conxn '_' varies(varied_param(1)).param];
options.variedField = replace(expVar,'->','_');

annotTable = createSimNotes(snn_out,simDataDir,options);

% save output spikes and varied params to struct
names = snn_out(1).varied; results = struct;

if any(strcmp({s.populations.name},'C'))
    output_spks = 'C_V_spikes';
    output_pop = 'C';
else
    output_spks = 'R2On_V_spikes';
    output_pop = 'R2On';
end

for i = 1:length(snn_out)
    results(i).(output_spks) = snn_out(i).(output_spks);
    for t = 1:length(names)
        results(i).(names{t}) = snn_out(i).(names{t});
    end
end
results(1).model = snn_out(1).model; save([simDataDir filesep output_pop '_results.mat'],'results');

%% calculate discriminability and firing rate at each configuration
data = struct();

% calculate number of parameter sets (excluding repeats for optogenetic
% trials)
nVaried = length(snn_out)/(20*nSims);

tic;
for i = 1:length(subz)
    trialStart = PPtrialStartTimes(i); trialEnd = PPtrialEndTimes(i);

    figName = [simDataDir filesep configName{subz(i)}];
    [data(subz(i)).perf , data(subz(i)).fr , data(subz(i)).spks] = ...
        postProcessData_new(snn_out,s,trialStart,trialEnd,figName,options);
    data(subz(i)).config = configName{subz(i)};

    % tree-plotting functions: makes figures for all units for each config

    %plotRasterTree(data(subz(i)),figName,options); %close;
    % plotPSTHTree(data(subz(i)),figName,options); %close; 

    % make PSTH from spks

    t_bin = 20; % in [ms]
    psth_vec = (300:t_bin:(300 + 3000))/dt;
    for vv = 1:nVaried
        for tid = 1:2
            raster = data(subz(i)).spks.(output_pop)(vv).channel1((1:10) + 10*(tid-1),:);
            [~,spk_inds] = find(raster);
            data(subz(i)).output_PSTH(tid,:,vv) = histcounts(spk_inds,psth_vec);
        end
    end

    if nVaried == 1
        data(subz(i)).output_PSTH = squeeze(data(subz(i)).output_PSTH);
    end
end
toc;

% calculate control and laser performance for optogenetic trials
if nSims >= 2
    calcMeanOptoPerf(results,nVaried,simDataDir);
end
%%
pc = struct;
for i = 1:length(subz)
    [perf,FR] = calcPerfsandFR( data(subz(i)).spks , nVaried ,options.dt , output_pop);
    pc(subz(i)).fr = FR;
    pc(subz(i)).perf = perf;
    pc(subz(i)).config = configName{subz(i)};
end
save([simDataDir filesep 'perf_fr_' output_pop '.mat'],'pc');

% make surface plot for 2D parameter searches
if nVaried >= 10
    for i = 1:length(subz)
        plotPerfvsParams(output_pop,data(subz(i)),varies,simDataDir,data(subz(i)).config)
        %close all;
    end
end

dist_measures = struct;

for i = 1:length(subz)
   % calculate trial similarity and RMS difference at output for each
    % config
    clearvars outputSpks_cell
    for vv = 1:nVaried
        for ns = 1:nSims
            % if the data(indexing) doesn't work replace it with the snn_out
            % version: {snn_out((vv + (nS-1)*nVaried) : nVaried*nSims : end).C_V_spikes(trialStart:trialEnd)}
            outputSpks = data(subz(i)).spks.(output_pop)(vv).channel1;
            for n = 1:20
                outputSpks_cell{n} = find(outputSpks(n,:))' * dt/1000 - 0.3;
            end
            outputSpks_cell = reshape(outputSpks_cell,10,2);
            [dist_measures(subz(i)).TS(vv,ns),dist_measures(subz(i)).RMS(vv,ns)] = calcTrialSim(outputSpks_cell);
        end
    end
    dist_measures(subz(i)).config = configName{subz(i)};
end

save([simDataDir filesep 'TS_RMS_' output_pop '.mat'],'dist_measures')

%% Make full grid if we ran all spatial grid configurations

% if the locNum field is empty, the simulation runs all configurations
% (including masker-only) by default

% the second condition is if we run all configurations that have a target
% playing (will need to implement cases where locNum is an array instead of
% a single value)
if numel(subz) > 1
    simOptions.subz = subz;
    simOptions.locationLabels = strtrim(cellstr(num2str(locs'))');
    simOptions.chanLabels = chanLabels;

    subPops = {'C','ROn'};
    targetIdx = 5:5:20;
    mixedIdx = setdiff(1:24,[1:4 targetIdx]);

    plotPerformanceGrids_v3(data,s,annotTable,subPops,targetIdx,mixedIdx,simOptions,expName)
end

disp('here')

%% peakDrv for spatial grid stimuli

% ignore this for making full spatial grids

% getPeakDrv_SpatialStim;
