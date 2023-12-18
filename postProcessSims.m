
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

% save C spikes and varied params to struct
names = snn_out(1).varied; results = struct;
for i = 1:length(snn_out)
    results(i).C_V_spikes = snn_out(i).C_V_spikes;
    for t = 1:length(names)
        results(i).(names{t}) = snn_out(i).(names{t});
    end
end
results(1).model = snn_out(1).model; save([simDataDir filesep 'C_results.mat'],'results');

%% calculate discriminability and firing rate at each configuration
data = struct();

tic;
for i = 1:length(subz)
    trialStart = PPtrialStartTimes(i); trialEnd = PPtrialEndTimes(i);

    figName = [simDataDir filesep configName{subz(i)}];
    [data(subz(i)).perf , data(subz(i)).fr , data(subz(i)).spks] = ...
        postProcessData_new(snn_out,s,trialStart,trialEnd,figName,options);
    data(subz(i)).config = configName{subz(i)};

    % tree-plotting functions: makes figures for all units for each config

    plotRasterTree(data(subz(i)),figName,options); close;
    plotPSTHTree(data(subz(i)),figName,options); close; 

    % make PSTH from spks

    t_bin = 20; % in [ms]
    psth_vec = (300:t_bin:(300 + 3000))/dt;
    for tid = 1:2
        raster = data(subz(i)).spks.C.channel1((1:10) + 10*(tid-1),:);
        [~,spk_inds] = find(raster);
        data(subz(i)).output_PSTH(tid,:) = histcounts(spk_inds,psth_vec);
    end    
end
toc;

% % calculate number of parameter sets (excluding repeats for optogenetic
% % trials)
% nVaried = length(snn_out)/(20*nSims);
% 
% % calculate control and laser performance for optogenetic trials
% if nSims == 5
%     calcMeanOptoPerf(results,nVaried,simDataDir);
% end
% 
% if ~isempty(options.locNum)
%     [pc,fr] = plotParamvsPerf_1D(results,nVaried,options.dt);
%     save([simDataDir filesep 'perf_fr_C.mat'],'pc','fr');
% end
% 
% % make surface plot for 2D parameter searches
% if nVaried >= 10
%     plotPerfvsParams('C',data,varies,simDataDir)
%     close all;
% end
% 
% % % calculate trial similarity and RMS difference
% % clearvars TS RMS
% % 
% % for nS = 1:nSims
% %     for nV = 1:nVaried
% %         outputSpks = {snn_out((nV + (nS-1)*nVaried) : nVaried*nSims : end).C_V_spikes};
% %         for n = 1:20
% %             outputSpks{n} = find(outputSpks{n})*dt/1000 - 0.3;
% %         end
% %         outputSpks = reshape(outputSpks,10,2);
% %         [TS(nV,nS),RMS(nV,nS)] = calcTrialSim(outputSpks);
% %     end
% % end
% % 
% % save([simDataDir filesep 'TS_RMS_C.mat'],'TS','RMS');
% 
% %% Make full grid if we ran all spatial grid configurations
% 
% % if the locNum field is empty, the simulation runs all configurations
% % (including masker-only) by default
% 
% % the second condition is if we run all configurations that have a target
% % playing (will need to implement cases where locNum is an array instead of
% % a single value)
% if numel(subz) > 1
%     simOptions.subz = subz;
%     simOptions.locationLabels = strtrim(cellstr(num2str(locs'))');
%     simOptions.chanLabels = chanLabels;
% 
%     subPops = {'C','ROn'};
%     targetIdx = 5:5:20;
%     mixedIdx = setdiff(1:24,[1:4 targetIdx]);
% 
%     plotPerformanceGrids_v3(data,s,annotTable,subPops,targetIdx,mixedIdx,simOptions,expName)
% end
% 
% %% peakDrv for spatial grid stimuli
% 
% % ignore this for making full spatial grids
% 
% % getPeakDrv_SpatialStim;
