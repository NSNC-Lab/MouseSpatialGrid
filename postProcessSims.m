nVaried = length(snn_out)/(20*nSims);

% load ICfiles struct just for the names of the configs
load('ICfiles.mat');

% trialStartTimes and trialEndTimes need to be the same cumsum as the trial
% lengths
trialStartTimes = zeros(1,length(subz)); %ms
for i = 1:length(subz) 
    trialStartTimes(i) = padToTime; %3500 ms
end

data = struct();

% start slicing data by trial lengths
options.time_end = padToTime; %ms
PPtrialStartTimes = [1 cumsum(trialStartTimes)/dt+1]; %units of samples
PPtrialEndTimes = PPtrialStartTimes(2:end)-(padToTime/dt-options.time_end/dt+1);

configName = cellfun(@(x) strsplit(x,'_'),{ICfiles.name}','UniformOutput',false);
configName = vertcat(configName{:}); configName = configName(:,1);
options.variedField = strrep(expVar,'-','_');

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

tic;
if ~isempty(options.locNum) && numel(options.locNum) == 1
    trialStart = 1; trialEnd = padToTime/dt;
    figName = [simDataDir filesep configName{options.locNum}(1:end-4)];
    [data.perf,data.fr] = postProcessData_new(snn_out,s,trialStart,trialEnd,figName,options);
    plotRasterTree(snn_out,s,trialStart,trialEnd,figName,options); close
    plotPSTHTree(snn_out,s,trialStart,trialEnd,figName,options); close;
else
    for i = 1:length(subz)

        trialStart = PPtrialStartTimes(i); trialEnd = PPtrialEndTimes(i);

        figName = [simDataDir filesep configName{subz(i)}(1:end-4)];
        [data(subz(i)).perf,data(subz(i)).fr] = postProcessData_new(snn_out,s,trialStart,trialEnd,figName,options);
        plotRasterTree(snn_out,s,trialStart,trialEnd,figName,options); close;
        plotPSTHTree(snn_out,s,trialStart,trialEnd,figName,options); close;
    end
end
toc;

% close all;

% for optogenetic trials
if nSims == 5
    [pc,fr]= plotParamvsPerf_1D(results,nVaried);

    % performance
    pc_trials = struct2cell(pc);

    ctrl_mean = cellfun(@(x) mean(x(:,1)),pc_trials);
    laser_mean = cellfun(@(x) mean(x(:,end)),pc_trials);

    ctrl_se = cellfun(@(x) std(x(:,1))/sqrt(numel(x(:,1))),pc_trials);
    laser_se = cellfun(@(x) std(x(:,end))/sqrt(numel(x(:,end))),pc_trials);

    figure('unit','inches','position',[5 5 3 3]);
    bar((1:4)-.2,ctrl_mean,0.4,'facecolor','none','linewidth',2); hold on;
    bar((1:4)+.2,laser_mean,0.4,'facecolor','k','linewidth',2);
    xlim([0.4 4.6]);

    errorbar((1:4)-.2,ctrl_mean,ctrl_se,'color','k','linestyle','none','linewidth',1); hold on;
    errorbar((1:4)+.2,laser_mean,laser_se,'color','k','linestyle','none','linewidth',1);

    p_vals = cellfun(@(x) ranksum(x(:,1),x(:,2)),pc_trials)
    groups = mat2cell([[1:4]'-0.2 [1:4]'+0.2],[ 1 1 1 1 ]);
    ylim([50 100]);

    sigstar(groups,p_vals);

    set(gca,'xticklabels',{'SPIKE','ISI','RI-SPIKE','Spike count'},'xtick',1:4,'fontsize',8);
    ytickformat('percentage');
    ylabel('Performance'); legend('Control','Laser');
    saveas(gcf,[simDataDir filesep 'opto_performance_results.fig']);

    % firing rate
    ctrl_mean = mean(fr(:,1));
    laser_mean = mean(fr(:,end));

    ctrl_se = std(fr(:,1))/sqrt(5);
    laser_se = std(fr(:,end))/sqrt(5);

    figure('unit','inches','position',[5 5 2 3]);
    bar(1-.2,ctrl_mean,0.4,'facecolor','none','linewidth',2); hold on;
    bar(1+.2,laser_mean,0.4,'facecolor','k','linewidth',2);
    xlim([0.4 1.6]);

    errorbar(1-.2,ctrl_mean,ctrl_se,'color','k','linestyle','none','linewidth',1); hold on;
    errorbar(1+.2,laser_mean,laser_se,'color','k','linestyle','none','linewidth',1);

    ylim([0 60]);
    set(gca,'xticklabels',{'Control','Laser'},'xtick',[],'fontsize',8);
    ylabel('Firing rate (Hz)'); legend('Control','Laser');
    saveas(gcf,[simDataDir filesep 'opto_FR_results.fig']);
end

if ~isempty(options.locNum)
    [pc,fr] = plotParamvsPerf_1D(results,nVaried)
    save([simDataDir filesep 'perf_fr_C.mat'],'pc','fr')
end

if nVaried >= 10
    plotPerfvsParams('C',data,varies,simDataDir)
    close all;
end

% trial similarity and RMS difference
clearvars spks TS RMS

for nS = 1:nSims
    for nV = 1:nVaried

        spks = {snn_out((nV + (nS-1)*nVaried) : nVaried*nSims : end).C_V_spikes};
        for n = 1:20
            spks{n} = find(spks{n})/10000 - 0.3;
        end
        spks = reshape(spks,10,2);
        [TS(nV,nS),RMS(nV,nS)] = calcTrialSim(spks);

    end
end

save([simDataDir filesep 'TS_RMS_C.mat'],'TS','RMS');

%% Make full grid if we ran all spatial grid configurations

% if the locNum field is empty, the simulation runs all configurations
% (including masker-only) by default

% the second condition is if we run all configurations that have a target
% playing (will need to implement cases where locNum is an array instead of
% a single value)
if isempty(options.locNum) | all(ismember(5:24,options.locNum))
    simOptions.subz = subz;
    simOptions.locationLabels = {'90','45','0','-90'};
    simOptions.chanLabels = {'contra','center','ipsi'};

    subPops = {'C'};
    targetIdx = [5:5:20];
    mixedIdx = setdiff(1:24,[1:4 targetIdx]);

    plotPerformanceGrids_v3(data,s,annotTable,subPops,targetIdx,mixedIdx,simOptions,expName)
end

%% peakDrv for spatial grid stimuli

% ignore this for making full spatial grids

% getPeakDrv_SpatialStim;
