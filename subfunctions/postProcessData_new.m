function [perf,fr,spks] = postProcessData_new(data,s,trialStart,trialEnd,configName,options)
% calculate performance and FR for *a single spot on the spatial grid*
% input:
%   data structure with length = #sims, containing voltage information of
%   [time x spatial channel]. Time dimension contains information from each
%   spot in the spatial grid. Extract this information with
%   options.trialStart and options.trialEnd.
% output:
%   data
%       spatial channels,trials,time
% no plotting function for now

trial_length = options.trial_length;
plot_rasters = options.plotRasters;
dt = options.dt;
fields = fieldnames(data);
ind = find(contains(fields,'_trial'),1);
jump = length(find([data.(fields{ind})]==1));  % number of variations/parameter sets
numTrials = length(data)/jump; %usually, 20 trials
% numChannels = 4;

% network properties
popNames = {s.populations.name};
popSizes = [data(1).model.specification.populations.size];
nPops = numel(popNames);
fieldNames = strcat(popNames,'_V_spikes');

% for this trial
tstart = trialStart;
tend = trialEnd;

[y1,fs_stim] = audioread('200k_target1.wav');
y2 = audioread('200k_target2.wav');
t = (0:(length(y1)-1))/fs_stim;

% visualize spikes for specified populations
if ~isfield(options,'subPops'), options.subPops = popNames; end
if plot_rasters, figure; end


for vv = 1:jump % for each varied parameter
    subData = data(vv:jump:length(data)); %grab data for this param variation
    try
        variedParamVal = mode([subData.(options.variedField)]); % should all be the same
    catch
        variedParamVal = 0;
    end
    
    for currentPop = 1:nPops
        
        if popSizes(currentPop) > 1 && isfield(options,'chansToPlot')
            chansToPlot = options.chansToPlot;
        else
            chansToPlot = 1:popSizes(currentPop);
        end

        % skip processing for current population if not within specified subpopulation.
        if ~contains(popNames(currentPop),options.subPops), continue; end
        
        % for each spatial channel
        for channelNum = 1:popSizes(currentPop)
            channel = struct();
            
            % for each trial
            for trial = 1:numTrials
                channel(channelNum).popSpks(trial,:) = subData(trial).(fieldNames{currentPop})(tstart:tend,channelNum);
            end
            spks.(popNames{currentPop})(vv).(['channel' num2str(channelNum)]) = channel(channelNum).popSpks;
            
            figName = sprintf('%s_%s%.03f_%s_channel%i_%s',configName,options.variedField,variedParamVal,popNames{currentPop},channelNum,num2str(vv));
            
            if (strcmp(popNames{currentPop},'On') || strcmp(popNames{currentPop},'Off'))  && vv > 1
                plot_rasters_final = 0;
            elseif ismember(channelNum,chansToPlot)
                plot_rasters_final = plot_rasters;
            else
                plot_rasters_final = 0;
            end
            
            [perf.(popNames{currentPop}).(['channel' num2str(channelNum)])(vv),...
                fr.(popNames{currentPop}).(['channel' num2str(channelNum)])(vv)] = ...
                calcPCandPlot(channel(channelNum).popSpks,trial_length,numTrials,dt,plot_rasters_final,y1,y2,t,figName);
            
        end
    end
end

end

function [pc,fr,psth] = calcPCandPlot(raster,trial_length,numTrials,dt,plot_rasters,y1,y2,t,figName)

% inputs:
% raster - 0s and 1s with size [trials x samples]
% trial_length - simulation time in [ms]
% calcPC - if performance is calculated == 1
% numTrials - # trials per target ( #rows in raster / 2 )
% dt - timestep [in ms]
% plot_rasters - if plotting config/unit, == 1
% y1 and y2 - target waveforms
% t - target time vector
% figName - figure name

PCstr = '';

% use dt to calculate indexes for stimulus response
start_time = 300; % in [ms]
end_time = start_time + 3000; % in [ms]

% spks to spiketimes in a cell array of 20x2
spkTimes = cell(numTrials,1);
for ii = 1:numTrials
    % convert raster spike indexes to ms
    spkTimes{ii} = find(raster(ii,:))*dt;
end
spkTimes = reshape(spkTimes,numTrials/2,2);
input = reshape(spkTimes,1,numTrials);
fr = round(mean(cellfun(@(x) sum(x >= start_time & x < end_time) / 3,input)));

% claculate performance
STS = SpikeTrainSet(input,start_time,end_time);
distMat = STS.SPIKEdistanceMatrix(start_time,end_time);

performance = calcpcStatic(distMat, numTrials/2, 2, 0);
pc = mean(max(performance));
PCstr = ['PC = ' num2str(round(pc)) '%'];

%plot
x = 0.86;
y_raster = 0.4;
y_stim = 0.08;
y_psth = 0.12;
dy = 0.01;
x0 = 0.1;

% plot rasters and PSTHs (in seconds)
if plot_rasters
    clf
    
    ypos = 0.91 - y_stim;
    subplot('position',[x0 ypos x y_stim]); % target 2
    plot(t,y2,'k'); xlim([0 trial_length/1000]);
    title({PCstr,['FR = ' num2str(fr)]}); set(gca,'xtick',[],'ytick',[])
    
    ypos = ypos - dy - y_psth;
    subplot('position',[x0 ypos x y_psth]);
    plotPSTH(raster(numTrials/2+1:end,:),dt); set(gca,'xtick',[])
    
    ypos = ypos - dy - y_raster;
    subplot('position',[x0 ypos x y_raster]);
    plotSpikeRasterFs(flipud(logical(raster)), 'PlotType','vertline');
    xlim([0 trial_length/dt]);
    line([0,trial_length/dt],[numTrials/2 + 0.5,numTrials/2 + 0.5],'color',[0.3 0.3 0.3]); set(gca,'xtick',[],'ytick',[])
    
    ypos = ypos - dy - y_psth;
    subplot('position',[x0 ypos x y_psth]);
    plotPSTH(raster(1:numTrials/2,:),dt); set(gca,'xtick',[])

    ypos = ypos - dy - y_stim;
    subplot('position',[x0 ypos x y_stim]); % target 1
    plot(t,y1,'k'); xlim([0 trial_length/1000]); set(gca,'ytick',[]); 
    
    saveas(gcf,[figName '.png'])
end
    
end

function plotPSTH(raster,dt)

% plots PSTH (in seconds)

t_vec = 0:0.02:3.5;
[~,inds] = find(raster);
spktimes = inds*dt/1000;

psth = histcounts(spktimes,t_vec);
psth(end+1) = 0;

bar(t_vec,psth,'k');

lim = 30;
if max(psth,[],'all') > lim
   lim = lim*ceil(max(psth,[],'all')/lim); 
end

xlim([t_vec(1) t_vec(end)]); ylim([0 lim]);
ylabel('Spike count');
set(gca,'ytick',0:lim/3:lim);

end