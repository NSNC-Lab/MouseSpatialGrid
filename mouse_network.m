function [perf, fr, annotstr, distMat, VR] = mouse_network(study_dir,time_end,...
    varies,plot_rasters,plot_distances,data_spks,data_tau,restrict_vary_flag)
% [performance, tauMax] = mouse_network(study_dir,time_end,varies,plot_rasters,data_spks)
% study_dir: location of IC spike files + directory for log and data files
% time_end: length of simulation in ms
% varies: vary parameter, a structure. e.g.
%   varies(1).conxn = '(IC->IC)';
%   varies(1).param = 'trial';
%   varies(1).range = 1:20;
%   *first set of parameters should always be "trial" for IC->IC cxn
% plot_rasters: 1 or 0
% data_spks: spikes from data we want to match 
%
% @Kenny F Chou, Boston Univ. 2019-06-20
% 2019-08-04 - added sharpening neurons
% 2019-08-14 - plotting now handles multiple varied parameters
% 2019-09-11 - removed redundant code. Return R performance in addition to
%              C performance
% 2020-02-05 - added netCons as parameter
% to do - fix plot_rasters option

% @Jio Nocon, Boston Univ., 2020-6-18
% 2020-6-18 - added subfunction to plot model rasters vs. data rasters and
%             calculate VR distance between the two

%% Input check
if ~strcmp(varies(1).param,'trial')
    error('first set of varied params should be ''trial''')
end

%% solver params
solverType = 'euler';
dt = 1; %ms % the IC input is currently dt=1

%% neuron populations
%tonic = bias = cells spontaneous firing

nCells = 4;
noise = 0.01; % low noise
s = struct();

s.populations(1).name = 'IC';
s.populations(1).equations = 'chouLIF';
s.populations(1).size = nCells;
s.populations(1).parameters = {'Itonic',0,'noise',0}; % 10-20 Hz spiking at rest

% s.populations(end+1).name = 'X';
% s.populations(end).equations = 'chouLIF';
% s.populations(end).size = nCells;
% s.populations(end).parameters = {'Itonic',0, 'noise',0}; % 10-20 Hz spiking at rest

% s.populations(end+1).name = 'S'; % S for sharpening
% s.populations(end).equations = 'chouLIF';
% s.populations(end).size = nCells;
% s.populations(end).parameters = {'Itonic',0, 'noise',0}; % 10-20 Hz spiking at rest

s.populations(end+1).name='R';
s.populations(end).equations = 'chouLIF';
s.populations(end).size = nCells;
s.populations(end).parameters = {'Itonic',0,'noise',0};

s.populations(end+1).name='C';
s.populations(end).equations = 'chouLIF';
s.populations(end).size = 1;
s.populations(end).parameters = {'noise',noise,'G_inc',0};

%% connections

s.connections(1).direction='IC->IC';
s.connections(1).mechanism_list={'IC_V2'};
s.connections(1).parameters={'g_postIC',0.25,'trial',5}; % 100 hz spiking

% s.connections(end+1).direction='IC->X';
% s.connections(end).mechanism_list={'synDoubleExp'};
% s.connections(end).parameters={'gSYN',0.21, 'tauR',0.4, 'tauD',2, 'netcon',diag(ones(1,nCells))}; 

% s.connections(end+1).direction='IC->S';
% s.connections(end).mechanism_list={'synDoubleExp'};
% s.connections(end).parameters={'gSYN',0.21, 'tauR',0.4, 'tauD',2, 'netcon',diag(ones(1,nCells))}; 

s.connections(end+1).direction='IC->R';
s.connections(end).mechanism_list={'synDoubleExp'};
s.connections(end).parameters={'gSYN',0.21, 'tauR',0.3, 'tauD',1.5, 'netcon', diag(ones(1,nCells)),'delay',0}; 

% s.connections(end+1).direction='X->R';
% s.connections(end).mechanism_list={'synDoubleExp'};
% s.connections(end).parameters={'gSYN',0.21, 'tauR',0.4, 'tauD',2, 'netcon',xrNetcon, 'ESYN',-80}; 
% 
% s.connections(end+1).direction='S->R';
% s.connections(end).mechanism_list={'synDoubleExp'};
% s.connections(end).parameters={'gSYN',0.21, 'tauR',0.4, 'tauD',5, 'netcon',srNetcon, 'ESYN',-80, 'delay',3}; 

s.connections(end+1).direction='R->C';
s.connections(end).mechanism_list={'synDoubleExp_V2'};
s.connections(end).parameters={'inputChan1',1,'inputChan2',1,'inputChan3',1,'inputChan4',1}; 

% if viz_network, vizNetwork; end

%% vary params
vary = cell(length(varies),3);
for i = 1:length(varies)
    vary{i,1} = varies(i).conxn;
    vary{i,2} = varies(i).param;
    vary{i,3} = varies(i).range;
end

numTrials = length(vary{1,3});

%% simulate
tic;
data = dsSimulate(s,'time_limits',[dt time_end], 'solver',solverType, 'dt',dt,...
  'downsample_factor',1, 'save_data_flag',0, 'save_results_flag',1,...
  'study_dir',study_dir, 'vary',vary, 'debug_flag',0, 'verbose_flag',0,...
  'restrict_vary_flag',restrict_vary_flag);


%% insert spikes
V_spike = 50;
for iData = 1:length(data)
    for pop = {s.populations.name}
        pop = pop{1};
        data(iData).([pop '_V'])(data(iData).([pop '_V_spikes']) == 1) = V_spike; % insert spike
    end
end

%%
jump = length(find([data.IC_IC_trial]==1));
for vv = 1:jump % for each varied parameter
    subData = data(vv:jump:length(data));

    %% visualize spikes
    ICspks = zeros(numTrials,4,time_end);
    %Sspks = zeros(numTrials,4,time_end);
    Rspks = zeros(numTrials,4,time_end);
    Cspks = zeros(numTrials,time_end);
    for i = 1:numTrials
        for j = 1:4
            ICspks(i,j,:) = subData(i).IC_V_spikes(:,j);
            %Sspks(i,j,:) = subData(i).S_V_spikes(:,j);
            Rspks(i,j,:) = subData(i).R_V_spikes(:,j);
        end
        Cspks(i,:) = subData(i).C_V_spikes;
    end

    % plot
    clf;
    locs = {'Cont. sigmoid','U-shaped','Gaussian','Ipsi. sigmoid'};
    for i = 1:4 %for each spatially directed neuron
        if i < 4
            ip = i;
        else
            ip = i + 1;
        end
        if plot_rasters, subplot(4,5,15+ip); end
        thisRaster = squeeze(ICspks(:,i,:));
        [perf.IC(i,vv),fr.IC(i,vv)] = calcPCandPlot(thisRaster,time_end,1,plot_rasters,numTrials);        

        if i==4, ylabel('IC'); end
        ax = get(gca,'position'); 
        annotation('textbox',[ax(1)+ax(3)/2 ax(2)-0.02 0 0],...
            'string',locs{i},...
            'HorizontalAlignment','center',...
            'LineStyle','none')
        
%         if plot_rasters, subplot(4,5,10+ip); end
% 
%         thisRaster = squeeze(Sspks(:,i,:));
%         calcPCandPlot(thisRaster,time_end,0,plot_rasters,numTrials);        
%         if i==4, ylabel('S'); end
%         xticklabels([])

        if plot_rasters, subplot(4,5,5+ip); end
        thisRaster = squeeze(Rspks(:,i,:));
        [perf.R(i,vv),fr.R(i,vv)] = calcPCandPlot(thisRaster,time_end,1,plot_rasters,numTrials);

        if i==4, ylabel('R'); end
        xticklabels([])
    end
    
    if plot_rasters, subplot(4,5,3); ylabel('C spikes'); xticklabels([]); end
    [perf.C(vv),fr.C(vv)] = calcPCandPlot(Cspks,time_end,1,plot_rasters,numTrials);     
    
    if plot_rasters, subplot(4,5,2); xticklabels([]); ylabel('Data spikes');
        tempspks = zeros(numel(data_spks),time_end);
        for tt = 1:numel(data_spks)
            temp = round(data_spks{tt}*1000);
            temp(temp < 0 | temp >= time_end) = [];
            tempspks(tt,temp+1) = 1;
        end
    calcPCandPlot(tempspks,time_end,1,plot_rasters,size(tempspks,1));  
    end
        
    % figure annotations
    annotstr(vv,:) = createAnnotStr(data(vv),varies(1:end-1));
    
    parts = strsplit(study_dir, filesep);
    DirPart = fullfile(parts{1:end-1});
    
    if plot_rasters
        annotation('textbox',[.675 .85 .2 .1],...
            'string',annotstr(vv,:),...
            'FitBoxToText','on',...
            'LineStyle','none')
        
        saveas(gca,[filesep DirPart filesep parts{end} '_v2_' num2str(vv) '.tiff'])
    end
    
    % plot comparison rasters and calculate VR distance
    
    [distMat(vv),VR(vv),taufig] = plotVRDists(Cspks,data_spks,data_tau,time_end,plot_distances);
        
    if plot_distances
        figure(taufig);
        annotation('textbox',[.75 .85 .2 .1],...
            'string',annotstr(vv,:),...
            'FitBoxToText','on',...
            'LineStyle','none')
        
        saveas(gca,[filesep DirPart filesep parts{end} '_' num2str(vv) '_VRdistances.tiff']) 
    end
    close(taufig);
end

toc;

end

function [pc,fr] = calcPCandPlot(raster,time_end,calcPC,plot_rasters,numTrials)
    PCstr = '';
    if calcPC
        % spks to spiketimes in a cell array of 20x2
        tau = linspace(1,30,100);
        spkTime = cell(numTrials,1);
        for ii = 1:numTrials, spkTime{ii} = find(raster(ii,:)); end
        spkTime = reshape(spkTime,numTrials/2,2);
        % calculate distance matrix & performance
        distMat = calcvr(spkTime, tau);
        [performance, ~] = calcpc(distMat, numTrials/2, 2, 1,[], 'new');
        pc = mean(max(performance));
        PCstr = ['PC = ' num2str(pc)];
    end
    
    fr = mean(sum(raster,2))/time_end*1000;
    %plot
    if plot_rasters
    plotSpikeRasterFs(flipud(logical(raster)), 'PlotType','vertline');
    title({PCstr,['FR = ' num2str(fr)]});
    xlim([0 time_end])
    line([0,time_end],[numTrials/2 + 0.5,numTrials/2 + 0.5],'color',[0.3 0.3 0.3])
    end
    
%     t_vec = 1:5:size(raster,2);
%     temp1 = mean(raster(1:numTrials/2,:));
%     temp2 = mean(raster(numTrials/2 + 1:end,:));
%     for t = 1:length(t_vec)-1
%         t1(t) = sum(temp1(t_vec(t):t_vec(t+1)));
%         t2(t) = sum(temp2(t_vec(t):t_vec(t+1)));
%     end
%     t1(end+1) = sum(temp1(t_vec(end):end));
%     t2(end+1) = sum(temp2(t_vec(end):end));
end

function [distMat,VR,taufig] = plotVRDists(raster,data_spks,data_tau,time_end,plot_distances)

% VR distances for model
numModel = size(raster,1);

% target 1
for tid = 1:2
    spkTime = cell(numModel/2,1);
    for ii = 1:numModel/2, spkTime{ii} = find(raster(numModel/2*(tid-1)+ii,:))/1000; end
    spkTime = spkTime(~cellfun('isempty',spkTime));
    % calculate distance matrix & performance
    distMat.model{tid} = calcvr(spkTime, data_tau/1000);
    
    temp_model{tid} = distMat.model{tid}(distMat.model{tid} ~= 0);
    VR.model(tid) = mean(temp_model{tid},'all');
end

% target 2

% VR distances for data
numExp = numel(data_spks);

for tid = 1:2
    spkTime = cell(size(data_spks(:,tid)));
    for ii = 1:length(data_spks(:,tid))
        spkTime{ii} = data_spks{ii,tid}';
        spkTime{ii}(spkTime{ii} < 0 | spkTime{ii} > time_end/1000) = []; 
    end
    
    spkTime = spkTime(~cellfun('isempty',spkTime));
    % calculate distance matrix & performance
    distMat.exp{tid} = calcvr(spkTime, data_tau/1000);
    
    temp_exp{tid} = distMat.exp{tid}(distMat.exp{tid} ~= 0);
    VR.exp(tid) = mean(temp_exp{tid},'all');
end

taufig = figure;

if plot_distances
allVRs = [distMat.exp{1}(:);distMat.exp{2}(:);distMat.model{1}(:);distMat.model{2}(:)];

maxVR = max(allVRs,[],'all');

binsize = 10;
bins = 0:binsize:maxVR+10;

for tid = 1:2
    [N_model{tid}] = histcounts(temp_model{tid},bins);
    [N_exp{tid}] = histcounts(temp_exp{tid},bins);
end
bins = bins + binsize/2;

subplot(2,1,1);
plot(bins(1:end-1),N_model{1},bins(1:end-1),N_model{2});
title(['Model VR: [' num2str(VR.model(1)) ', ' num2str(VR.model(2)) '], #trials = 40']);
xlim([0 maxVR]);
ylim([0 max([N_model{1},N_exp{1},N_model{2},N_exp{2}])]);

subplot(2,1,2);
plot(bins(1:end-1),N_exp{1},bins(1:end-1),N_exp{2});
xlabel('VR distance');
title(['Experimental VR: [' num2str(VR.exp(1)) ',' num2str(VR.exp(2)) '], #trials = ' num2str(numExp)]);
xlim([0 maxVR]);
ylim([0 max([N_model{1},N_exp{1},N_model{2},N_exp{2}])]);
end

end

function annotstr = createAnnotStr(data,varies)

paramstr = {data(1).varied{2:end}};
RCs = []; rc = 1;
gSYNs = []; gs = 1;
i = 1;
for aa = 1:length(varies)-1
    if contains(['data.' paramstr{aa}],'R_C_inputChan')
        RCs = cat(2,RCs,eval(['data.' paramstr{aa}]));
        rc = rc + 1;
    elseif contains(['data.' paramstr{aa}],'R_C_gSYN')
        gSYNs = cat(2,gSYNs,eval(['data.' paramstr{aa}]));
        gs = gs + 1;
    else
        annotstr{:,i} = sprintf('%s = %.3f',paramstr{aa},...
            eval(['data.' paramstr{aa}]));
        i = i + 1;
    end
end
annotstr{:,end+1} = ['RC_{gSYN} = ' mat2str(gSYNs)];

end