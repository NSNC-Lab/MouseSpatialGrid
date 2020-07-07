clearvars;
close all;

addpath('mechs')
addpath('dependencies')
addpath('eval_scripts')
addpath('genlib')
addpath(genpath('dynasim'))

%% define parameters for gradient-descent

nIterations = 5;   % number of parameters searcches
nSims = 1; % number of simulations to run for average MSE
learningRate = 0.00001;
beta = 0.9; % 0 = without momentum

loss = zeros(1,nIterations);
convGoal = 0.1;

%param(1) - rcNetcon weights

param = [0.1;0;0;0.7]; %add this as input to mouse_network

%% load expeimental data to optimize model to

dataCh = 25;

load('9-21-2016_0dB_removed_trialscleaned(-1,4).mat','Spks_clean','Spks_masked');
load('9-21-2016_0dB_removed_trials_performance.mat','Max','max_masked');
data_perf = [Max{dataCh}(:);max_masked{dataCh}(:)];

%% load STRF spikes

researchDrive = 'MiceSpatialGrids/';
ICdir = [researchDrive 'ICStim/Mouse/full_grids//BW_0.009 BTM_3.8 t0_0.1 phase0.4985//s30_STRFgain3.00_20200624-160717'];
ICdirPath = [ICdir filesep];
ICstruc = dir([ICdirPath '*.mat']);
if isempty(ICstruc), error('empty data directory'); end

%% Initialize variables
plot_rasters = 0;
calcDistances = 0;
data_tau = 10;

datetime = datestr(now,'yyyymmdd-HHMMSS');

set(0, 'DefaultFigureVisible', 'off');
h = figure('Position',[50,50,850,690]);

del = 0;
for it = 1:nIterations

disp(['Starting iteration #', num2str(it)]);
    
%% varied parameters
varies(1).conxn = '(IC->IC)';
varies(1).param = 'trial';
varies(1).range = 1:40;

varies(2).conxn = 'R->C';
varies(2).param = 'gSYN';
varies(2).range = 0.21;

varies(3).conxn = 'C';
varies(3).param = 'noise';
varies(3).range = 1.6;

nvaried = {varies(2:end).range};
nvaried = prod(cellfun(@length,nvaried));

%% netcons
nCells = 4; %synchronise this variable with mouse_network

% % -90, 0, 45, 90º
% % x-channel inhibition
% irNetcon = zeros(nCells);
% 
% srNetcon = zeros(nCells);


% netCons.irNetcon = irNetcon;
% netCons.srNetcon = srNetcon;
rcNetcon = param;
netCons.rcNetcon = rcNetcon;

inputChans = find(rcNetcon ~= 0);

subz = find(contains({ICstruc.name},'m0.mat')); % sXm0 (target only) cases
%subz = 5:length({ICstruc.name});    % all cases except masker-only
for z = subz
    % restructure IC spikes
    load([ICdirPath ICstruc(z).name],'t_spiketimes');
    temp = cellfun(@max,t_spiketimes,'UniformOutput',false);
    tmax = max([temp{:}]);
    spks = zeros(20,4,tmax); %I'm storing spikes in a slightly different way...
    for j = 1:size(t_spiketimes,1) %trials [1:10]
        for k = 1:size(t_spiketimes,2) %neurons [(1:4),(1:4)]
            if k < 5 %song 1
                spks(j,k,round(t_spiketimes{j,k})) = 1;
            else
                spks(j+size(t_spiketimes,1),k-4,round(t_spiketimes{j,k})) = 1;
            end
        end
    end
    
    % save spk file
    spatialConfig = strsplit(ICstruc(z).name,'.');
    study_dir = fullfile(pwd, 'run', 'gradient-descent',datetime, spatialConfig{1});
    if exist(study_dir, 'dir')
      rmdir(study_dir, 's');
    end
    mkdir(fullfile(study_dir, 'solve'));
    save(fullfile(study_dir, 'solve','IC_spks.mat'),'spks');

    % call network
    h.Name = ICstruc(z).name;
    time_end = size(spks,3);
    
    % load spikes from data
    if ~contains(spatialConfig{1},'s0')
        if strcmp(spatialConfig{1}(4),'0')
            data_spks = squeeze(Spks_clean{dataCh}(:,5-str2double(spatialConfig{1}(2)),:));
        else
            data_spks = squeeze(Spks_masked{dataCh}(:,5-str2double(spatialConfig{1}(2)),5-str2double(spatialConfig{1}(4)),:));
        end
        data_FR(z) = mean(mean(cellfun(@(x) sum(x >= 0 & x < time_end/1000),data_spks)))/time_end*1000;
    else
        data_spks = [];
    end
    
    for ns = 1:nSims
    [temp_perf(ns), temp_fr(ns)] = ...
        mouse_network(study_dir,time_end,varies,netCons,plot_rasters,...
        calcDistances,data_spks,data_tau);
    end
    data(z).perf.R = mean([temp_perf.R],2);
    data(z).perf.C = mean([temp_perf.C]);
    
    data(z).fr.R = mean([temp_fr.R],2);
    data(z).fr.C = mean([temp_fr.C]);
    
    data(z).name = ICstruc(z).name;
end
save([pwd filesep 'run' filesep 'gradient-descent' filesep...
    datetime filesep 'summary_results_iteration' num2str(it) '.mat'],'data')
close(h);

%% performance grids

set(0,'defaultfigurevisible','on');
    
% performance vector has dimensions [numSpatialChan,nvaried]
neurons = {'Ipsi. sigmoid','Gaussian','U','Cont. sigmoid'};

temp = {data.name};
temp(cellfun('isempty',temp)) = {'empty'}; %label empty content

targetIdx = find(contains(temp,'m0') & ~strcmp(temp,'s0m0.mat'));
maskerIdx = find(contains(temp,'s0') & ~strcmp(temp,'s0m0.mat'));

%mixedIdx = setdiff(1:length(temp),[targetIdx,maskerIdx]);
mixedIdx = find(~contains(temp,'m0') & ~contains(temp,'s0') & ~contains(temp,'empty'));
textColorThresh = 70;
numSpatialChan = 4;

figstr = 'CleanGrid vary ';

% if more than one parameter is varied
if sum(cellfun(@length,{varies(2:end).range}) > 1) > 1
    
    multiFlag = 1;
    
    temp = find((cellfun(@length,{varies.range}) > 1) == 1);
    temp(1) = [];
    variedConxns = {varies(temp).conxn};
    variedParams = {varies(temp).param};
    variedRanges = {varies(temp).range};
    
    for i = 1:length(temp)
        figstr = cat(2,figstr,[variedConxns{i},variedParams{i},'%0.3f, ']);
    end
    figstr(end-1:end) = []; % delete extra ', '
    
    [A,B] = meshgrid(variedRanges{1},variedRanges{2});
    c = cat(2,A',B');
    paramPairs = reshape(c,[],2);
else
    multiFlag = 0;
    %figstr = cat(2,figstr,variedParam,'%0.3f');
end

width=8; hwratio=1;
x0=.08; y0=.08;
dx=.04; dy=.04;
lx=.125; ly=.125/hwratio;

x=-108:108;
tuningcurve=zeros(4,length(x));
ono = load('ono_curves_V2.mat','sigmoid','gauss','ushaped');
tuningcurve(1,:) = ono.sigmoid * rcNetcon(1);
tuningcurve(2,:) = ono.gauss * rcNetcon(2);
tuningcurve(3,:) = ono.ushaped * rcNetcon(3);
tuningcurve(4,:) = fliplr(ono.sigmoid) * rcNetcon(4);

for vv = 1:nvaried
    
    h = figure('visible','on');
    figuresize(width, width*hwratio,h, 'inches')
    
    subplot(2,2,1); 
    plot(x,tuningcurve','b','linewidth',1);
    hold on;
    plot(x,sum(tuningcurve',2),'k','linewidth',2);
    ylim([0 max(sum(tuningcurve',2))+0.2]);
    xlim([x(1) x(end)]);
    set(gca,'xdir','reverse');
    xticks([-90,0:45:90]);
    xlabel('Azimuth');
    
    % C neuron; target or masker only cases
    perf.CT = zeros(1,4);
    fr.CT = zeros(1,4);
    if ~isempty(targetIdx)
        for i = 1:length(targetIdx)
            perf.CT(i) = data(targetIdx(i)).perf.C(vv);
            fr.CT(i) = data(targetIdx(i)).fr.C(vv);
            fr.R(:,i) = data(targetIdx(i)).fr.R(:,vv);
        end
    end    
    
    % flip CT 
    perf.CT = fliplr(perf.CT);
    fr.R = fliplr(fr.R);
    
    subplot('Position',[0.6 0.6 0.2 0.2/4])
    plotPerfGrid(perf.CT,[],[],textColorThresh);
           
    % show data grid next to model grid    
    subplot('Position',[0.6 0.5 0.2 0.2/4])
    plotPerfGrid(data_perf(1:4)',[],[],textColorThresh);

    % calculate error and correlation with data
    [cc_clean,MSE_clean] = calcModelPerf(perf.CT',data_perf(1:4));
    % [cc_masked,MSE_masked] = calcModelPerf(perf.C,data_perf(5:end));

    str = {sprintf('Clean C.C. = %0.3f',cc_clean),...
        sprintf('Clean MSE = %0.1f ± %0.1f',MSE_clean),...
        ['learningRate = ' num2str(learningRate)],...
        ['beta = ' num2str(beta)]};
    for rr = 1:length(inputChans)
       str{end+1} = ['param(' num2str(inputChans(rr)) ') = ' num2str(param(inputChans(rr)))]; 
    end

    annotation('textbox',[0.6 .35 0.2 0.1],...
           'string',str,...
           'FitBoxToText','on',...
           'LineStyle','none')
       
    % Show FR vs azimuth
    subplot(2,2,3)
    plot([-90 0 45 90],data_FR(targetIdx),'-b',[-90 0 45 90],fr.CT,'-r','linewidth',2);
    hold on
    plot([-90 0 45 90],ones(1,4)*mean(data_FR(targetIdx)),'--b',...
        [-90 0 45 90],ones(1,4)*mean(fr.CT),'--r','linewidth',2);
    legend('Data','Model');
    xlabel('Azimuth');
    ylabel('Clean FR (Hz)')
    set(gca,'xdir','reverse');
    ylim([0 70]);
    xticks([-90,0:45:90]);
    % title(['RC weights = ' mat2str(round(100*rcNetcon)/100)])

    % save grid
    Dirparts = strsplit(study_dir, filesep);
    DirPart = fullfile(Dirparts{1:end-1});
    % if multiFlag == 1
    %    saveas(gca,[filesep DirPart filesep sprintf(figstr,paramPairs(vv,1),paramPairs(vv,2)) '.tiff'])
    % else
        saveas(gca,[filesep DirPart filesep 'iteration' num2str(it) '.tiff'])
    % end
    % clf

end

% update params and loss function

paramHistory(it,:) = param;
perfHistory(it,:) = perf.CT;
if it > 1 && abs(MSE_clean(1) - loss(it-1)) < convGoal
    break;
else
    loss(it) = MSE_clean(1);
    for rr = 1:length(inputChans)
        del = beta*del + learningRate*2*(perf.CT-data_perf(1:4)')*fr.R(rr,:)';
        param(inputChans(rr)) = param(inputChans(rr)) - del;
        if param(inputChans(rr)) < 0
           param(inputChans(rr)) = 0; 
        end
    end
end

end

figure; subplot(2,1,1); plot(loss);
ylabel('MSE');
subplot(2,1,2);
plot(paramHistory(:,inputChans)); xlabel('iteration');
legend('rcNetcon contra')
ylim([0 1]);

saveas(gcf,[filesep DirPart filesep 'MSE.png'])
save([filesep DirPart filesep 'optimization_results.mat'],'loss','paramHistory','perfHistory');
