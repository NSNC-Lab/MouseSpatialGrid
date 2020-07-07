%% Main script for calling mouse data simulation network
clearvars;
close all;

addpath('mechs')
addpath('dependencies')
addpath('eval_scripts')
addpath('genlib')
addpath(genpath('dynasim'))

makeGrids = 1;  % plot spatial grids
calcDistances = 1;  % plot VR distances
plot_rasters = 1;   % plot rasters

%% load expeimental data to optimize model to

dataCh = 25;

load('9-21-2016_0dB_removed_trialscleaned(-1,4).mat','Spks_clean','Spks_masked');
load('9-21-2016_0dB_removed_trials_performance.mat','Max','max_masked','opt_tau','opt_tau_masked');
data_perf = [Max{dataCh}(:);max_masked{dataCh}(:)];

% use tau of best clean location
[~,best_loc] = max(Max{dataCh});
data_tau = opt_tau{dataCh}(best_loc);

%% load STRF spikes

researchDrive = 'MiceSpatialGrids/';
ICdir = [researchDrive 'ICStim/Mouse/full_grids//BW_0.009 BTM_3.8 t0_0.1 phase0.4985//s30_STRFgain3.00_20200624-160717'];
ICdirPath = [ICdir filesep];
ICstruc = dir([ICdirPath '*.mat']);
if isempty(ICstruc), error('empty data directory'); end

%% varied parameters
varies(1).conxn = '(IC->IC)';
varies(1).param = 'trial';
varies(1).range = 1:40;

varies(end+1).conxn = 'R->C';
varies(end).param = 'gSYN';
varies(end).range = 0.03:0.03:0.21;

variedParam = 'R->CgSYN';

%% netcons
nCells = 4; %synchronise this variable with mouse_network

% -90, 0, 45, 90º
% x-channel inhibition
irNetcon = zeros(nCells);

srNetcon = zeros(nCells);

rcNetcon = [0;0;0;1]; %add this as input to mouse_network
% make rnNetcon have variable weights (instead of zeros)

netCons.irNetcon = irNetcon;
netCons.srNetcon = srNetcon;
netCons.rcNetcon = rcNetcon;

%% Initialize variables

nvaried = {varies(2:end).range};
nvaried = prod(cellfun(@length,nvaried));
datetime = datestr(now,'yyyymmdd-HHMMSS');

set(0, 'DefaultFigureVisible', 'off');
h = figure('Position',[50,50,850,690]);

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
    study_dir = fullfile(pwd, 'run', datetime, spatialConfig{1});
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
    
    [data(z).perf, data(z).fr, data(z).annot, ~, ~] = ...
        mouse_network(study_dir,time_end,varies,netCons,plot_rasters,calcDistances,...
        data_spks,data_tau);
    data(z).name = ICstruc(z).name;
end
save([pwd filesep 'run' filesep, datetime filesep 'summary_results.mat'],'data')
close(h);

%% performance grids
if makeGrids

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

width=13; hwratio=0.8;
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

figstr = 'SpatialGrid vary ';

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
    figstr = cat(2,figstr,variedParam,'%0.3f');
end

for vv = 1:nvaried
    
    h = figure('visible','on');
    figuresize(width, width*hwratio,h, 'inches')

    if ~isempty(mixedIdx)
        % mixed config cases
        for i = 1:length(mixedIdx)
            perf.IC(i,:) = data(mixedIdx(i)).perf.IC(:,vv);
            fr.IC(i,:) = data(mixedIdx(i)).fr.IC(:,vv);
            
            perf.R(i,:) = data(mixedIdx(i)).perf.R(:,vv);
            fr.R(i,:) = data(mixedIdx(i)).fr.R(:,vv);
            
            perf.C(i) = data(mixedIdx(i)).perf.C(vv);
            fr.C(i) = data(mixedIdx(i)).fr.C(vv);
        end

        % IC and relay neurons
        for nn = 1:numSpatialChan
            if nn > 1
                subplotloc = nn+1;
            else
                subplotloc = nn;
            end
            posVec = [1-(x0+subplotloc*(dx+lx)) y0+dy+ly lx ly];
            subplot('Position',posVec)
            plotPerfGrid(perf.IC(end:-1:1,nn),fr.IC(end:-1:1,nn),[],textColorThresh);
            
            if nn == 4
                ylabel('IC')
            else
                ylabel('');
            end
            xlabel(neurons{nn});
                            
            posVec = [1-(x0+subplotloc*(dx+lx)) y0+2*(dy+ly) lx ly];
            subplot('Position',posVec)
            plotPerfGrid(perf.R(end:-1:1,nn),fr.R(end:-1:1,nn),[],textColorThresh);
            if nn == 4
                ylabel('R')
            else
                ylabel('');
            end
            xlabel('');   
        end
        
        % C neuron

        posVec = [1-(x0+3*(dx+lx)) y0+3*(dy+ly) lx ly];
        subplot('Position',posVec);
        plotPerfGrid(perf.C(end:-1:1)',fr.C(end:-1:1)',[],textColorThresh);
        xlabel('')
    end
    
    subplot('position',[0.3 y0+dy/2 0.4 4*ly/5]); 
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
        end
    end    
    
    if ~exist('posVec','var')
       posVec = [x0+2*(lx+dx) y0+3*(dy+ly)]; 
    end
    posVec = [posVec(1) posVec(2)+ly+dy/4 lx ly/4];
    subplot('Position',posVec)
    plotPerfGrid(perf.CT(end:-1:1),[],[],textColorThresh);
    
    % simulation info
    annotation('textbox',[.75 .75 .2 .1],...
           'string',data(z).annot(vv,3:end),...
           'FitBoxToText','on',...
           'LineStyle','none')
       
    % show data grid next to model grid    
    posVec = [posVec(1)+dx+lx posVec(2) lx ly/4];
    subplot('Position',posVec)
    plotPerfGrid(data_perf(1:4)',[],[],textColorThresh);
    
    temp = flipud(reshape(data_perf(5:end),[4 4]));
    
    posVec = [posVec(1) posVec(2)-ly-dy/4 lx ly];
    subplot('Position',posVec)
    plotPerfGrid(temp,[],[],textColorThresh);
    xlabel('');
    ylabel('');
    
    % calculate error and correlation with data
    [cc_clean,MSE_clean] = calcModelPerf(perf.CT',data_perf(1:4));
    if exist('perf.C','var')
    [cc_masked,MSE_masked] = calcModelPerf(perf.C,data_perf(5:end));
    end

    str = {sprintf('Clean C.C. = %0.3f',cc_clean),sprintf('Clean MSE = %0.1f ± %0.1f',MSE_clean)};

    annotation('textbox',[.75 .7 .2 .1],...
           'string',str,...
           'FitBoxToText','on',...
           'LineStyle','none')
       
    % Show FR vs azimuth
    posVec = [x0 y0+3*ly+4*dy 2*lx+dx 5*ly/4];
    subplot('position',posVec)
    plot([-90 0 45 90],data_FR(targetIdx),[-90 0 45 90],fr.CT,'linewidth',2);
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
    if multiFlag == 1
        saveas(gca,[filesep DirPart filesep sprintf(figstr,paramPairs(vv,1),paramPairs(vv,2)) '.tiff'])
    else
        saveas(gca,[filesep DirPart filesep sprintf(figstr,varies(end).range(vv)) '.tiff'])
    end
    % clf

end

end


