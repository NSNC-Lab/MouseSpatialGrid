clearvars;
close all;

addpath('mechs')
addpath('dependencies')
addpath('eval_scripts')
addpath('genlib')
addpath(genpath('dynasim'))

%%%%%%%% start of user inputs

makeGrids = 1;  % plot spatial grids
calcDistances = 1;  % plot VR distances
plot_rasters = 0;   % plot rasters

%% define parameters for optimization

nSims = 3; % number of simulations to run for average MSE

temp = 1;
weights = [zeros(size(temp,2),3),temp']; %add this as input to mouse_network

%%%%%%%% end of user inputs

inputChans = find(weights(1,:) ~= 0);

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
ICdir = [researchDrive 'ICStim/Mouse/full_grids//BW_0.009 BTM_3.8 t0_0.1 phase0.4985//s30_STRFgain4.50_20200622-150507'];
ICdirPath = [ICdir filesep];
ICstruc = dir([ICdirPath '*.mat']);
if isempty(ICstruc), error('empty data directory'); end

%% Initialize variables

datetime = datestr(now,'yyyymmdd-HHMMSS');

set(0, 'DefaultFigureVisible', 'off');
h = figure('Position',[50,50,850,690]);

for it = 1:size(weights,1)

disp(['Starting iteration #', num2str(it)]);
    
%% varied parameters
varies(1).conxn = '(IC->IC)';
varies(1).param = 'trial';
varies(1).range = 1:40;

varies(2).conxn = 'R->C';
varies(2).param = 'gSYN';
varies(2).range = 0.18;

varies(3).conxn = 'C';
varies(3).param = 'noise';
varies(3).range = 1.3:0.1:1.4;

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
rcNetcon = weights(it,:)';
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
    study_dir = fullfile(pwd, 'run', 'grid-search', datetime, spatialConfig{1});
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
    [temp_perf(ns), temp_fr(ns), data(z).annot,~, data(z).VR] = ...
        mouse_network(study_dir,time_end,varies,netCons,plot_rasters,...
        calcDistances,data_spks,data_tau);
    end
    
    % average across multiple simulations if needed
    for vv = 1:nvaried
        inds = 1:nvaried:nvaried*nSims + (vv-1);
        temp = [temp_perf.R]; data(z).perf.R(:,vv) = mean(temp(:,inds),2);
        temp = [temp_perf.C]; data(z).perf.C(vv) = mean(temp(inds));
        temp = [temp_fr.R]; data(z).fr.R(:,vv) = mean(temp(:,inds),2);
        temp = [temp_fr.C]; data(z).fr.C(vv) = mean(temp(inds));
    end
    
    data(z).name = ICstruc(z).name;
end
save([pwd filesep 'run' filesep 'grid-search' filesep...
    datetime filesep 'summary_results_iteration' num2str(it) '.mat'],'data')
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
        ['RC weights = ',mat2str(rcNetcon)]};
        
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
        
    perfHistory(it,:,vv) = perf.CT;
    frHistory(it,:,vv) = mean(fr.CT);
    loss(it,vv) = MSE_clean(1);
end

end

end

%% need to manually code this when you change the parameter space for grid search

names = {'RC_{ipsi}','RC_{G}','RC_{U}','RC_{contra}'};

dims = length(inputChans);

% plot loss vs param
temp = weights(:,inputChans);
figure;
if dims < 2
    [temp,I] = sort(temp);
    plot(temp,loss(I));
    xlabel(names{inputChans}); ylabel('Loss');
elseif dims == 2
    sz1 = length(unique(temp(:,1)));
    sz2 = length(unique(temp(:,2)));
    
    X = reshape(temp(:,1),sz1,sz2)';
    Y = reshape(temp(:,2),sz1,sz2)';
    Z = reshape(loss,sz1,sz2)';
    surf(X,Y,Z);
    xlabel(names{inputChans(1)}); ylabel(names{inputChans(2)}); zlabel('Loss');
end

saveas(gcf,[filesep DirPart filesep 'MSE grid search.png'])
save([filesep DirPart filesep 'grid_search_results.mat'],'loss','weights','perfHistory');

% plot FR vs param
figure;
if dims < 2
    [temp,I] = sort(temp);
    plot(temp,squeeze(frHistory(:,:,I)));
    xlabel(names{inputChans}); ylabel('Mean clean FR');
    hold on; plot(temp,ones(size(temp))*mean(data_FR(targetIdx)))
elseif dims == 2
    sz1 = length(unique(temp(:,1)));
    sz2 = length(unique(temp(:,2)));
    
    X = reshape(temp(:,1),sz1,sz2)';
    Y = reshape(temp(:,2),sz1,sz2)';
    Z = reshape(frHistory,sz1,sz2)';
    surf(X,Y,Z);
    hold on;
    surf(X,Y,ones(size(Z))*mean(data_FR(targetIdx)))
    xlabel(names{inputChans(1)}); ylabel(names{inputChans(2)}); zlabel('Mean clean FR');
end
legend('Model','Data')

saveas(gcf,[filesep DirPart filesep 'MSE FR.png'])
