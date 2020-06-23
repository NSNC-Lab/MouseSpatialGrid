%% Main script for calling mouse data simulation network
% 2019-08-14 plotting now handles multiple varied parameters
% 2019-09-11 moved plotting to a separate script plotPerfGrid.m
%            mouseNetwork returns both R and C performance
%            now plot both R and C performance grids
%
% to do:
%   add ability to adjust RC netcon in main code
clearvars;
close all;

addpath('mechs')
addpath('dependencies')
addpath('eval_scripts')
addpath('genlib')
addpath(genpath('dynasim'))

ch = 28;

load('10_21_17_0dB_rawData(-1,4).mat','Spks_clean','Spks_masked');

load('10_21_17_0dB_performance_raw(-1,4).mat','Max','max_masked');
data_perf = [Max{ch}(:);max_masked{ch}(:)];

researchDrive = 'MiceSpatialGrids/';
ICdir = [researchDrive 'ICStim/Mouse/full_grids//BW_0.009 BTM_3.8 t0_0.1 phase0.4985//s30_STRFgain4.50_20200622-150507'];

% TemporalBW 0.009s/s30_gain0.5_20200602-143614 - sharpened sigmoids

ICdirPath = [ICdir filesep];
ICstruc = dir([ICdirPath '*.mat']);
if isempty(ICstruc), error('empty data directory'); end

%% varied parameters
varies(1).conxn = '(IC->IC)';
varies(1).param = 'trial';
varies(1).range = 1:40;

varies(end+1).conxn = 'C';
varies(end).param = 'noise';
varies(end).range = 0;

% varies(end+1).conxn = 'IC->X';
% varies(end).param = 'gSYN';
% varies(end).range = 0.6;

% rcCases = ones(4,4);
% rcCases = rcCases - diag(ones(4,1));
%rcCases = rcCases';
rcCases = [1;1;1;1];
for nc = 1:size(rcCases,2)
%% netcons
nCells = 4; %synchronise this variable with mouse_network

% -90, 0, 45, 90º
% x-channel inhibition
% irNetcon = diag(ones(1,nCells))*0.1;
irNetcon = zeros(nCells);
% irNetcon(4,2) = 0.3;
% irNetcon(1,4) = 1;
% irNetcon(4,2) = 0.2;

srNetcon = zeros(nCells);

rcNetcon = rcCases(:,nc); %add this as input to mouse_network
% make rnNetcon have variable weights (instead of zeros)

netCons.irNetcon = irNetcon;
netCons.srNetcon = srNetcon;
netCons.rcNetcon = rcNetcon;
%% Initialize variables
plot_rasters = 1;

nvaried = {varies(2:end).range};
nvaried = prod(cellfun(@length,nvaried));
diagConfigs = [6,12,18,24];
datetime=datestr(now,'yyyymmdd-HHMMSS');

%set(0, 'DefaultFigureVisible', 'off')
h = figure('Position',[50,50,850,690]);

subz = find(contains({ICstruc.name},'m0.mat')); % sXm0 (target only) cases
%subz = 5:length({ICstruc.name});    % all cases except masker-only
for z = subz(1)
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

    % load spikes from data
    if ~contains(spatialConfig{1},'s0')
        if strcmp(spatialConfig{1}(4),'0')
            data_spks = squeeze(Spks_clean{ch}(:,str2double(spatialConfig{1}(2)),:));
        else
            data_spks = squeeze(Spks_masked{ch}(:,str2double(spatialConfig{1}(2)),str2double(spatialConfig{1}(4)),:));
        end
    else
        data_spks = [];
    end
    
    % call network
    h.Name = ICstruc(z).name;
    time_end = size(spks,3);
    [data(z).perf, data(z).fr, data(z).optTau, data(z).annot] = ...
        mouse_network(study_dir,time_end,varies,netCons,plot_rasters,data_spks);
    data(z).name = ICstruc(z).name;
end
% save([pwd filesep 'run' filesep, datetime filesep 'summary_results.mat'],'data')
% 
% save(fullfile(pwd,'run',datetime,filesep,'data.mat'),'data');
% 
% %% performance grids
% % performance vector has dimensions [numSpatialChan,nvaried]
% neurons = {'Ipsi. sigmoid','Gaussian','U','Cont. sigmoid'};
% 
% temp = {data.name};
% temp(cellfun('isempty',temp)) = {'empty'}; %label empty content
% 
% targetIdx = find(contains(temp,'m0') & ~strcmp(temp,'s0m0.mat'));
% maskerIdx = find(contains(temp,'s0') & ~strcmp(temp,'s0m0.mat'));
% mixedIdx = setdiff(1:length(temp),[targetIdx,maskerIdx,1]);
% %mixedIdx = find(~contains(temp,'m0') & ~contains(temp,'s0') & ~contains(temp,'empty'));
% textColorThresh = 70;
% numSpatialChan = 4;
% 
% width=13; hwratio=0.8;
% x0=.08; y0=.08;
% dx=.04; dy=.04;
% lx=.125; ly=.125/hwratio;
% 
% sigma = 30;
% 
% x=-108:108;
% tuningcurve=zeros(4,length(x));
% % tuningcurve(1,:)= (sigmf(-x,[0.07 30]))*rcNetcon(1); %flipped sigmodial
% % tuningcurve(2,:)= (gaussmf(x,[sigma,0]))*rcNetcon(2); %guassian
% % tuningcurve(3,:)= (1-gaussmf(x,[sigma,0]))*rcNetcon(3); %U shaped gaussian
% % tuningcurve(4,:)= (sigmf(x,[0.07 30]))*rcNetcon(4); %sigmodial
% % tuningcurve(sum(tuningcurve,2)==0,:) = [];
% 
%         ono = load('ono_curves.mat');
%         s = ono.s;
%         u = ono.u;
%         g = ono.g;
%         tuningcurve(1,:) = s(ono.By1E,x)/100 * rcNetcon(1);
%         tuningcurve(2,:) = fliplr(u(ono.By2E,x))/100 * rcNetcon(2);
%         tuningcurve(3,:) = fliplr(g(ono.By3E,x))/100 * rcNetcon(3);
%         tuningcurve(4,:) = fliplr(s(ono.By1E,x))/100 * rcNetcon(4);
% 
% figstr = 'SpatialGrid vary ';
% 
% % if more than one parameter is varied
% if sum(cellfun(@length,{varies(2:end).range}) > 1) > 1
%     
%     multiFlag = 1;
%     
%     temp = find((cellfun(@length,{varies.range}) > 1) == 1);
%     temp(1) = [];
%     variedConxns = {varies(temp).conxn};
%     variedParams = {varies(temp).param};
%     variedRanges = {varies(temp).range};
%     
%     for i = 1:length(temp)
%         figstr = cat(2,figstr,[variedConxns{i},variedParams{i},'%0.3f, ']);
%     end
%     figstr(end-1:end) = []; % delete extra ', '
%     
%     [A,B] = meshgrid(variedRanges{1},variedRanges{2});
%     c = cat(2,A',B');
%     paramPairs = reshape(c,[],2);
% else
%     multiFlag = 0;
%     figstr = cat(2,figstr,variedParam,'%0.3f');
% end
% 
% h = figure('position',[200 200 1600 600]);
% 
% for vv = 1:nvaried
%     if ~isempty(mixedIdx)
%         % mixed config cases
%         for i = 1:length(mixedIdx)
%             perf.IC(i,:) = data(mixedIdx(i)).perf.IC(:,vv);
%             fr.IC(i,:) = data(mixedIdx(i)).fr.IC(:,vv);
%             
%             perf.R(i,:) = data(mixedIdx(i)).perf.R(:,vv);
%             fr.R(i,:) = data(mixedIdx(i)).fr.R(:,vv);
%             
%             perf.C(i) = data(mixedIdx(i)).perf.C(vv);
%             fr.C(i) = data(mixedIdx(i)).fr.C(vv);
%         end
%         h = figure;
%         figuresize(width, width*hwratio,h, 'inches')
% 
%         % IC and relay neurons
%         for nn = 1:numSpatialChan
%             if nn > 1
%                 subplotloc = nn+1;
%             else
%                 subplotloc = nn;
%             end
%             posVec = [1-(x0+subplotloc*(dx+lx)) y0+dy+ly lx ly];
%             subplot('Position',posVec)
%             plotPerfGrid(perf.IC(end:-1:1,nn),fr.IC(end:-1:1,nn),[],textColorThresh);
%             
%             if nn == 4
%                 ylabel('IC')
%             else
%                 ylabel('');
%             end
%             xlabel(neurons{nn});
%                             
%             posVec = [1-(x0+subplotloc*(dx+lx)) y0+2*(dy+ly) lx ly];
%             subplot('Position',posVec)
%             plotPerfGrid(perf.R(end:-1:1,nn),fr.R(end:-1:1,nn),[],textColorThresh);
%             if nn == 4
%                 ylabel('R')
%             else
%                 ylabel('');
%             end
%             xlabel('');   
%         end
%         
%         % C neuron
% 
%         posVec = [1-(x0+3*(dx+lx)) y0+3*(dy+ly) lx ly];
%         subplot('Position',posVec);
%         plotPerfGrid(perf.C(end:-1:1)',fr.C(end:-1:1)',[],textColorThresh);
%         xlabel('')
%         ylabel('C')
%     end
%     
%     subplot('position',[0.3 y0+dy/2 0.4 4*ly/5]); 
%     plot(x,tuningcurve','b','linewidth',1);
%     hold on;
%     plot(x,sum(tuningcurve',2),'k','linewidth',2);
%     ylim([0 max(sum(tuningcurve',2))+0.2]);
%     xlim([x(1) x(end)]);
%     set(gca,'xdir','reverse');
%     xticks([-90,0:45:90]);
%     xlabel('azimuth');
%     
%     % C neuron; target or masker only cases
%     perf.CT = zeros(1,4);
%     perf.CM = zeros(1,4);
% 
%     fr.CT = zeros(1,4);
%     fr.CM = zeros(1,4);
%     if ~isempty(targetIdx)
%         for i = 1:length(targetIdx)
%             perf.CT(i) = data(targetIdx(i)).perf.C(vv);
%             fr.CT(i) = data(targetIdx(i)).fr.C(vv);
%         end
%     end    
%     if ~isempty(maskerIdx)
%         for i = 1:length(targetIdx)
%             perf.CM(i) = data(maskerIdx(i)).perf.C(vv);
%             fr.CM(i) = data(maskerIdx(i)).fr.C(vv);
%         end
%     end
%     
% %     posVec = [posVec(1) posVec(2)+ly+dy/4 lx ly/2];
% %     subplot('Position',posVec)
% %     plotPerfGrid([perf.CT(end:-1:1);perf.CM(end:-1:1)],[fr.CT(end:-1:1);fr.CM(end:-1:1)],[],textColorThresh);
% 
%     posVec = [posVec(1) posVec(2)+ly+dy/4 lx ly/4];
%     subplot('Position',posVec)
%     plotPerfGrid(perf.CT(end:-1:1),[fr.CT(end:-1:1);fr.CM(end:-1:1)],[],textColorThresh);
%     
% %     % clean firing rate info
% %     for nn = 1:4
% %         str = cellstr(num2str(round(fr.CT(end+1-nn))));
% %         annotation('textbox',[posVec(1)-0.005+(lx/4+0.005)*(nn-1) posVec(2)+ly/3+0.02 lx/4 ly/4],...
% %            'string',str,...
% %            'FitBoxToText','on',...
% %            'LineStyle','none')
% %     end
%     
%     % simulation info
%     annotation('textbox',[.65 .8 .2 .1],...
%            'string',data(1).annot(vv,3:end),...
%            'FitBoxToText','on',...
%            'LineStyle','none')
% 
% %     % calculate error and correlation with data
% %     model_perf = fliplr(reshape(perf.C,[4 4]));
% %     model_perf = [perf.CT';model_perf(:)];
% %     [model_corr,model_error] = calcModelPerf(model_perf,data_perf);
% % 
% %     str = {sprintf('cc = %0.3f',model_corr),sprintf('deviation = %0.1f ± %0.1f',model_error)};
% %        
% %     % correlation and error stuff
% %     annotation('textbox',[.65 .7 .2 .1],...
% %            'string',str,...
% %            'FitBoxToText','on',...
% %            'LineStyle','none')
%        
%     % save grid
%     Dirparts = strsplit(study_dir, filesep);
%     DirPart = fullfile(Dirparts{1:end-1});
%     if multiFlag == 1
%         saveas(gca,[filesep DirPart filesep sprintf(figstr,paramPairs(vv,1),paramPairs(vv,2)) '.tiff'])
%     else
%         saveas(gca,[filesep DirPart filesep sprintf(figstr,varies(end).range(vv)) '.tiff'])
%     end
%     %save([filesep DirPart filesep 'performance.mat'],'perf');
%     clf
% end
% 
% mouseNetwork_plotFiring;

end
