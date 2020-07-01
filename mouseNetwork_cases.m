clearvars;
close all;

addpath('mechs')
addpath('dependencies')
addpath('eval_scripts')
addpath('genlib')
addpath(genpath('dynasim'))

ch = 25;

load('9-21-2016_0dB_removed_trialscleaned(-1,4).mat','Spks_clean','Spks_masked');
load('9-21-2016_0dB_removed_trials_performance.mat','Max','max_masked');
data_perf = [Max{ch}(:);max_masked{ch}(:)];

researchDrive = 'MiceSpatialGrids/';
ICdir = [researchDrive 'ICStim/Mouse/full_grids//BW_0.009 BTM_3.8 t0_0.1 phase0.4985//s30_STRFgain4.50_20200623-171359'];

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

varies(end+1).conxn = 'R->C';
varies(end).param = 'gSYN';
varies(end).range = 0.21;

variedParam = 'R->CgSYN';

% enter each case here as rows
rcCases = eye(4);
rcCases([3 4],:) = rcCases([3 4],:)*0.6;

% manually edit these too
cases = {'ipsi','gauss','ushaped','contra'};
caseType = 'single cases';

datetime = datestr(now,'yyyymmdd-HHMMSS');

% transpose rcCases for dynasim runs
rcCases = rcCases';

for nc = 3:size(rcCases,2)
%% netcons
nCells = 4; %synchronise this variable with mouse_network

% -90, 0, 45, 90º
irNetcon = zeros(nCells);

srNetcon = zeros(nCells);

rcNetcon = rcCases(:,nc);

netCons.irNetcon = irNetcon;
netCons.srNetcon = srNetcon;
netCons.rcNetcon = rcNetcon;
%% Initialize variables
plot_rasters = 0;

nvaried = {varies(2:end).range};
nvaried = prod(cellfun(@length,nvaried));
diagConfigs = [6,12,18,24];

% set(0, 'DefaultFigureVisible', 'off')
h = figure('Position',[50,50,850,690]);

subz = find(contains({ICstruc.name},'m0.mat')); % sXm0 (target only) cases
% subz = 5:length({ICstruc.name});    % all cases except masker-only
loc = 1;
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
            data_spks = squeeze(Spks_clean{ch}(:,5-str2double(spatialConfig{1}(2)),:));
            dataFR(z) = mean(mean(cellfun(@(x) sum(x >= 0 & x < time_end/1000),data_spks)))/time_end*1000;
        else
            data_spks = squeeze(Spks_masked{ch}(:,str2double(spatialConfig{1}(2)),str2double(spatialConfig{1}(4)),:));
        end
    else
        data_spks = [];
    end
    
    [~, temp, ~, ~] = ...
        mouse_network(study_dir,time_end,varies,netCons,plot_rasters,data_spks);
    fr(nc,loc) = temp.C;
    loc = loc + 1;
end
% save([pwd filesep 'run' filesep, datetime filesep 'summary_results.mat'],'data')
end

%% plots

figure;

colors = {'k','r','g','b'};
    
% Show FR vs azimuth
hold on
for nc = 1:size(rcCases,2)
    plot([-90 0 45 90],fr(nc,:),colors{nc},'linewidth',2);
end

legend(cases)
for nc = 1:size(rcCases,2)
    plot([-90 0 45 90],mean(fr(nc,:))*ones(1,4),[':',colors{nc}],'linewidth',2);
end
xlabel('Azimuth');
ylabel('Clean FR (Hz)')
set(gca,'xdir','reverse');
ylim([0 50]);
xticks([-90,0:45:90]);
%title(['RC weights = ' mat2str(round(100*rcNetcon)/100)])

Dirparts = strsplit(study_dir, filesep);
DirPart = fullfile(Dirparts{1:end-1});

saveas(gca,[filesep DirPart filesep 'clean FR ' caseType '.png']);

%clf;

