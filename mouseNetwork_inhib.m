% Run this part after mouseNetwork_premanual

% Mess with this code to manually fit models to data

clearvars;
close all;

addpath('mechs')
addpath('dependencies')
addpath('eval_scripts')
addpath('genlib')
addpath(genpath('dynasim'))
addpath(genpath('ICSimStim'))

%%%%%%%% 

spks_file = '9-21-2016_0dB_removed_trials_cleaned(-1,4).mat';
perf_file = '9-21-2016_0dB_removed_trials_performance.mat';
dataCh = 25;

%%%%%%%%

subject = [extractBefore(perf_file,'_performance') '-Ch' num2str(dataCh)];
folder = ['Data-fitting' filesep subject];

load(fullfile(folder,'default_parameters.mat'));

ICdirPath = [ICdir filesep];

%% Initiate varies

varies = struct;

varies(1).conxn = '(IC->IC)';
varies(1).param = 'trial';
varies(1).range = 1:40;%[1:10 21:30]; %

varies(end+1).conxn = 'C';
varies(end).param = 'noise';
varies(end).range = 2;

Cnoise2 = 0.5; % additional noise, so colocated noise = Cnoise2 + varies(2).range

varies(end+1).conxn = 'R->C';
varies(end).param = 'gSYN1';
varies(end).range = 0.18;

varies(end+1).conxn = 'R->C';
varies(end).param = 'gSYN2';
varies(end).range = 0;

varies(end+1).conxn = 'R->C';
varies(end).param = 'gSYN3';
varies(end).range = 0;

varies(end+1).conxn = 'R->C';
varies(end).param = 'gSYN4';
varies(end).range = 0;

xrNetcon = zeros(4);
xrNetcon(1,4) = 1;

varies(end+1).conxn = 'X->R';
varies(end).param = 'gSYN';
varies(end).range = 0.14;

plot_rasters = 1;

ICstruc = dir([ICdirPath '*.mat']);
fitFlag = 0;
if fitFlag == 1
    subz = find(cellfun(@(x) strcmp(x(2),x(4)),{ICstruc.name})); % co-located cases
    subz = cat(2,subz,find(contains({ICstruc.name},'m0.mat'))); % target only cases
    detail = '-manual-fit-inhib';
else
    subz = find(~contains({ICstruc.name},'s0'));  % all spots
    detail = '-full-grid-inhib';
end

[simdata,DirPart] = mouseNetwork_inhib_initialize(varies,xrNetcon,Cnoise2,ICdirPath,...
    spks_file,dataCh,plot_rasters,folder,detail,subz,[]);

%% performance grids for 4D search

temp = {simdata.name};
temp(cellfun('isempty',temp)) = {'empty'}; %label empty content

targetIdx = find(contains(temp,'m0') & ~strcmp(temp,'s0m0.mat'));
colocIdx = find((cellfun(@(x) x(2),temp) == cellfun(@(x) x(4),temp)));
perf_clean = [];

if plot_rasters
    rasters = dir([DirPart filesep '*.fig']);
    temp = {rasters.name};
    rasters(strcmp(cellfun(@(x) x(4),temp,'uniformoutput',false),'0')) = [];
    
    locs = {'90°','45°','0°','-90°'};
    
    % plot PSTHs of clean trials over co-located trials
    for i = 1:length(targetIdx)
        for m = 1:length(targetIdx)
            cleanPSTH = simdata(targetIdx(i)).PSTH;
            mixedPSTH = simdata(targetIdx(i)+m).PSTH;
            t_vec1 = 0:0.02:(length(cleanPSTH)-1)*0.02;
            t_vec2 = 0:0.02:(length(mixedPSTH)-1)*0.02;
            openfig([DirPart filesep rasters((i-1)*4 + m).name],'visible');
            subplot(2,3,6);
            hold off;
            plot(t_vec2,50*mixedPSTH/size(varies(1).range,2));
            hold on;
            plot(t_vec1,50*cleanPSTH/size(varies(1).range,2));
            legend('Mixed','Clean');
            ylim([0 150]);
            savefig([DirPart filesep rasters((i-1)*4 + m).name]);
            saveas(gcf,[DirPart filesep rasters((i-1)*4 + m).name(1:end-4) '.png']);
            close;
        end
    end
end


% clean performance and FR
for i = 1:length(targetIdx)
    perf_clean(:,i) = simdata(targetIdx(i)).perf.C;
end

perf_masked = [];
% mixed performance and FR
if ~isempty(colocIdx)
for i = 1:length(colocIdx)
    perf_masked(:,i) = simdata(colocIdx(i)).perf.C;
end
end

[~,MSE_clean_perf] = calcModelLoss(perf_clean,data_perf(1:4)');
[~,MSE_masked_perf] = calcModelLoss(perf_masked,data_perf([5:5:end])');
loss = MSE_clean_perf(:,1) + MSE_masked_perf(:,1);

% make spatial grid
if length(subz) ~= 20
    makeGrids_fitting(simdata,DirPart,data_perf,data_FR,data_FR_masked(1:5:end));
else
    makeGrids(simdata,DirPart,data_perf,data_FR,data_FR_masked(1:5:end),loss)
end

%[] = find(xrNetcon == 1);
