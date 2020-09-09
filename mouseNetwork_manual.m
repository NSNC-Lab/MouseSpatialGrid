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
dataCh = 31;

fitFlag = 1;    % if only want to generate clean and co-lcoated spots, =1

%%%%%%%%

subject = [extractBefore(perf_file,'_performance') '-Ch' num2str(dataCh)];
folder = ['Data-fitting' filesep subject];

load(fullfile(folder,'default_parameters.mat'));

ICdirPath = [ICdir filesep];

%% Initiate varies

varies = struct;

varies(1).conxn = '(IC->IC)';
varies(1).param = 'trial';
varies(1).range = 1:40;

varies(end+1).conxn = 'C';
varies(end).param = 'noise';
varies(end).range = 3;

Cnoise2 = 4; % additional noise, so colocated noise = Cnoise2 + varies(2).range

varies(end+1).conxn = 'R->C';
varies(end).param = 'gSYN1';
varies(end).range = 0;

varies(end+1).conxn = 'R->C';
varies(end).param = 'gSYN2';
varies(end).range = 0;

varies(end+1).conxn = 'R->C';
varies(end).param = 'gSYN3';
varies(end).range = 0.12;

varies(end+1).conxn = 'R->C';
varies(end).param = 'gSYN4';
varies(end).range = 0;

nvaried = {varies(2:end).range};
nvaried = prod(cellfun(@length,nvaried));

plot_rasters = 1;

ICstruc = dir([ICdirPath '*.mat']);

if fitFlag == 1
    subz = find(cellfun(@(x) strcmp(x(2),x(4)),{ICstruc.name})); % co-located cases
    subz = cat(2,subz,find(contains({ICstruc.name},'m0.mat'))); % sXm0 (target only) cases
    detail = '-manual-fit';
else
    subz = find(~contains({ICstruc.name},'s0'));  % all spots
    detail = '-full-grid';
end

[simdata,DirPart] = mouseNetwork_initialize(varies,Cnoise2,ICdirPath,spks_file,dataCh,plot_rasters,...
    folder,detail,subz,[]);

%% performance grids for 4D search

temp = {simdata.name};
temp(cellfun('isempty',temp)) = {'empty'}; %label empty content

targetIdx = find(contains(temp,'m0') & ~strcmp(temp,'s0m0.mat'));
colocIdx = find((cellfun(@(x) x(2),temp) == cellfun(@(x) x(4),temp)));
perf_clean = []; model_FR = [];

% clean performance and FR
for i = 1:length(targetIdx)
    perf_clean(:,i) = simdata(targetIdx(i)).perf.C;
    model_FR(:,i) = simdata(targetIdx(i)).fr.C;
end

perf_masked = []; model_FR_masked = [];
% mixed performance and FR
if ~isempty(colocIdx)
for i = 1:length(colocIdx)
    perf_masked(:,i) = simdata(colocIdx(i)).perf.C;
    model_FR_masked(:,i) = simdata(colocIdx(i)).fr.C;
end
end

[~,MSE_clean_perf] = calcModelLoss(perf_clean,data_perf(1:4)');
% calculate loss for FR as a percentage of model FR
[~,MSE_clean_FR] = calcModelLoss(100*(model_FR)./data_FR,100*data_FR./data_FR);

if ~isempty(colocIdx)
    [~,MSE_masked_perf] = calcModelLoss(perf_masked,data_perf([5 10 15 20])');
    % calculate loss for FR as a percentage of model FR
    [~,MSE_masked_FR] = calcModelLoss(100*(model_FR_masked)./data_FR_masked([1 6 11 16]),100*data_FR_masked([1 6 11 16])./data_FR_masked([1 6 11 16]));
    loss = MSE_clean_perf(:,1) + MSE_clean_FR(:,1) + MSE_masked_perf(:,1) + MSE_masked_FR(:,1);
else
    loss = MSE_clean_perf(:,1) + MSE_clean_FR(:,1);
end

% make spatial grid
if length(subz) ~= 20
    makeGrids_fitting(simdata,DirPart,data_perf,data_FR,data_FR_masked(1:5:end));
else
    makeGrids(simdata,DirPart,data_perf,data_FR,data_FR_masked(1:5:end),loss)
end

