clearvars;
close all;

addpath('mechs')
addpath('dependencies')
addpath('eval_scripts')
addpath('genlib')
addpath(genpath('dynasim'))

%%%%%%%% start of user inputs

plot_grids = 0;  % plot spatial grids
plot_distances = 0;  % plot VR distances
plot_rasters = 1;   % plot rasters

data_spks_file = '9-21-2016_0dB_removed_trialscleaned(-1,4).mat';
data_perf_file = '9-21-2016_0dB_removed_trials_performance.mat';
dataCh = 25;

%%%%%%%% end of user inputs

subject = [extractBefore(data_perf_file,'_performance') '-Ch' num2str(dataCh)];
folder = 'Data-fitting';

%% load experimental data to optimize model to

nCells = 4; %synchronise this variable with mouse_network

load(data_spks_file,'Spks_clean','Spks_masked');
load(data_perf_file,'Max','max_masked','opt_tau','opt_tau_masked');
data_perf = [Max{dataCh}(:);max_masked{dataCh}(:)];

% use tau of best clean location
[~,best_loc] = max(Max{dataCh});
data_tau = opt_tau{dataCh}(best_loc);

%% calculate cortical noise based on spontaneous FR

% data_FR arranged ipsi->contra
for s = 1:4
    data_spks = squeeze(Spks_clean{dataCh}(:,5-s,:));
    FR_r0(s) = mean(cellfun(@(x) sum(x < 0 | x >= 3),data_spks),'all')/2;
    data_FR(s) = mean(cellfun(@(x) sum(x >= 0 & x < 3),data_spks),'all')/3;
end
load('Cnoise_vs_FR0.mat','fit');
Cnoise = (mean(FR_r0)-fit(1))/fit(2);

%% calculate best STRF gain data based on single channel case

[~,best_chan] = max(data_perf(1:4));
best_chan = 5-best_chan;    % azimuth indexes are flipped b/w model and data

varies = struct;

varies(1).conxn = '(IC->IC)';
varies(1).param = 'trial';
varies(1).range = 1:40;

varies(end+1).conxn = 'C';
varies(end).param = 'noise';
varies(end).range = Cnoise;

for c = 1:4
    varies(end+1).conxn = 'R->C';
    varies(end).param = ['gSYN' num2str(c)];
    if c == best_chan
        varies(end).range = 0.21;
    else
        varies(end).range = 0;
    end
end

% extract all STRF gains
gains = dir('/Users/jionocon/Documents/MATLAB/Spatial-grid-simulations/MiceSpatialGrids/ICStim/Mouse/full_grids/BW_0.009 BTM_3.8 t0_0.1 phase0.4985');
gains(~contains({gains.name},'s30')) = [];

model_FR = [];
for g = 1:length(gains)-1
    ICdir = [gains(g).folder filesep gains(g).name];
    ICdirPath = [ICdir filesep];
    ICstruc = dir([ICdirPath '*.mat']);
    
    data = mouseNetwork_initialize(varies,ICstruc,ICdirPath,Spks_clean,...
        Spks_masked,dataCh,data_tau,plot_distances,plot_rasters,folder,subject,'-STRF-fitting');
    
    temp = {data.name};
    temp(cellfun('isempty',temp)) = {'empty'}; %label empty content
    
    targetIdx = find(contains(temp,'m0') & ~strcmp(temp,'s0m0.mat'));
    %maskerIdx = find(contains(temp,'s0') & ~strcmp(temp,'s0m0.mat'));
    %mixedIdx = find(~contains(temp,'m0') & ~contains(temp,'s0') & ~contains(temp,'empty'));
    for i = 1:length(targetIdx)
        model_FR(g,i) = data(targetIdx(i)).fr.C;
    end
end
% find best gain closest to average clean trial FR
[~,ind] = min(abs(mean(model_FR,2) - mean(data_FR)));

ICdir = [gains(ind).folder filesep gains(ind).name];
ICdirPath = [ICdir filesep];
ICstruc = dir([ICdirPath '*.mat']);

%% single channel cases
gsyn_range = 0.09:0.03:0.21;

single_channel_folders = [];

for cases = 4
    
ranges = cell(1,4);
for n = 1:nCells
    if n == cases
        ranges{n} = gsyn_range;
    else
        ranges{n} = 0;
    end
end
    
%% varied parameters

varies = struct;

varies(1).conxn = '(IC->IC)';
varies(1).param = 'trial';
varies(1).range = 1:40;

varies(end+1).conxn = 'C';
varies(end).param = 'noise';
varies(end).range = Cnoise;

varies(end+1).conxn = 'R->C';
varies(end).param = 'gSYN1';
varies(end).range = ranges{1};

varies(end+1).conxn = 'R->C';
varies(end).param = 'gSYN2';
varies(end).range = ranges{2};

varies(end+1).conxn = 'R->C';
varies(end).param = 'gSYN3';
varies(end).range = ranges{3};

varies(end+1).conxn = 'R->C';
varies(end).param = 'gSYN4';
varies(end).range = ranges{4};

nvaried = {varies(2:end).range};
nvaried = prod(cellfun(@length,nvaried));

[data,~,single_channel_folders{cases}] = mouseNetwork_initialize(varies,ICstruc,ICdirPath,Spks_clean,...
    Spks_masked,dataCh,data_tau,plot_distances,plot_rasters,folder,subject,'-one-channel');

%% performance grids

if plot_grids
    makeGrids(data,varies,single_channel_folders{cases},nvaried,data_perf,data_FR);
end

temp = {data.name};
temp(cellfun('isempty',temp)) = {'empty'}; %label empty content

targetIdx = find(contains(temp,'m0') & ~strcmp(temp,'s0m0.mat'));
%maskerIdx = find(contains(temp,'s0') & ~strcmp(temp,'s0m0.mat'));
%mixedIdx = find(~contains(temp,'m0') & ~contains(temp,'s0') & ~contains(temp,'empty'));
perf = []; model_FR = [];
for i = 1:length(targetIdx)
    perf(:,5-i) = data(targetIdx(i)).perf.C;
    model_FR(:,5-i) = data(targetIdx(i)).fr.C;
end

[~,MSE_clean] = calcModelPerf(perf,data_perf(1:4)');

loss = MSE_clean(:,1);
frdiffs = abs(mean(model_FR,2) - mean(data_FR(data_FR ~= 0)));
within_thresh = frdiffs <= 5;

% no wts are within FR thresh, just get the absolute minimum
if sum(within_thresh) == 0
   within_thresh = frdiffs == min(frdiffs);
end

[minloss(cases),i] = min(loss(within_thresh));
temp = find(within_thresh);
best_iteration(cases) = temp(i);

best_wt(cases) = gsyn_range(best_iteration);
best_perf(cases,:) = perf(best_iteration,:);

end

% determine best of first channel cases

[onechan_loss,best_channel] = min(minloss);
best_gSYN = best_wt(best_channel);

onechan_perf = best_perf(best_channel,:);

onechan_results.loss = onechan_loss;
onechan_results.best_channel = best_channel;
onechan_results.best_gSYN = best_gSYN;
onechan_results.perf = onechan_perf;
onechan_results.folders = single_channel_folders;
onechan_results.best_iteration = best_iteration;
onechan_results.minloss = minloss;

%% second channel cases
minloss = [];
best_wt = [];
twochan_perf = [];

two_channel_folders = [];
for cases = find(1:4 ~= best_channel)
     
% % -90, 0, 45, 90º

best_gsyn_range = best_gSYN-0.02:0.01:best_gSYN+0.02;

second_gsyn_range = 0.03:0.03:best_gSYN;

ranges = cell(1,4);
for n = 1:nCells
    if n == cases
        ranges{n} = second_gsyn_range;
    elseif n == best_channel
        ranges{n} = best_gsyn_range;
    else
        ranges{n} = 0;
    end
end

%% varied parameters
varies = struct;

varies(1).conxn = '(IC->IC)';
varies(1).param = 'trial';
varies(1).range = 1:40;

varies(end+1).conxn = 'C';
varies(end).param = 'noise';
varies(end).range = Cnoise;

varies(end+1).conxn = 'R->C';
varies(end).param = 'gSYN1';
varies(end).range = ranges{1};

varies(end+1).conxn = 'R->C';
varies(end).param = 'gSYN2';
varies(end).range = ranges{2};

varies(end+1).conxn = 'R->C';
varies(end).param = 'gSYN3';
varies(end).range = ranges{3};

varies(end+1).conxn = 'R->C';
varies(end).param = 'gSYN4';
varies(end).range = ranges{4};

nvaried = {varies(2:end).range};
nvaried = prod(cellfun(@length,nvaried));

[data,~,two_channel_folders{cases}] = mouseNetwork_initialize(varies,ICstruc,ICdirPath,Spks_clean,...
    Spks_masked,dataCh,data_tau,plot_distances,plot_rasters,folder,subject,'-two-channel');

%% performance grids

if plot_grids
    makeGrids(data,varies,two_channel_folders{cases},nvaried,data_perf,data_FR);
end

temp = {data.name};
temp(cellfun('isempty',temp)) = {'empty'}; %label empty content

targetIdx = find(contains(temp,'m0') & ~strcmp(temp,'s0m0.mat'));
%maskerIdx = find(contains(temp,'s0') & ~strcmp(temp,'s0m0.mat'));
%mixedIdx = find(~contains(temp,'m0') & ~contains(temp,'s0') & ~contains(temp,'empty'));
perf = []; model_FR = [];
for i = 1:length(targetIdx)
    perf(:,5-i) = data(targetIdx(i)).perf.C;
    model_FR(:,5-i) = data(targetIdx(i)).fr.C;
end

[~,MSE_clean] = calcModelPerf(perf,data_perf(1:4)');

loss = MSE_clean(:,1);
frdiffs = abs(mean(model_FR,2) - mean(data_FR(data_FR ~= 0)));
within_thresh = frdiffs <= 5;

% no wts are within FR thresh, just get the absolute minimum
if sum(within_thresh) == 0
   within_thresh = frdiffs == min(frdiffs);
end

[minloss(cases),i] = min(loss(within_thresh));
temp = find(within_thresh);
best_iteration(cases) = temp(i);

[A,B] = meshgrid(best_gsyn_range,second_gsyn_range);
c = cat(2,A',B');
all_wt_combos = reshape(c,[],2);

best_two_wts(cases,:) = all_wt_combos(best_iteration,:);
twochan_perf(cases,:) = perf(best_iteration,:);

end

minloss(minloss == 0) = [];
best_two_wts(best_two_wts == 0) = [];
best_two_wts = reshape(best_two_wts,3,2);
two_channel_folders(cellfun(@isempty,two_channel_folders)) = [];

twochan_perf(twochan_perf == 0) = [];
twochan_perf = reshape(twochan_perf,3,4);

% determine if adding second channel improves one-channel case

decreased_loss = (minloss - onechan_loss) < 0;
for cases = 1:3
   sig_perf(cases) = ttest(onechan_perf,twochan_perf(cases,:),'Tail','left'); 
end

if sum(sig_perf) == 0 && sum(decreased_loss) == 0
    disp('One-channel case is best');
    singleFlag = 1;
else
    disp(['Two-channel case ' num2str(find(decreased_loss)) ' outperformed single channel']);
    singleFlag = 0;
end

twochan_results.best_gSYNs = best_two_wts;
twochan_results.decreased_loss = decreased_loss;
twochan_results.sig_perf = sig_perf;
twochan_results.folders = two_channel_folders;
twochan_results.best_iteration = best_iteration;
twochan_results.minloss = minloss;

save([folder filesep subject 'fitting_results.mat'],...
    'twochan_results','Cnoise','onechan_results');

%% generate best grid

if singleFlag
    best_case_folder = single_channel_folders{best_channel};
else
    best_case_folder = two_channel_folders{decreased_loss};
end

addpath(genpath([cd filesep folder]));

load([best_case_folder filesep 'summary_results.mat'])
nvaried = {varies(2:end).range};
nvaried = prod(cellfun(@length,nvaried));

makeGrids(data,varies,[cd filesep folder filesep subject],nvaried,data_perf,data_FR)