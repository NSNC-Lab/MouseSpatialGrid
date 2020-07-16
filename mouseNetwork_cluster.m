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
plot_rasters = 0;   % plot rasters

data_spks_file = '03_30_18_0dB_cleaned(-1,4).mat';
data_perf_file = '03_30_18_0dB_performance.mat';
dataCh = 23;

%%%%%%%% end of user inputs

%% load experimental data to optimize model to

nCells = 4; %synchronise this variable with mouse_network

load(data_spks_file,'Spks_clean','Spks_masked');
load(data_perf_file,'Max','max_masked','opt_tau','opt_tau_masked');
data_perf = [Max{dataCh}(:);max_masked{dataCh}(:)];

% use tau of best clean location
[~,best_loc] = max(Max{dataCh});
data_tau = opt_tau{dataCh}(best_loc);

%% calculate which STRF gain to use

for s = 1:4
    data_spks = squeeze(Spks_clean{dataCh}(:,5-s,:));
    time_end = 2971;
    data_FR(s) = mean(mean(cellfun(@(x) sum(x >= 0 & x < time_end/1000),data_spks)))/time_end*1000;
end
data_FR = mean(data_FR);
% from strfFR code: all_fit is linear fit to STRF gain vs FR
all_fit = [9.35093660665019; 5.90086696599368];

% extract all STRF gains
gains = dir('/Users/jionocon/Documents/MATLAB/Spatial-grid-simulations/MiceSpatialGrids/ICStim/Mouse/full_grids/BW_0.009 BTM_3.8 t0_0.1 phase0.4985');
gains(~contains({gains.name},'s30')) = [];
gain_nums = extractBetween({gains.name},'STRFgain','_20');
gain_nums = cellfun(@str2num,gain_nums);

gain1 = (data_FR-all_fit(1))/all_fit(2);
[~,n] = min(abs(gain_nums-gain1));
ICdir = [gains(n).folder filesep gains(n).name];

%% calculate cortical noise parameter using data VR distance
for s = 1:4
    data_spks = squeeze(Spks_clean{dataCh}(:,5-s,1));
    time_end = 2971;
    data_FR(s) = mean(mean(cellfun(@(x) sum(x >= 0 & x < time_end/1000),data_spks)))/time_end*1000;
    spkTime = cell(size(data_spks));
    for ii = 1:length(data_spks)
        spkTime{ii} = data_spks{ii}';
        spkTime{ii}(spkTime{ii} < 0 | spkTime{ii} > time_end/1000) = []; 
    end
    spkTime = spkTime(~cellfun('isempty',spkTime));
    % calculate distance matrix & performance
    distMat = calcvr(spkTime, data_tau/1000);
    VR(s) = mean(distMat(distMat ~= 0),'all');
end


%% load STRF spikes

ICdirPath = [ICdir filesep];
ICstruc = dir([ICdirPath '*.mat']);

if isempty(ICstruc), error('empty data directory'); end

gsyn_range = 0.03:0.03:0.21;

for cases = 1:4
    
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
varies(end).range = 1.5;

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

[data,data_FR,DirPart] = mouseNetwork_initialize(varies,ICstruc,ICdirPath,Spks_clean,...
    Spks_masked,dataCh,data_tau,plot_distances,plot_rasters);

%% performance grids

if plot_grids
    makeGrids(data,varies,DirPart,nvaried,data_perf,data_FR);
end

temp = {data.name};
temp(cellfun('isempty',temp)) = {'empty'}; %label empty content

targetIdx = find(contains(temp,'m0') & ~strcmp(temp,'s0m0.mat'));
%maskerIdx = find(contains(temp,'s0') & ~strcmp(temp,'s0m0.mat'));
%mixedIdx = find(~contains(temp,'m0') & ~contains(temp,'s0') & ~contains(temp,'empty'));
perf = []; fr = [];
for i = 1:length(targetIdx)
    perf(:,5-i) = data(targetIdx(i)).perf.C;
    fr(:,5-i) = data(targetIdx(i)).fr.C;
end

[~,MSE_clean] = calcModelPerf(perf,data_perf(1:4));

loss = MSE_clean(:,1);
frdiffs = abs(mean(fr,2) - mean(data_FR(data_FR ~= 0)));
within_thresh = frdiffs <= 5;

% no wts are within FR thresh, just get the absolute minimum
if sum(within_thresh) == 0
   within_thresh = frdiffs == min(frdiffs);
end

[minloss(cases),i] = min(loss(within_thresh));
temp = find(within_thresh);
best_iteration = temp(i);

best_wt(cases) = gsyn_range(best_iteration);
best_perf(cases,:) = perf(best_iteration,:);

end

% determine best of first channel cases

[onechan_loss,best_channel] = min(minloss);
best_gSYN = best_wt(best_channel);

onechan_perf = best_perf(best_channel,:);

%% second channel cases
minloss = [];
best_wt = [];
twochan_perf = [];
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
varies(end).range = 1.5;

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

[data,data_FR,DirPart] = mouseNetwork_initialize(varies,ICstruc,ICdirPath,Spks_clean,...
    Spks_masked,dataCh,data_tau,plot_distances,plot_rasters);

%% performance grids

if plot_grids
    makeGrids(data,varies,DirPart,nvaried,data_perf,data_FR);
end

temp = {data.name};
temp(cellfun('isempty',temp)) = {'empty'}; %label empty content

targetIdx = find(contains(temp,'m0') & ~strcmp(temp,'s0m0.mat'));
%maskerIdx = find(contains(temp,'s0') & ~strcmp(temp,'s0m0.mat'));
%mixedIdx = find(~contains(temp,'m0') & ~contains(temp,'s0') & ~contains(temp,'empty'));
perf = []; fr = [];
for i = 1:length(targetIdx)
    perf(:,5-i) = data(targetIdx(i)).perf.C;
    fr(:,5-i) = data(targetIdx(i)).fr.C;
end

[~,MSE_clean] = calcModelPerf(perf,data_perf(1:4));

loss = MSE_clean(:,1);
frdiffs = abs(mean(fr,2) - mean(data_FR(data_FR ~= 0)));
within_thresh = frdiffs <= 5;

% no wts are within FR thresh, just get the absolute minimum
if sum(within_thresh) == 0
   within_thresh = frdiffs == min(frdiffs);
end

[minloss(cases),i] = min(loss(within_thresh));
temp = find(within_thresh);
best_iteration = temp(i);

[A,B] = meshgrid(best_gsyn_range,second_gsyn_range);
c = cat(2,A',B');
all_wt_combos = reshape(c,[],2);

best_two_wts(cases,:) = all_wt_combos(best_iteration,:);
twochan_perf(cases,:) = perf(best_iteration,:);

end

minloss(minloss == 0) = [];
best_two_wts(best_two_wts == 0) = [];
twochan_perf(twochan_perf == 0) = [];

% determine if adding second channel improves one-channel case

decreased_loss = (minloss - onechan_loss) < 0;
for cases = 1:2
   sig_perf(cases) = ttest(onechan_perf,twochan_perf(cases,:),'Tail','left'); 
end

if sum(sig_perf) == 0 && sum(decreased_loss) == 0
    disp('One-channel case is best');
else
    disp(['Two-channel case ' num2str(find(decreased_loss)) ' outperformed single channel']);
end

best_two_wts(decreased_loss);