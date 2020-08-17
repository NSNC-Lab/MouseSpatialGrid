clearvars;
close all;

addpath('mechs')
addpath('dependencies')
addpath('eval_scripts')
addpath('genlib')
addpath(genpath('dynasim'))

%%%%%%%% start of user inputs

plot_rasters = 0;   % plot rasters

data_spks_file = '9-21-2016_0dB_removed_trials_cleaned(-1,4).mat';
data_perf_file = '9-21-2016_0dB_removed_trials_performance.mat';
dataCh = 31;

%%%%%%%% end of user inputs

subject = [extractBefore(data_perf_file,'_performance') '-Ch' num2str(dataCh)];
folder = 'Data-fitting';

%% load experimental data to optimize model to

nCells = 4; %synchronise this variable with mouse_network

load(data_spks_file,'Spks_clean','Spks_masked','n_cum_clean','n_cum_masked');
load(data_perf_file,'Max','max_masked','opt_tau','opt_tau_masked');

data_perf = Max{dataCh}(:); % clean performance
for t = 1:4
    data_perf(end+1:end+4) = flipud(max_masked{dataCh}(:,t)); % mixed performance
end
% data_perf(1:4) - clean perf contra->ipsi
% data_perf(5:end) - target(m m m m) target (m m m m) etc.
        % target - contra->ipsi
        % masker - contra->ipsi

% use tau of best clean location
[~,best_loc] = max(Max{dataCh});
data_tau = opt_tau{dataCh}(best_loc);

%% calculate cortical noise based on spontaneous FR

% data_spks arranged contra->ipsi
FR_r0 = zeros(1,4);
data_FR = zeros(1,4);
for s = 1:4
    data_spks = squeeze(Spks_clean{dataCh}(:,s,:));
    FR_r0(s) = mean(cellfun(@(x) sum(x < 0 | x >= 3),data_spks),'all')/2;
    data_FR(s) = mean(cellfun(@(x) sum(x >= 0 & x < 3),data_spks),'all')/3;
end

% data_FR_masked arranged: 4 inds for target at -90, 4 inds at 0, etc.
% similar arrangement as data_FR for clean
for s = 1:16
    ti = ceil(s/4);
    mi = mod(s,4);
    if mi == 0
       mi = 4; 
    end
   data_spks = squeeze(Spks_masked{dataCh}(:,ti,mi,:));
   data_FR_masked(s) = mean(cellfun(@(x) sum(x >= 0 & x < 3),data_spks),'all')/3;
end
load('Cnoise_vs_FR0.mat','fit');

Cnoise = (mean(FR_r0)-fit(1))/fit(2);

%% calculate best STRF gain data based on best channel

% extract all STRF gains
gains = dir(fullfile(pwd,'STRFs/BW_0.009 BTM_3.8 t0_0.1 phase0.4985'));
gains(~contains({gains.name},'s30')) = [];

onset = 0.3; % s

STRF_peak = zeros(length(gains),2);
STRF_avg = zeros(length(gains),2);
for g = 1:length(gains)
    ICdir = [gains(g).folder filesep gains(g).name];
    ICdirPath = [ICdir filesep];
    ICstruc = dir([ICdirPath '*.mat']);
    
    subz(1) = find(strcmp({ICstruc.name},['s' num2str(best_loc) 'm0.mat'])); % sXm0 (target only) case
    subz(2) = find(strcmp({ICstruc.name},['s' num2str(best_loc) 'm' num2str(best_loc) '.mat'])); % sXmX (colocated) case
    i = 1;
    for z = subz
        load([ICdirPath ICstruc(z).name],'fr','avgSpkRate');
        temp = [];
        % fr in STRFs is stored like [trials x (target 1+target2)]
        for t = 1:20
            temp(t,:) = fr{t,best_loc};
            temp(t+20,:) = fr{t,best_loc+4};
        end
        temp(isnan(temp)) = 0;
        temp = mean(temp);  % find average FR across trials
        
        % change FR trace from per second to per every 20 ms
        t_vec = 1:20:length(temp);
        for t = 1:length(t_vec)-1
            STRF_PSTH{g,i}(t) = mean(temp(t_vec(t):t_vec(t+1)-1));
        end
        % find peak FR during user-defined onset period
        STRF_peak(g,i) = max(STRF_PSTH{g,i}(t_vec >= 0 & t_vec <= onset*1000));
        
        % get avg FR during stimulus for calculating synaptic weights
        % for model
        STRF_avg(g,i) = avgSpkRate(best_loc);

        i = i + 1;
    end
end
% % find best gain closest to average clean trial FR
% [~,ind] = min(abs(mean(STRF_FR,2) - mean(data_FR)));

% fit STRF based on peak inital FR in data
t_vec = -1:1/50:4;
t_inds = t_vec >= 0 & t_vec < 3;

scale_data = 50/size(Spks_clean{dataCh},1);

% calculate average data PSTH across all clean trials in best channel
PSTH_clean_data = scale_data*mean([n_cum_clean{dataCh}{best_loc,1}(t_inds);n_cum_clean{dataCh}{best_loc,2}(t_inds)]);
clean_peak = max(PSTH_clean_data(1:(50*onset)+1));

PSTH_mixed_data = scale_data*mean([n_cum_masked{dataCh}{best_loc,best_loc,1}(t_inds);n_cum_masked{dataCh}{best_loc,best_loc,2}(t_inds)]);
mixed_peak = max(PSTH_mixed_data(1:(50*onset)+1));

% best fit based on clean data
[~,ind] = min(abs(STRF_peak(:,1) - clean_peak));

gain_ratio = clean_peak/mean(STRF_peak(ind,1),2);

% plot data vs STRF psths - clean location
figure;
subplot(2,1,1); plot(0:0.02:0.02*(length(PSTH_clean_data)-1),PSTH_clean_data)
title('PSTH of clean data at best location')
line([onset onset],[0 300],'color','r');
subplot(2,1,2); plot(0:0.02:0.02*(length(STRF_PSTH{ind,1})-1),STRF_PSTH{ind,1});
line([onset onset],[0 300],'color','r');
title('PSTH of STRF at best location')

% plot data vs STRF psths - masked location
figure;
subplot(2,1,1); plot(0:0.02:0.02*(length(PSTH_mixed_data)-1),PSTH_mixed_data)
title('PSTH of co-located data at best location')
line([onset onset],[0 300],'color','r');
subplot(2,1,2); plot(0:0.02:0.02*(length(STRF_PSTH{ind,2})-1),STRF_PSTH{ind,2});
line([onset onset],[0 300],'color','r');
title('PSTH of co-located STRF at best location')

ICdir = [gains(ind).folder filesep gains(ind).name];
ICdirPath = [ICdir filesep];
ICstruc = dir([ICdirPath '*.mat']);

% back of envelope calculation of upper limit of sum of synaptic
% conductances
gsyn_sum = gain_ratio*0.21*(mean(data_FR) - mean(FR_r0))./(STRF_avg(ind,1));

%% Adaptation conductance search

% for 2d adaptation search, use single-channel case
ranges = cell(1,4);
for c = 1:4
    if c == best_loc
        ranges{c} = gsyn_sum;
    else
        ranges{c} = 0; 
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

varies(end+1).conxn = 'C';
varies(end).param = 'G_inc';
varies(end).range = [0.001:0.001:0.005];
    
varies(end+1).conxn = 'C';
varies(end).param = 'tau_ad';
varies(end).range = [45:45:270];

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

adaptFlag = 1;
allFlag = 0;
restrict_vary_flag = 0;

[data,DirPart] = mouseNetwork_initialize(varies,ICstruc,ICdirPath,Spks_clean,Spks_masked,dataCh,plot_rasters,folder,subject,'-adaptation-search',best_loc,adaptFlag,allFlag,restrict_vary_flag);

%% compare model and data PSTHs

temp = {data.name};
temp(cellfun('isempty',temp)) = {'empty'}; %label empty content

targetIdx = find(contains(temp,'m0') & ~strcmp(temp,'s0m0.mat'));
mixedIdx = find(~contains(temp,'m0') & ~contains(temp,'s0') & ~contains(temp,'empty'));

nvaried = numel(data(mixedIdx).G_inc);

t_vec = -1:1/50:4;
startInd = find(t_vec == 0);
for vv = 1:nvaried

    endInd = startInd + size(data(targetIdx).PSTH,2) - 1;
    
    [temp,lags] = xcorr(data(targetIdx).PSTH(vv,:,1),...
        n_cum_clean{dataCh}{best_loc,1}(startInd:endInd),'normalized');
    fit(vv,1) = temp(lags == 0);
    
    [temp,lags] = xcorr(data(targetIdx).PSTH(vv,:,1),...
        n_cum_clean{dataCh}{best_loc,2}(startInd:endInd),'normalized');
    fit(vv,2) = temp(lags == 0);
    
    endInd = startInd + size(data(mixedIdx).PSTH,2) - 1;
    
    [temp,lags] = xcorr(data(mixedIdx).PSTH(vv,:,1),...
        n_cum_masked{dataCh}{best_loc,best_loc,1}(startInd:endInd),'normalized');
    fit(vv,3) = temp(lags == 0);
    
    [temp,lags] = xcorr(data(mixedIdx).PSTH(vv,:,2),...
        n_cum_masked{dataCh}{best_loc,best_loc,2}(startInd:endInd),'normalized');
    fit(vv,4) = temp(lags == 0);
end

% plot grid of adaptation fits

X = unique(data(mixedIdx).tau_ad);
Y = unique(data(mixedIdx).G_inc);

figure('visible','on');
imagesc(X,Y,reshape(mean(fit,2),numel(X),numel(Y))');
xlabel('Time constant (ms)');
ylabel('Adaptation strength');
xticks(X);
yticks(Y);
title('Adaptation search by PSTH fit');
h = colorbar;
ylabel(h,'Cross-correlation at 0 lag');

% find adaptation values within FR threshold

scale_model = (1/20)/3;
scale_data = (1/size(Spks_clean{dataCh},1))/3;

startInd = find(t_vec == 0);
for vv = 1:nvaried

    endInd = startInd + size(data(mixedIdx).PSTH,2) - 1;
    
    FR_threshold(vv,1) = abs(scale_model*sum(data(mixedIdx).PSTH(vv,:,1),2)-(scale_data*sum(n_cum_masked{dataCh}{best_loc,best_loc,1}(startInd:endInd)))) <= 5;
    FR_threshold(vv,2) = abs(scale_model*sum(data(mixedIdx).PSTH(vv,:,2),2)-(scale_data*sum(n_cum_masked{dataCh}{best_loc,best_loc,2}(startInd:endInd)))) <= 5;
end

FR_threshold = sum(FR_threshold,2);

% plot scatter of adaptation values w/in FR threshold on top of 2D grid
hold on;
scatter3(data(mixedIdx).tau_ad(FR_threshold == 1),data(mixedIdx).G_inc(FR_threshold == 1),...
    FR_threshold(FR_threshold == 1),'filled','r');

saveas(gcf,fullfile([filesep DirPart],'adaptation_search.tiff'));

temp_best = max(mean(fit(FR_threshold == 1,:),2));
ind_ad = find(mean(fit,2) == temp_best);
G_inc = data(mixedIdx).G_inc(ind_ad);
tau_ad = data(mixedIdx).tau_ad(ind_ad);

% plot best model PSTHs vs. data PSTHs
scale_model = 50/20;
scale_data = 50/size(Spks_clean{dataCh},1);
for t = 1:2
    
    figure('visible','on');
    mixedT = 0:0.02:0.02*(length(data(mixedIdx).PSTH(ind_ad,:,t))-1);
    
    endInd = startInd + size(data(mixedIdx).PSTH,2) - 1;
    
    subplot(2,1,1); plot(mixedT,data(mixedIdx).PSTH(ind,:,t)*scale_model);ylim([0 150]) % multiply by 20 due to
    xlabel('Time (s)'); ylabel('FR'); title(['Model PSTH, masked target ' num2str(t)]);
    subplot(2,1,2); plot(mixedT,n_cum_masked{dataCh}{best_loc,best_loc,t}(startInd:endInd)*scale_data);ylim([0 150])
    xlabel('Time (s)'); ylabel('FR'); title(['Experiment PSTH, masked target ' num2str(t)]);
    
    saveas(gcf,fullfile(filesep,DirPart,['PSTHs_target_' num2str(t) '.tiff']));
    %close;
end

save(fullfile([filesep DirPart],'fitting_results.mat'),'G_inc','tau_ad','Cnoise',...
    'ICdirPath','ICstruc','gsyn_sum','data_FR','data_FR_masked','data_perf');

%% 4D search
gsyn_range = 0:gsyn_sum/5:gsyn_sum;

ranges = cell(1,4);
for c = 1:4
    ranges{c} = gsyn_range; 
end

%% varied parameters

varies = struct;

varies(1).conxn = '(IC->IC)';
varies(1).param = 'trial';
varies(1).range = 1:40;

varies(end+1).conxn = 'C';
varies(end).param = 'noise';
varies(end).range = Cnoise;

varies(end+1).conxn = 'C';
varies(end).param = 'G_inc';
varies(end).range = G_inc;

varies(end+1).conxn = 'C';
varies(end).param = 'tau_ad';
varies(end).range = tau_ad;

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

% for now: last varies field has to be restriction
varies(end+1).conxn = 'R->C';
varies(end).param = {'gSYN1','gSYN2','gSYN3','gSYN4'};
varies(end).range = gsyn_sum;

nvaried = {varies(2:end).range};
nvaried = prod(cellfun(@length,nvaried));

allFlag = 0;
restrict_vary_flag = 1;
adaptFlag = 0;
[data,DirPart] = mouseNetwork_initialize(varies,ICstruc,ICdirPath,Spks_clean,Spks_masked,...
    dataCh,plot_rasters,folder,subject,'-4D-search',[],adaptFlag,allFlag,restrict_vary_flag);

%% performance grids for 4D search

temp = {data.name};
temp(cellfun('isempty',temp)) = {'empty'}; %label empty content

targetIdx = find(contains(temp,'m0') & ~strcmp(temp,'s0m0.mat'));
%maskerIdx = find(contains(temp,'s0') & ~strcmp(temp,'s0m0.mat'));
mixedIdx = find(~contains(temp,'m0') & ~contains(temp,'s0') & ~contains(temp,'empty'));
perf_clean = []; model_FR = [];

% clean performance and FR
for i = 1:length(targetIdx)
    perf_clean(:,i) = data(targetIdx(i)).perf.C;
    model_FR(:,i) = data(targetIdx(i)).fr.C;
end

perf_masked = []; model_FR_masked = [];
% mixed performance and FR
if ~isempty(mixedIdx)
for i = 1:length(mixedIdx)
    perf_masked(:,i) = data(mixedIdx(i)).perf.C;
    model_FR_masked(:,i) = data(mixedIdx(i)).fr.C;
end
end

[~,MSE_clean_perf] = calcModelLoss(perf_clean,data_perf(1:4)');
% calculate loss for FR as a percentage of model FR
[~,MSE_clean_FR] = calcModelLoss(100*(model_FR)./data_FR,100*data_FR./data_FR);

if ~isempty(mixedIdx)
    [~,MSE_masked_perf] = calcModelLoss(perf_masked,data_perf([5 10 15 20])');
    % calculate loss for FR as a percentage of model FR
    [~,MSE_masked_FR] = calcModelLoss(100*(model_FR_masked)./data_FR_masked([1 6 11 16]),100*data_FR_masked([1 6 11 16])./data_FR_masked([1 6 11 16]));
    loss = MSE_clean_perf(:,1) + MSE_clean_FR(:,1) + MSE_masked_perf(:,1) + MSE_masked_FR(:,1);
else
    loss = MSE_clean_perf(:,1) + MSE_clean_FR(:,1);
end

frdiffs = abs(mean(model_FR,2) - mean(data_FR(data_FR ~= 0)));
within_thresh = frdiffs <= 5;

% no wts are within FR thresh, just get the absolute minimum
if sum(within_thresh) == 0
   within_thresh = frdiffs == min(frdiffs);
end

% plotFRvsgSYN(data,within_thresh);


% find sets of parameters that are within threshold
[temp,inds] = sort(loss(within_thresh),'ascend');

if sum(within_thresh) >= 5
    best_iterations = zeros(1,5);
    makeParallelPlot(data,within_thresh,loss);
else
    best_iterations = zeros(1,sum(within_thresh));
end

for i = 1:length(best_iterations)
    best_iterations(i) = find(loss == temp(i));
end

% make grid of best iterations
makeGrids_bestIteration(data,DirPart,data_perf,data_FR,best_iterations,loss);

%% rerun dynasim and obtain rasters for best iteration

gsyn_strs = cellfun(@str2num,extractAfter({data(targetIdx(1)).annot{:,end}},'RC_{gSYN} = '),'UniformOutput',false);
best_gsyns = gsyn_strs{best_iterations(1)};

varies = struct;

varies(1).conxn = '(IC->IC)';
varies(1).param = 'trial';
varies(1).range = 1:40;

varies(end+1).conxn = 'C';
varies(end).param = 'noise';
varies(end).range = Cnoise;

varies(end+1).conxn = 'C';
varies(end).param = 'G_inc';
varies(end).range = G_inc;

varies(end+1).conxn = 'C';
varies(end).param = 'tau_ad';
varies(end).range = tau_ad;

varies(end+1).conxn = 'R->C';
varies(end).param = 'gSYN1';
varies(end).range = best_gsyns(1);

varies(end+1).conxn = 'R->C';
varies(end).param = 'gSYN2';
varies(end).range = best_gsyns(2);

varies(end+1).conxn = 'R->C';
varies(end).param = 'gSYN3';
varies(end).range = best_gsyns(3);

varies(end+1).conxn = 'R->C';
varies(end).param = 'gSYN4';
varies(end).range = best_gsyns(4);

allFlag = 1; % run all spots on grid
restrict_vary_flag = 0;
plot_rasters = 1;
adaptFlag = 0;

[data,DirPart] = mouseNetwork_initialize(varies,ICstruc,ICdirPath,Spks_clean,...
    Spks_masked,dataCh,plot_rasters,folder,subject,'-best-iteration',[],adaptFlag,allFlag,restrict_vary_flag);

%% make full spatial grid

temp = {data.name};
temp(cellfun('isempty',temp)) = {'empty'}; %label empty content

targetIdx = find(contains(temp,'m0') & ~strcmp(temp,'s0m0.mat'));
%maskerIdx = find(contains(temp,'s0') & ~strcmp(temp,'s0m0.mat'));
mixedIdx = find(~contains(temp,'m0') & ~contains(temp,'s0') & ~contains(temp,'empty'));
perf_clean = []; model_FR = [];

% clean performance and FR
for i = 1:length(targetIdx)
    perf_clean(:,i) = data(targetIdx(i)).perf.C;
    model_FR(:,i) = data(targetIdx(i)).fr.C;
end

perf_masked = []; model_FR_masked = [];
% mixed performance and FR
if ~isempty(mixedIdx)
for i = 1:length(mixedIdx)
    perf_masked(:,i) = data(mixedIdx(i)).perf.C;
    model_FR_masked(:,i) = data(mixedIdx(i)).fr.C;
end
end

[~,MSE_clean_perf] = calcModelLoss(perf_clean,data_perf(1:4)');
% calculate loss for FR as a percentage of model FR
[~,MSE_clean_FR] = calcModelLoss(100*(model_FR)./data_FR,100*data_FR./data_FR);

if ~isempty(mixedIdx)
    [~,MSE_masked_perf] = calcModelLoss(perf_masked,data_perf(5:end)');
    % calculate loss for FR as a percentage of model FR
    [~,MSE_masked_FR] = calcModelLoss(100*(model_FR_masked)./data_FR_masked,100*data_FR_masked./data_FR_masked);
    loss = MSE_clean_perf(:,1) + MSE_clean_FR(:,1) + MSE_masked_perf(:,1) + MSE_masked_FR(:,1);
else
    loss = MSE_clean_perf(:,1) + MSE_clean_FR(:,1);
end

makeGrids(data,DirPart,data_perf,data_FR,loss);
