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

%% define parameters for optimization

nSims = 1; % number of simulations to run for average MSE

ipsi = [0:0.05:0.2];
contra = [0.3:0.05:0.6];
[w1,w2] = meshgrid(ipsi,contra);
c=cat(2,w1',w2');
d=reshape(c,[],2);

weights(1,:) = d(:,1); %add this as input to mouse_network
weights(2,:) = zeros(size(d,1),1);
weights(3,:) = zeros(size(d,1),1);
weights(4,:) = d(:,2);

%%%%%%%% end of user inputs

inputChans = find(sum(weights,2) ~= 0);

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

for it = 1:size(weights,2)

disp(['Starting iteration #', num2str(it)]);
    
%% varied parameters
varies = struct;

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
% xrNetcon = zeros(nCells);

%srNetcon = diag(ones(1,nCells));
rcNetcon = weights(:,it);

% netCons.xrNetcon = xrNetcon;
% netCons.srNetcon = srNetcon;
netCons.rcNetcon = rcNetcon;

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
    
%     for ns = 1:nSims
%     [temp_perf(ns), temp_fr(ns), data(z).annot,~, data(z).VR] = ...
%         mouse_network(study_dir,time_end,varies,netCons,plot_rasters,...
%         plot_distances,data_spks,data_tau);
%     end
%     
%     % average across multiple simulations if needed
%     for vv = 1:nvaried
%         if nSims > 1
%             inds = 1:nvaried:nvaried*nSims + (vv-1);
%         else
%             inds = vv;
%         end
%         temp = [temp_perf.R]; data(z).perf.R(:,vv) = mean(temp(:,inds),2);
%         temp = [temp_perf.C]; data(z).perf.C(vv) = mean(temp(inds));
%         temp = [temp_fr.R]; data(z).fr.R(:,vv) = mean(temp(:,inds),2);
%         temp = [temp_fr.C]; data(z).fr.C(vv) = mean(temp(inds));
%     end

    [data(z).perf, data(z).fr, data(z).annot,~, data(z).VR] = ...
        mouse_network(study_dir,time_end,varies,netCons,plot_rasters,...
        plot_distances,data_spks,data_tau);

    data(z).name = ICstruc(z).name;
end
save([pwd filesep 'run' filesep 'grid-search' filesep...
    datetime filesep 'summary_results_iteration' num2str(it) '.mat'],'data','varies')
close(h);

%% performance grids

if plot_grids
    makeGrids(data,varies);
end

temp = {data.name};
temp(cellfun('isempty',temp)) = {'empty'}; %label empty content

targetIdx = find(contains(temp,'m0') & ~strcmp(temp,'s0m0.mat'));
% maskerIdx = find(contains(temp,'s0') & ~strcmp(temp,'s0m0.mat'));
% mixedIdx = find(~contains(temp,'m0') & ~contains(temp,'s0') & ~contains(temp,'empty'));

if ~isempty(targetIdx)
    for i = 1:length(targetIdx)
        perf(i) = data(targetIdx(i)).perf.C(vv);
        fr(i) = data(targetIdx(i)).fr.C(vv);
    end
end

[~,MSE_clean] = calcModelPerf(perf',data_perf(1:4));

perfHistory(it,:) = perf;
frHistory(it) = mean(fr);
loss(it) = MSE_clean(1);
within_thresh(it) = abs(mean(fr) - mean(data_FR(data_FR ~= 0))) <= 5;

[temp_model,temp_exp] = findmeanVR(data,targetIdx);
VRdiff(it) = mean((temp_model-temp_exp).^2);

end

%% need to manually code this when you change the parameter space for grid search

names = {'RC_{ipsi}','RC_{G}','RC_{U}','RC_{contra}'};

temp = weights(inputChans,:);
% check how many dimensions in parameter space is used
dims = 0;
for r = 1:length(inputChans)
    dims = dims + double(length(unique(temp(r,:))) ~= 1);
    if length(unique(temp(r,:))) == 1
        inputChans(r) = [];
        temp(r,:) = [];
    end
end

% plot loss vs param
figure;
if dims < 2
    [temp,I] = sort(temp);
    temp = temp;
    plot(temp,loss(I),'linewidth',1);
    xlabel(names{inputChans}); ylabel('Loss');
    hold on;
    scatter(temp(within_thresh),loss(within_thresh),80,[1 0 0],'filled');
    yyaxis right
    plot(temp,frHistory(I),'--k','linewidth',1);
    scatter(temp(within_thresh),frHistory(I),80,'filled');
    ylabel('Mean clean FR');
    title('Weights vs. Loss and Mean FR');
    plot(temp,ones(size(temp))*mean(data_FR(targetIdx)),'linewidth',1);
    legend('Loss','Loss within FR threshold','Model FR','Data FR','location','best');
    xlim([temp(1) temp(end)]);
saveas(gcf,[filesep DirPart filesep 'MSE and FR grid search.png'])

elseif dims == 2
    sz1 = length(unique(temp(1,:)));
    sz2 = length(unique(temp(2,:)));
    
    X = reshape(temp(1,:),sz1,sz2);
    Y = reshape(temp(2,:),sz1,sz2);
    Z = reshape(loss,sz1,sz2);
    surf(X,Y,Z);
    
    hold on;
    scatter3(X(within_thresh),Y(within_thresh),Z(within_thresh),80,[1 0 0],'filled')
    xlabel(names{inputChans(1)}); ylabel(names{inputChans(2)}); zlabel('Loss');
    zlim([0 200]);
    title('Weights vs. Loss');
    saveas(gcf,[filesep DirPart filesep 'MSE grid search 2d space.fig'])
    legend('Model','Models within threshold','location','best')
end

save([filesep DirPart filesep 'grid_search_results.mat'],'loss','weights',...
    'perfHistory','frHistory','within_thresh','VRdiff');

% plot FR vs param for 2d parameter space
figure;
if dims == 2
    sz1 = length(unique(temp(1,:)));
    sz2 = length(unique(temp(2,:)));
    
    X = reshape(temp(1,:),sz1,sz2);
    Y = reshape(temp(2,:),sz1,sz2);
    Z = reshape(frHistory,sz1,sz2);
    surf(X,Y,Z);
    hold on;
    scatter3(X(within_thresh),Y(within_thresh),Z(within_thresh),80,[1 0 0],'filled')
    %zlim([35.5-6 35.5+6]);
    title('Weights vs. Mean FR');
    xlabel(names{inputChans(1)}); ylabel(names{inputChans(2)}); zlabel('Mean clean FR');
    legend('Model','Models within threshold','location','best')

    saveas(gcf,[filesep DirPart filesep 'MSE FR 2d space.fig'])
end

[minloss,i] = min(loss(within_thresh));
good = find(within_thresh);
disp(['Weight with min. loss and within FR threshold: '...
    mat2str(weights(:,good(i))) ', loss = ' num2str(minloss),...
    ', iteration ' num2str(good(i))]);

load([filesep DirPart filesep 'summary_results_iteration' num2str(good(i)) '.mat')
makeGrids(data,varies);

%% plot weights vs VR MSE 
figure;
if dims < 2
    [temp,I] = sort(temp);
    temp = temp;
    plot(temp,VRdiff(I),'linewidth',1);
    xlabel(names{inputChans}); ylabel('Sq. VR difference');
    hold on;
    scatter(temp(within_thresh),VR(within_thresh),80,[1 0 0],'filled');
    title('Weights vs. Mean VR^{2} Difference');
    legend('Model','Models within threshold','location','best');
    xlim([temp(1) temp(end)]);
    saveas(gcf,[filesep DirPart filesep 'MSE and VR grid search.png'])
elseif dims == 2
    sz1 = length(unique(temp(1,:)));
    sz2 = length(unique(temp(2,:)));
    
    X = reshape(temp(1,:),sz1,sz2);
    Y = reshape(temp(2,:),sz1,sz2);
    Z = reshape(VRdiff,sz1,sz2);
    surf(X,Y,Z);
    
    hold on;
    scatter3(X(within_thresh),Y(within_thresh),Z(within_thresh),80,[1 0 0],'filled')
    xlabel(names{inputChans(1)}); ylabel(names{inputChans(2)}); zlabel('Sq. VR difference');
    title('Weights vs. Mean VR^{2} Difference');
    %zlim([0 200]);
    saveas(gcf,[filesep DirPart filesep 'MSE and VR grid search 2d space.fig'])
    legend('Model','Models within threshold','location','best')
end

