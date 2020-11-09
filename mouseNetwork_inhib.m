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

fitFlag = 0;    % if only want to generate clean and co-located spots, =1
plot_rasters = 1;

%%%%%%%%

subject = [extractBefore(perf_file,'_performance') '-Ch' num2str(dataCh)];
folder = ['Data-fitting' filesep subject];

load(fullfile(folder,'default_parameters.mat'));

% manually choose strf instead of using the saved ICdir!
ICdir = uigetdir('STRFs/BW_0.009_BTM_3.8_t0_0.1_phase0.4985/');
if ICdir == 0
   load(fullfile(folder,'default_parameters.mat'),'ICdir');
end
ICdirPath = [ICdir filesep];

%% Initiate varies

varies = struct;

varies(1).conxn = '(IC->IC)';
varies(1).param = 'trial';
varies(1).range = 1:40;%[1:10 21:30]; %

varies(end+1).conxn = 'C';
varies(end).param = 'noise';
varies(end).range = 1.5;

Cnoise2 = 1; % additional noise, so colocated noise = Cnoise2 + varies(2).range

varies(end+1).conxn = 'R->C';
varies(end).param = 'gSYN1';
varies(end).range = 0;

varies(end+1).conxn = 'R->C';
varies(end).param = 'gSYN2';
varies(end).range = 0.14;

varies(end+1).conxn = 'R->C';
varies(end).param = 'gSYN3';
varies(end).range = 0;

varies(end+1).conxn = 'R->C';
varies(end).param = 'gSYN4';
varies(end).range = 0;

% row: origin (does the inhibition)
% column: destination (to be inhibited)
xrNetcons = zeros(4);
xrNetcons(3,2) = 1;

varies(end+1).conxn = 'X->R';
varies(end).param = 'gSYN';
varies(end).range = 0.05;

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

[simdata,DirPart] = mouseNetwork_inhib_initialize(varies,xrNetcons,Cnoise2,ICdirPath,...
    spks_file,dataCh,plot_rasters,folder,detail,subz,[]);

%% performance grids for 4D search

temp = {simdata.name};
temp(cellfun('isempty',temp)) = {'empty'}; %label empty content

targetIdx = find(contains(temp,'m0') & ~strcmp(temp,'s0m0.mat'));
colocIdx = find((cellfun(@(x) x(2),temp) == cellfun(@(x) x(4),temp)));

% makeParallelPlot_adapt(simdata,data_perf(1:4),data_perf(5:end),DirPart);

nvaried = size(simdata(targetIdx(1)).annot,1);

perf_clean = [];
% clean performance and FR
for i = 1:length(targetIdx)
    perf_clean(:,i) = simdata(targetIdx(i)).perf.C;
end

perf_masked = [];
% mixed performance and FR
for i = 1:length(colocIdx)
    perf_masked(:,i) = simdata(colocIdx(i)).perf.C;
end

[~,MSE_clean_perf] = calcModelLoss(perf_clean,data_perf(1:4)');
[~,MSE_masked_perf] = calcModelLoss(perf_masked,data_perf([5:5:end])');
loss = MSE_clean_perf(:,1) + MSE_masked_perf(:,1);

% plot PSTHs below rasters
[y1,~] = audioread('stimuli-fixed-V2/200k_target1.wav');
[y2,fs] = audioread('stimuli-fixed-V2/200k_target2.wav');
mask = audioread('stimuli-fixed-V2/200k_masker1.wav');

y1(1:round(0.250*fs)) = []; y2(1:round(0.250*fs)) = []; mask(1:round(0.250*fs)) = [];

if plot_rasters
    
    rasters = dir([DirPart filesep '*.fig']);
    temp = {rasters.name};
    rasters(strcmp(cellfun(@(x) x(4),temp,'uniformoutput',false),'0')) = [];
        
    % plot PSTHs of clean trials over mixed trials
    for i = 1:length(targetIdx)
        
        % find all varied
        varied = find(contains(temp,'m0') & (cellfun(@(x) x(2),temp) == num2str(i)));
        
        for n = varied
            openfig([DirPart filesep temp{n}],'invisible');
            
            % plot stimuli
            subplot('position',[0.08 0.8 0.22 0.15]);hold off;
            plot((0:length(y1)-1)/fs,y1); hold on;
            plot((0:length(y2)-1)/fs,y2 + 2);
            plot([0 (length(y1)-1)/fs],[1 1],'k');
            set(gca,'xtick',[]);
            set(gca,'ytick',[]);
            
            subplot('position',[0.38 0.8 0.22 0.15]);hold off;
            plot((0:length(y1)-1)/fs,y1); hold on;
            plot((0:length(y2)-1)/fs,y2 + 2);
            plot([0 (length(y1)-1)/fs],[1 1],'k');
            set(gca,'xtick',[]);
            set(gca,'ytick',[]);
            
            saveas(gcf,[DirPart filesep temp{n}(1:end-4) '.png']);
            close;
        end
        
        % find mixed spots with the same target location
        for m = 1:4
%             
%             cleanPSTH = simdata(targetIdx(i)).PSTH;
%             mixedPSTH = simdata(targetIdx(i)+m).PSTH;
%             t_vec1 = 0:0.02:(size(cleanPSTH,2)-1)*0.02;
%             t_vec2 = 0:0.02:(size(mixedPSTH,2)-1)*0.02;   
            
            % find all rasters for given spatial grid spot
            varied = find((cellfun(@(x) x(4),temp) == num2str(m) & (cellfun(@(x) x(2),temp) == num2str(i))));
            
            for n = varied

            openfig([DirPart filesep temp{n}],'invisible');
            
            subplot('position',[0.08 0.8 0.22 0.15]); hold off;
            plot((0:length(y1)-1)/fs,y1 + mask); hold on;
            plot((0:length(y2)-1)/fs,y2 + mask + 2);
            plot([0 (length(y1)-1)/fs],[1 1],'k');
            set(gca,'xtick',[]);
            set(gca,'ytick',[]);
            
            subplot('position',[0.38 0.8 0.22 0.15]); hold off;
            plot((0:length(y1)-1)/fs,y1 + mask); hold on; 
            plot((0:length(y2)-1)/fs,y2 + mask + 2);
            plot([0 (length(y1)-1)/fs],[1 1],'k');
            set(gca,'xtick',[]);
            set(gca,'ytick',[]);
            
            % plot clean and masked PSTHs on top of each other in separate
            % subplot
            
%             subplot('position',[0.68 0.05 0.22 0.25]);
%             hold off;
%             mixed = plot(t_vec2,2*50*mixedPSTH(1,:)/size(varies(1).range,2) + 150); hold on;
%             plot(t_vec2,2*50*mixedPSTH(2,:)/size(varies(1).range,2));
%             
%             plot([0 3],[150 150],'k');
%             
%             clean = plot(t_vec1,(2*50*cleanPSTH(1,:))/size(varies(1).range,2) + 150);
%             plot(t_vec1,2*50*cleanPSTH(2,:)/size(varies(1).range,2));
%             legend(mixed,clean,'Mixed','Clean');
%             ylim([0 300]);

            savefig([DirPart filesep temp{n}]);
            saveas(gcf,[DirPart filesep temp{n}(1:end-4) '.png']);
            close;
            
            end
            
        end
    end
end



% make spatial grid
if length(subz) ~= 20
    makeGrids_fitting(simdata,DirPart,data_perf,data_FR,data_FR_masked(1:5:end));
else
    makeGrids(simdata,DirPart,data_perf,data_FR,data_FR_masked(1:5:end),loss)
end

%[] = find(xrNetcons == 1);
