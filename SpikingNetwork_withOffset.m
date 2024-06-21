%% Initialize
profile on;
warning('off','all');
warning

% change current directory to folder where this script is stored
mfileinfo = mfilename('fullpath');
mfiledir = strsplit(mfileinfo,filesep);
% cd(fullfile(mfiledir{1:end-1}));

dynasimPath = '../DynaSim';

addpath('mechs');
addpath('resampled-stimuli');
addpath(genpath('ICSimStim'));
addpath('genlib');
addpath(genpath(dynasimPath));
addpath('cSPIKE'); InitializecSPIKE;
addpath('plotting');
addpath('subfunctions');


%% Make ICfiles.mat if it's not in your directory

if ~isfile('ICfiles.mat'), makeICfiles; end

%% user inputs
dt = 0.5; %ms, should be a multiple of 0.1 ms

% study_dir: folder under 'run' where m files and input spikes for simulations are written and saved
study_dir = fullfile(pwd,'run','4-channel-PV-inputs');

if exist(study_dir, 'dir'), msg = rmdir(study_dir, 's'); end
mkdir(fullfile(study_dir, 'solve'));

addpath(fullfile(study_dir, 'solve'));


% expName: folder under 'simData' where results are saved
expName = '12-15-23 very phasic, higher PV potential';
simDataDir = [pwd filesep 'simData' filesep expName];
if ~exist(simDataDir,'dir'), mkdir(simDataDir); end

%% Run .m file to generate options and varies structs for simulations
addpath('params-4-channel');
% params_4channel_cleanonly;
params_4channel_cleanonly_PVinputs;

% addpath('params-3-channel');
% params_3channel;

% To re-create figures in paper, look at 'params' directory
% (other 'params' m files are for masked configs)
% addpath('params');

%% Creating TuningCurves

% spatial tuning at inputs
[azi,spatialCurves,chanLabels,bestLocs] = genSpatiallyTunedChans(options.nCells);


%% Handle Netcons

NetconHandler;

%% load input stimuli (targets and maskers) from ICSimStim
load('default_STRF_with_offset_200k.mat');

% firing rates were generated with sampling rate of 10000 Hz to match old
% simulation time step, downsample if dt's don't match
if dt ~= 0.1
    dsamp_fac = dt/0.1;
    for m = 1:10
        fr_masker{m} = downsample(fr_masker{m},dsamp_fac);
    end
    for t = 1:2
        fr_target_on{t} = downsample(fr_target_on{t},dsamp_fac);
        fr_target_off{t} = downsample(fr_target_off{t},dsamp_fac);
    end
end

% edit strfGain if you want to rescale firing rate at inputs
newStrfGain = strfGain;

%% create input spikes from STRFs

padToTime = 3500; % [ms]

% ICfiles.mat contains names of spatial grid configs: s[targetloc]m[maskerloc]
% See 'config_idx_reference.JPG' for indexes

% options.locNum is defined in params .m file
% if it's empty, default to running all 24 configs (including masker-only
% trials)

if ~isempty(options.locNum)
    subz = options.locNum;
else
    subz = 1:24;
end

% concatenate spike-time matrices, save to study_dir
prepInputData;

%% run simulation

options.strfGain = newStrfGain;
options.dt = dt;

if isempty(options.locNum), options.time_end = size(spks,1)*dt; % [ms];
else, options.time_end = padToTime*numel(options.locNum); end

%In theory the thing below should have to go 4 times
%The question is, how do we get snn_out to reflect a certain target
%direction

%[snn_out,s] = columnNetwork_paper(study_dir,varies,options,netcons);
%[snn_out,s] = columnNetwork_simpler(study_dir,varies,options,netcons);
[snn_out,s] = columnNetwork_simpler_onoff(study_dir,varies,options,netcons);
%[snn_out,s] = columnNetwork_alternative(study_dir,varies,options,netcons);
%% post-process for performance and firing results

postProcessSims;
save('current_run_data.mat');

profile off;
profile viewer;

%4x4   * 1x4   * 4x4
