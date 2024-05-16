%% Initialize

% change current directory to folder where this script is stored
mfileinfo = mfilename('fullpath');
mfiledir = strsplit(mfileinfo,filesep);
cd(fullfile(mfiledir{1:end-1}));

dynasimPath = '../DynaSim';

addpath('mechs');
addpath('resampled-stimuli');
addpath(genpath('ICSimStim'));
addpath('genlib');
addpath(genpath(dynasimPath));
addpath('cSPIKE'); InitializecSPIKE;
addpath('plotting');
addpath('subfunctions');

%% user inputs
dt = 0.5; %ms

% study_dir: folder under 'run' where m files and input spikes for simulations are written and saved
study_dir = fullfile(pwd,'run','3-channel');

if exist(study_dir, 'dir'), msg = rmdir(study_dir, 's'); end
mkdir(fullfile(study_dir, 'solve'));

% expName: folder under 'simData' where results are saved
expName = '12-06-23 blah'
simDataDir = [pwd filesep 'simData' filesep expName];
if ~exist(simDataDir,'dir'), mkdir(simDataDir); end

%% Run .m file to generate options and varies structs for simulations
%% Options struct

options = struct;
options.nCells = 3; % <----- 
options.opto = 0;

options.mex_flag = 0;
options.parfor_flag = 0;
options.plotRasters = 0;

% locNum should be empty for full grids

% (see config_idx_reference.JPG for location numbers)
options.locNum = 15;
options.SpatialAttention = 0;

% use a separate struct for connectivity matrices (netcons) between populations
netcons = struct; % row = source, column = target
netcons.XRnetcon = zeros(options.nCells,options.nCells);
netcons.XRnetcon([2 2],[1 3]) = 1;

netcons.RCnetcon = ones(options.nCells,1);

%% define network parameters
clear varies

if options.opto, nSims = 5; else, nSims = 1; end

trialInds = repmat(1:20,nSims,1);

% % % DO NOT CHANGE THIS % % %
varies(1).conxn = 'On->On';
varies(1).param = 'trial';
varies(1).range =  1;
% % % DO NOT CHANGE THIS % % %

%% create spatially-tuned channels based on options.nCells

[azi,spatialCurves,chanLabels] = genSpatiallyTunedChans(options.nCells);

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

% edit this if you want to rescale firing rate at inputs
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

%% solver params
solverType = 'euler';
dt = options.dt; %ms
time_end = options.time_end;

%% neuron populations
% tonic = bias = cells spontaneous firing

nCells = options.nCells;

s = struct();

% onset column

s.populations(1).name = 'On';
s.populations(1).equations = 'noconLIF';
s.populations(1).size = nCells;
s.populations(1).parameters = {'t_ref',1};

s.populations(end+1).name='R1On';
s.populations(end).equations = 'noconLIF';
s.populations(end).size = nCells;
s.populations(end).parameters = {'g_inc',0.0003,'tau_ad',100,'t_ref',1};

s.populations(end+1).name='S1On';
s.populations(end).equations = 'noconLIF';
s.populations(end).size = nCells;
s.populations(end).parameters = {'g_L',1/100,'E_L',-57,'V_reset',-52,'t_ref',0.5};

%% connections

% ms
EE_rise = 0.7;  EE_fall = 1.5;   % E->E
IE_rise = 1;    IE_fall = 4.5;   % PV->E
EI_rise = 0.55; EI_fall = 1;     % E->PV
XE_rise = 2.5;  XE_fall = 6;     % SOM->E

% note: all rise times must be greater than dt
if any([EE_rise,IE_rise,EI_rise,XE_rise,EE_fall,IE_fall,EI_fall,XE_fall]) <= dt
    error('PSC time constants must be greater than simulation timestep.')
end

s.connections(1).direction='On->On';
s.connections(1).mechanism_list={'IC'};
s.connections(1).parameters={'g_postIC',0.265,'label','on','trial',1,'locNum',options.locNum,'netcon',eye(nCells,nCells),'t_ref',1,'t_ref_rel',1,'rec',2};

% excitatory inputs
s.connections(end+1).direction='On->R1On';
s.connections(end).mechanism_list={'PSC'};
s.connections(end).parameters={'gSYN',0,'tauR',EE_rise,'tauD',EE_fall,'fP',0.1,'tauP',30};

s.connections(end+1).direction='On->S1On';
s.connections(end).mechanism_list={'PSC'};
s.connections(end).parameters={'gSYN',0.03,'tauR',EI_rise,'tauD',EI_fall,'fP',0.2,'tauP',80};

s.connections(end+1).direction='S1On->R1On';
s.connections(end).mechanism_list={'PSC'};
s.connections(end).parameters={'gSYN',0.03,'tauR',IE_rise,'tauD',IE_fall,'ESYN',-80,'fP',0.5,'tauP',120}; 

%% simulate

% vary params
vary = cell(length(varies),3);
for i = 1:length(varies)
    vary{i,1} = varies(i).conxn;
    vary{i,2} = varies(i).param;
    vary{i,3} = varies(i).range;
end

tic;

snn_out = dsSimulate(s,'tspan',[dt time_end], 'solver',solverType, 'dt',dt,...
  'downsample_factor',1, 'save_data_flag',0, 'save_results_flag',1,...
  'study_dir',study_dir, 'vary',vary, 'debug_flag', 1, 'verbose_flag',0,...
  'mex_flag',options.mex_flag);

toc;

