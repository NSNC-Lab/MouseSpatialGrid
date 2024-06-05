function [simdata,s] = columnNetwork_simpler_onoff(study_dir,varies,options,netcons)
NetconHandler;
% Generates and simulates a network featuring columns of excitatory cells 
% that respond to onsets and offsets in auditory stimuli

% Features of circuit:

% Separate PV cells for onsets and offsets in stimuli; the PV cells will
% inhibit both onset and offset relay units in the same layer

% Output layer features single neuron ('C') where all spatially-tuned
% channels converge

% This is the closest combination of the AIM network and single-channel
% model to use

% INPUTS:
% study_dir: location of IC spike files + directory for log and data files
% varies: vary parameter, a structure. e.g.
%   varies(1).conxn = '(IC->IC)';
%   varies(1).param = 'trial';
%   varies(1).range = 1:20;
%   *first set of parameters should always be "trial" for IC->IC cxn
% options: struct with fields
%   nCells: number of channels in network
%   opto: simulating optogenetic suppression
%   SpatialAttention: =1 if simulating SpatialAttention
%   dt: timestep of trial (in ms)
%   time_end: length of simulations

% OUTPUTS:
% simdata - simulated voltages and spikes for all units
% s - model architecture (populations and connections)

% @Jio Nocon, BU 2023-12-02
% New version based on columnNetwork_V2 without L4 relays

%% Input check
if ~strcmp(varies(1).param,'trial')
    error('first set of varied params should be ''trial''')
end

%% solver params
solverType = 'euler';
dt = options.dt; %ms
time_end = options.time_end;

%% neuron populations
% tonic = bias = cells spontaneous firing

nCells = options.nCells;

s = struct();

XRnetcon = netcons.XRnetcon;
RCnetcon = netcons.RCnetcon;
PEnetcon = netcons.PEnetcon;

% onset column

s.populations(1).name = 'On';
s.populations(1).equations = 'noconLIF';
s.populations(1).size = nCells;
s.populations(1).parameters = {'t_ref',1};

s.populations(end+1).name = 'Off';
s.populations(end).equations = 'noconLIF';
s.populations(end).size = nCells;
s.populations(end).parameters = {'t_ref',1};

s.populations(end+1).name='ROn';
s.populations(end).equations = 'noconLIF';
s.populations(end).size = nCells;
s.populations(end).parameters = {'g_inc',0.0003,'tau_ad',100,'t_ref',1};

s.populations(end+1).name='ROff';
s.populations(end).equations = 'noconLIF';
s.populations(end).size = nCells;
s.populations(end).parameters = {'g_inc',0.0003,'tau_ad',100,'t_ref',1};

%s.populations(end+1).name='SOn';
%s.populations(end).equations = 'noconLIF';
%s.populations(end).size = nCells;
%s.populations(end).parameters = {'g_L',1/100,'E_L',-57,'V_reset',-52,'t_ref',0.5};

%s.populations(end+1).name='SOff';
%s.populations(end).equations = 'noconLIF';
%s.populations(end).size = nCells;
%s.populations(end).parameters = {'g_L',1/100,'E_L',-57,'V_reset',-52,'t_ref',0.5};

s.populations(end+1).name='SOnOff';
s.populations(end).equations = 'noconLIF';
s.populations(end).size = nCells;
s.populations(end).parameters = {'g_L',1/100,'E_L',-57,'V_reset',-52,'t_ref',0.5};

% TD unit
s.populations(end+1).name='TD';
s.populations(end).equations = 'noconLIF_currentOnly';
s.populations(end).size = nCells;
s.populations(end).parameters = {'Itonic',0.1,'numLocs',numel(options.locNum)};

% cross-channel X units, modeled as SOM units
s.populations(end+1).name='X';
s.populations(end).equations = 'noconLIF';
s.populations(end).size = nCells;
s.populations(end).parameters = {'g_L',1/275,'t_ref',1,'V_reset',-55};

% output unit
s.populations(end+1).name='C';
s.populations(end).equations = 'noconLIF';
s.populations(end).size = 1;
s.populations(end).parameters = {'g_inc',0.0003,'tau_ad',100,'t_ref',1};

%% connections

% ms
EE_rise = 0.7;  EE_fall = 1.5;   % E->E
IE_rise = 1;    IE_fall = 4.5;   % PV->E
EI_rise = 0.55; EI_fall = 1;     % E->PV
XE_rise = 2;    XE_fall = 8;    % SOM->E

% note: all rise times must be greater than dt
if any([EE_rise,IE_rise,EI_rise,XE_rise,EE_fall,IE_fall,EI_fall,XE_fall] <= dt)
    error('PSC time constants must be greater than simulation timestep.')
end

% % % Input layer % % %

s.connections(1).direction='On->On';
s.connections(1).mechanism_list={'IC'};
s.connections(1).parameters={'g_postIC',0.265,'label','on','trial',1,'locNum',options.locNum,'netcon',eye(nCells,nCells),'t_ref',1,'t_ref_rel',1,'rec',2};

s.connections(end+1).direction='Off->Off';
s.connections(end).mechanism_list={'IC'};
s.connections(end).parameters={'g_postIC',0.265,'label','off','trial',1,'locNum',options.locNum,'netcon',eye(nCells,nCells),'t_ref',1,'t_ref_rel',1,'rec',2};

% % % Intermediate layer % % %

% excitatory inputs
s.connections(end+1).direction='On->ROn';
s.connections(end).mechanism_list={'PSC'};
s.connections(end).parameters={'gSYN',0.02,'tauR',EE_rise,'tauD',EE_fall,'fP',0.1,'tauP',30,'netcon',eye(nCells,nCells)};

s.connections(end+1).direction='On->SOnOff';
s.connections(end).mechanism_list={'PSC'};
s.connections(end).parameters={'gSYN',0.03,'tauR',EI_rise,'tauD',EI_fall,'fP',0.2,'tauP',80,'netcon',eye(nCells,nCells)};

s.connections(end+1).direction='SOnOff->ROn';
s.connections(end).mechanism_list={'PSC'};
s.connections(end).parameters={'gSYN',0.03,'tauR',IE_rise,'tauD',IE_fall,'ESYN',-80,'fP',0.4,'tauP',100,'netcon',PEnetcon}; 

s.connections(end+1).direction='SOnOff->ROff';
s.connections(end).mechanism_list={'PSC'};
s.connections(end).parameters={'gSYN',0.03,'tauR',IE_rise,'tauD',IE_fall,'ESYN',-80,'fP',0.4,'tauP',100,'netcon',eye(nCells,nCells)}; 

% offset channels
s.connections(end+1).direction='Off->ROff';
s.connections(end).mechanism_list={'PSC'};
s.connections(end).parameters={'gSYN',0.02,'tauR',EE_rise,'tauD',EE_fall,'fP',0.1,'tauP',30,'netcon',eye(nCells,nCells)};

s.connections(end+1).direction='Off->SOnOff';
s.connections(end).mechanism_list={'PSC'};
s.connections(end).parameters={'gSYN',0.03,'tauR',EI_rise,'tauD',EI_fall,'fP',0.2,'tauP',80,'netcon',eye(nCells,nCells)};

% noise at relays
s.connections(end+1).direction='ROn->ROn';
s.connections(end).mechanism_list={'iNoise_V3'};
s.connections(end).parameters={'nSYN',0.015,'tauR_N',EE_rise,'tauD_N',EE_fall,'simlen',time_end/dt,'netcon',eye(nCells,nCells)}; 

% cross-channel inhibition
s.connections(end+1).direction='ROn->X';
s.connections(end).mechanism_list={'PSC'};
s.connections(end).parameters={'gSYN',0.012,'tauR',EE_rise,'tauD',EE_fall,'netcon',eye(nCells,nCells)};

s.connections(end+1).direction='X->ROn';
s.connections(end).mechanism_list={'PSC'};
s.connections(end).parameters={'gSYN',0.01,'tauR',XE_rise,'tauD',XE_fall,'ESYN',-80,'netcon',XRnetcon};

% apply TD->E and TD->S inhibition
s.connections(end+1).direction='TD->ROn';
s.connections(end).mechanism_list={'PSC'};
s.connections(end).parameters={'gSYN',0.015,'tauR',IE_rise,'tauD',IE_fall,'ESYN',-80,'netcon',eye(nCells,nCells)};

s.connections(end+1).direction='TD->ROff';
s.connections(end).mechanism_list={'PSC'};
s.connections(end).parameters={'gSYN',0.015,'tauR',IE_rise,'tauD',IE_fall,'ESYN',-80,'netcon',eye(nCells,nCells)};

s.connections(end+1).direction='TD->X';
s.connections(end).mechanism_list={'PSC'};
s.connections(end).parameters={'gSYN',0.015,'tauR',IE_rise,'tauD',IE_fall,'ESYN',-80,'netcon',eye(nCells,nCells)};

% convergence at output from onset-responding neurons
s.connections(end+1).direction='ROn->C';
s.connections(end).mechanism_list={'PSC'};
s.connections(end).parameters={'gSYN',0.02,'tauR',EE_rise,'tauD',EE_fall,'fP',0.1,'tauP',30,'netcon',RCnetcon};

%% vary params
vary = cell(length(varies),3);
for i = 1:length(varies)
    vary{i,1} = varies(i).conxn;
    vary{i,2} = varies(i).param;
    vary{i,3} = varies(i).range;
end

% do dsVary2Modifications here to save time and reduce # of opto
% simulations
I_ind = find(strcmp(vary(:,2),'Itonic')); FR_ind = find(matches(vary(:,2),'FR'));

if ~isempty(I_ind) && ~isempty(FR_ind)
if numel(vary{I_ind,3}) > 1 && size(vary{FR_ind,3},2) > 1
    numSets = numel(vary{I_ind,3})*size(vary{FR_ind,3},2);
    
    sqs = (1:numel(vary{I_ind,3})).^2;
    nonsqs = setdiff(1:numSets,sqs);
    
    vary = dsVary2Modifications(vary);
    
    temp = [];
    for n = nonsqs, temp = cat(2,temp,n:numSets:numel(vary)); end
    
    vary(temp) = [];
end
end

%% simulate
tic;

simdata = 0;
simdata = dsSimulate(s,'tspan',[dt time_end], 'solver',solverType, 'dt',dt,...
  'downsample_factor',1, 'save_data_flag',0, 'save_results_flag',1,...
  'study_dir',study_dir, 'vary',vary, 'debug_flag', 1, 'verbose_flag',0,...
  'mex_flag',options.mex_flag);

toc;

end