function [simdata,s] = columnNetwork_paper(study_dir,varies,options,netcons)

% Generates and simulates a network featuring columns of excitatory cells 
% that respond to onsets and offsets in auditory stimuli

% Features of circuit:

% Separate PV cells for onsets and offsets in stimuli; the PV cells will
% inhibit both onset and offset relay units in the same layer

% This is the original version of the single-channel model published in
% bioRxiv (Nocon et al 2022)

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

% onset column

s.populations(1).name = 'On';
s.populations(1).equations = 'noconLIF';
s.populations(1).size = nCells;
s.populations(1).parameters = {'t_ref',1};

s.populations(end+1).name = 'Off';
s.populations(end).equations = 'noconLIF';
s.populations(end).size = nCells;
s.populations(end).parameters = {'t_ref',1};

s.populations(end+1).name='R1On';
s.populations(end).equations = 'noconLIF';
s.populations(end).size = nCells;
s.populations(end).parameters = {'g_inc',0.0003,'tau_ad',100,'t_ref',1};

s.populations(end+1).name='R1Off';
s.populations(end).equations = 'noconLIF';
s.populations(end).size = nCells;
s.populations(end).parameters = {'g_inc',0.0003,'tau_ad',100,'t_ref',1};

s.populations(end+1).name='S1On';
s.populations(end).equations = 'noconLIF';
s.populations(end).size = nCells;
s.populations(end).parameters = {'g_L',1/100,'E_L',-57,'V_reset',-52,'t_ref',0.5};

s.populations(end+1).name='S1Off';
s.populations(end).equations = 'noconLIF';
s.populations(end).size = nCells;
s.populations(end).parameters = {'g_L',1/100,'E_L',-57,'V_reset',-52,'t_ref',0.5};

s.populations(end+1).name='R2On';
s.populations(end).equations = 'noconLIF';
s.populations(end).size = nCells;
s.populations(end).parameters = {'g_inc',0.0003,'tau_ad',100,'t_ref',1};

s.populations(end+1).name='R2Off';
s.populations(end).equations = 'noconLIF';
s.populations(end).size = nCells;
s.populations(end).parameters = {'g_inc',0.0003,'tau_ad',100,'t_ref',1};

s.populations(end+1).name='S2On';
s.populations(end).equations = 'noconLIF';
s.populations(end).size = nCells;
s.populations(end).parameters = {'g_L',1/100,'E_L',-57,'V_reset',-52,'t_ref',0.5};

s.populations(end+1).name='S2Off';
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

s.connections(end+1).direction='Off->Off';
s.connections(end).mechanism_list={'IC'};
s.connections(end).parameters={'g_postIC',0.265,'label','off','trial',1,'locNum',options.locNum,'netcon',eye(nCells,nCells),'t_ref',1,'t_ref_rel',1,'rec',2};

% % % L4 % % %

% excitatory inputs
s.connections(end+1).direction='On->R1On';
s.connections(end).mechanism_list={'PSC'};
s.connections(end).parameters={'gSYN',0.02,'tauR',EE_rise,'tauD',EE_fall,'fP',0.1,'tauP',30};

s.connections(end+1).direction='On->S1On';
s.connections(end).mechanism_list={'PSC'};
s.connections(end).parameters={'gSYN',0.03,'tauR',EI_rise,'tauD',EI_fall,'fP',0.2,'tauP',80};

s.connections(end+1).direction='S1On->R1On';
s.connections(end).mechanism_list={'PSC'};
s.connections(end).parameters={'gSYN',0.03,'tauR',IE_rise,'tauD',IE_fall,'ESYN',-80,'fP',0.5,'tauP',120}; 

s.connections(end+1).direction='S1On->R1Off';
s.connections(end).mechanism_list={'PSC'};
s.connections(end).parameters={'gSYN',0.03,'tauR',IE_rise,'tauD',IE_fall,'ESYN',-80,'fP',0.5,'tauP',120}; 

% offset channels
s.connections(end+1).direction='Off->R1Off';
s.connections(end).mechanism_list={'PSC'};
s.connections(end).parameters={'gSYN',0.02,'tauR',EE_rise,'tauD',EE_fall,'fP',0.1,'tauP',30};

s.connections(end+1).direction='Off->S1Off';
s.connections(end).mechanism_list={'PSC'};
s.connections(end).parameters={'gSYN',0.03,'tauR',EI_rise,'tauD',EI_fall,'fP',0.2,'tauP',80};

s.connections(end+1).direction='S1Off->R1On';
s.connections(end).mechanism_list={'PSC'};
s.connections(end).parameters={'gSYN',0.015,'tauR',IE_rise,'tauD',IE_fall,'ESYN',-80,'fP',0.5,'tauP',120}; 

s.connections(end+1).direction='S1Off->R1Off';
s.connections(end).mechanism_list={'PSC'};
s.connections(end).parameters={'gSYN',0.015,'tauR',IE_rise,'tauD',IE_fall,'ESYN',-80,'fP',0.5,'tauP',120}; 

% % %  L2/3  % % %

% excitatory inputs
s.connections(end+1).direction='R1On->R2On';
s.connections(end).mechanism_list={'PSC'};
s.connections(end).parameters={'gSYN',0.02,'tauR',EE_rise,'tauD',EE_fall,'fP',0.1,'tauP',30};

s.connections(end+1).direction='R1On->S2On';
s.connections(end).mechanism_list={'PSC'};
s.connections(end).parameters={'gSYN',0.03,'tauR',EI_rise,'tauD',EI_fall,'fP',0.2,'tauP',80};

s.connections(end+1).direction='S2On->R2On';
s.connections(end).mechanism_list={'PSC'};
s.connections(end).parameters={'gSYN',0.025,'tauR',IE_rise,'tauD',IE_fall,'ESYN',-80,'fP',0.5,'tauP',120}; 

s.connections(end+1).direction='S2On->R2Off';
s.connections(end).mechanism_list={'PSC'};
s.connections(end).parameters={'gSYN',0.01,'tauR',IE_rise,'tauD',IE_fall,'ESYN',-80,'fP',0.5,'tauP',120}; 

% offset channels
s.connections(end+1).direction='R1Off->R2Off';
s.connections(end).mechanism_list={'PSC'};
s.connections(end).parameters={'gSYN',0.02,'tauR',EE_rise,'tauD',EE_fall,'fP',0.1,'tauP',30};

s.connections(end+1).direction='R1Off->S2Off';
s.connections(end).mechanism_list={'PSC'};
s.connections(end).parameters={'gSYN',0.03,'tauR',EI_rise,'tauD',EI_fall,'fP',0.2,'tauP',80};

s.connections(end+1).direction='S2Off->R2On';
s.connections(end).mechanism_list={'PSC'};
s.connections(end).parameters={'gSYN',0.025,'tauR',IE_rise,'tauD',IE_fall,'ESYN',-80,'fP',0.5,'tauP',120}; 

s.connections(end+1).direction='S2Off->R2Off';
s.connections(end).mechanism_list={'PSC'};
s.connections(end).parameters={'gSYN',0.01,'tauR',IE_rise,'tauD',IE_fall,'ESYN',-80,'fP',0.5,'tauP',120}; 

s.connections(end+1).direction='R2On->R2On';
s.connections(end).mechanism_list={'iNoise_V3'};
s.connections(end).parameters={'nSYN',0.015,'tauR_N',EE_rise,'tauD_N',EE_fall,'simlen',options.time_end / dt}; 

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

simdata = dsSimulate(s,'tspan',[dt time_end], 'solver',solverType, 'dt',dt,...
  'downsample_factor',1, 'save_data_flag',0, 'save_results_flag',1,...
  'study_dir',study_dir, 'vary',vary, 'debug_flag', 1, 'verbose_flag',0,...
  'mex_flag',options.mex_flag);

toc;

end