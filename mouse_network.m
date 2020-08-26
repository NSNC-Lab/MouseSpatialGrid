function simdata = mouse_network(study_dir,time_end,varies,restrict_vary_flag)
% [perf, fr, annotstr, distMat, VR] = mouse_network(study_dir,time_end,...
%    varies,plot_rasters,plot_distances,data_spks,data_tau,restrict_vary_flag)

% [performance, tauMax] = mouse_network(study_dir,time_end,varies,plot_rasters,data_spks)
% study_dir: location of IC spike files + directory for log and data files
% time_end: length of simulation in ms
% varies: vary parameter, a structure. e.g.
%   varies(1).conxn = '(IC->IC)';
%   varies(1).param = 'trial';
%   varies(1).range = 1:20;
%   *first set of parameters should always be "trial" for IC->IC cxn
% plot_rasters: 1 or 0
% data_spks: spikes from data we want to match 
%
% @Kenny F Chou, Boston Univ. 2019-06-20
% 2019-08-04 - added sharpening neurons
% 2019-08-14 - plotting now handles multiple varied parameters
% 2019-09-11 - removed redundant code. Return R performance in addition to
%              C performance
% 2020-02-05 - added netCons as parameter
% to do - fix plot_rasters option

% @Jio Nocon, Boston Univ., 2020-6-18
% 2020-6-18 - added subfunction to plot model rasters vs. data rasters and
%             calculate VR distance between the two

%% Input check
if ~strcmp(varies(1).param,'trial')
    error('first set of varied params should be ''trial''')
end

%% solver params
solverType = 'euler';
dt = 1; %ms % the IC input is currently dt=1

%% neuron populations
% tonic = bias = cells spontaneous firing

nCells = 4;
noise = 0.01; % low noise
s = struct();

s.populations(1).name = 'IC';
s.populations(1).equations = 'chouLIF';
s.populations(1).size = nCells;
s.populations(1).parameters = {'Itonic',0,'noise',0}; % 10-20 Hz spiking at rest

% s.populations(end+1).name = 'X';
% s.populations(end).equations = 'chouLIF';
% s.populations(end).size = nCells;
% s.populations(end).parameters = {'Itonic',0, 'noise',0}; % 10-20 Hz spiking at rest

% s.populations(end+1).name = 'S'; % S for sharpening
% s.populations(end).equations = 'chouLIF';
% s.populations(end).size = nCells;
% s.populations(end).parameters = {'Itonic',0, 'noise',0}; % 10-20 Hz spiking at rest

s.populations(end+1).name='R';
s.populations(end).equations = 'chouLIF';
s.populations(end).size = nCells;
s.populations(end).parameters = {'Itonic',0,'noise',0};

s.populations(end+1).name='C';
s.populations(end).equations = 'chouLIF';
s.populations(end).size = 1;
s.populations(end).parameters = {'noise',noise};

%% connections

s.connections(1).direction='IC->IC';
s.connections(1).mechanism_list={'IC_V2'};
s.connections(1).parameters={'g_postIC',0.25,'trial',5}; % 100 hz spiking

% s.connections(end+1).direction='IC->X';
% s.connections(end).mechanism_list={'synDoubleExp'};
% s.connections(end).parameters={'gSYN',0.21, 'tauR',0.4, 'tauD',2, 'netcon',diag(ones(1,nCells))}; 

% s.connections(end+1).direction='IC->S';
% s.connections(end).mechanism_list={'synDoubleExp'};
% s.connections(end).parameters={'gSYN',0.21, 'tauR',0.4, 'tauD',2, 'netcon',diag(ones(1,nCells))}; 

s.connections(end+1).direction='IC->R';
s.connections(end).mechanism_list={'synDoubleExp'};
s.connections(end).parameters={'gSYN',0.21, 'tauR',0.3, 'tauD',1.5, 'netcon', diag(ones(1,nCells)),'delay',0}; 

% s.connections(end+1).direction='X->R';
% s.connections(end).mechanism_list={'synDoubleExp'};
% s.connections(end).parameters={'gSYN',0.21, 'tauR',0.4, 'tauD',2, 'netcon',xrNetcon, 'ESYN',-80}; 
% 
% s.connections(end+1).direction='S->R';
% s.connections(end).mechanism_list={'synDoubleExp'};
% s.connections(end).parameters={'gSYN',0.21, 'tauR',0.4, 'tauD',5, 'netcon',srNetcon, 'ESYN',-80, 'delay',3}; 

s.connections(end+1).direction='R->C';
s.connections(end).mechanism_list={'synDoubleExp_V2'};
s.connections(end).parameters={'inputChan1',1,'inputChan2',1,'inputChan3',1,'inputChan4',1}; 

% if viz_network, vizNetwork; end

%% vary params
vary = cell(length(varies),3);
for i = 1:length(varies)
    vary{i,1} = varies(i).conxn;
    vary{i,2} = varies(i).param;
    vary{i,3} = varies(i).range;
end

%% simulate
tic;

simdata = dsSimulate(s,'tspan',[dt time_end], 'solver',solverType, 'dt',dt,...
  'downsample_factor',1, 'save_data_flag',0, 'save_results_flag',1,...
  'study_dir',study_dir, 'vary',vary, 'debug_flag', 0, 'verbose_flag',0,...
  'restrict_vary_flag',restrict_vary_flag);

simdata = rmfield(simdata,{'IC_V','IC_g_ad','R_V','R_g_ad','labels','simulator_options','model'});

% save(fullfile(study_dir,'simulation_results.mat'),'simdata','-v7.3');

toc;

end
