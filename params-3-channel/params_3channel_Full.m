%% Options struct

options = struct;
options.nCells = 3;
options.opto = 0;

options.mex_flag = 0;
options.parfor_flag = 0;
options.plotRasters = 0;

% locNum should be empty for full grids
options.locNum = 5:24;
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
varies(1).conxn = '(On->On,Off->Off)';
varies(1).param = 'trial';
varies(1).range =  trialInds(:)';
% % % DO NOT CHANGE THIS % % %

% E->E connections
varies(end+1).conxn = '(On->R1On,R1On->R2On,Off->R1Off,R1Off->R2Off)';
varies(end).param = 'gSYN';
varies(end).range = 0.02;

% R2On->C connections
varies(end+1).conxn = 'R2On->C';
varies(end).param = 'gSYN';
varies(end).range = 0.015;

% onset pv inhibition
varies(end+1).conxn = '(S1On->R1On,S1On->R1Off,S2On->R2On,S2On->R2Off)';
varies(end).param = 'gSYN';
varies(end).range = 0.025;

% offset pvs inhibition
varies(end+1).conxn = '(S1Off->R1On,S1Off->R1Off,S2Off->R2On,S2Off->R2Off)';
varies(end).param = 'gSYN';
varies(end).range = 0.01;

% cross-channel inhibition
varies(end+1).conxn = '(X1->R1On,X2->R2On)';
varies(end).param = 'gSYN';
varies(end).range = 0.013;

% noise at intermediate neurons
varies(end+1).conxn = 'R2On->R2On';
varies(end).param = 'FR';
varies(end).range = 3;

