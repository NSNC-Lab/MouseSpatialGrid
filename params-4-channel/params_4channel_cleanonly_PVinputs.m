%% Options struct

options = struct;
options.nCells = 4;
options.opto = 0;

options.mex_flag = 0;
options.parfor_flag = 0;
options.plotRasters = 0;

% locations of speakers: [90 45 0 -90];

% locNum should be empty for full grids
%options.locNum = [5 10 15 20];
options.locNum = [];
options.SpatialAttention = 0;

%% define network parameters
clear varies

if options.opto, nSims = 5; else, nSims = 1; end

trialInds = repmat(1:20,nSims,1);

% % % DO NOT CHANGE THIS % % %
varies(1).conxn = '(On->On,Off->Off)';
varies(1).param = 'trial';
varies(1).range =  trialInds(:)';
% % % DO NOT CHANGE THIS % % %

%Lets just start with changing the input gsyn E->E connections. 1 for each
%channel 4 total



% (On->ROn) connections
% varies(end+1).conxn = '(On->ROn)';
% varies(end).param = 'gSYN';
% varies(end).range = varied_params;



% noise at intermediate neurons
varies(end+1).conxn = 'ROn->ROn';
varies(end).param = 'FR';
varies(end).range = 8/options.nCells;
