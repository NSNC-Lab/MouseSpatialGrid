function [options,varies] = params_4channel(gsynvec)

%% Options struct

options = struct;
options.nCells = 4;
options.opto = 0;

options.mex_flag = 0;
options.parfor_flag = 0;
options.plotRasters = 0;

% locations of speakers: [90 45 0 -90];

% locNum should be empty for full grids
options.locNum = [5 10 15 20];
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

% E->E connections
varies(end+1).conxn = '(On->ROn,Off->ROff)';
varies(end).param = 'gSYN';
varies(end).range = gsynvec(1);

% E->E connections
varies(end+1).conxn = 'ROn->C';
varies(end).param = 'gSYN';
varies(end).range = gsynvec(2);

% pv inputs
varies(end+1).conxn = '(SOn->ROn,SOn->ROff)';
varies(end).param = 'gSYN';
varies(end).range = gsynvec(3);

% offset pv inhibition
varies(end+1).conxn = '(SOff->ROn,SOff->ROff)';
varies(end).param = 'gSYN';
varies(end).range = gsynvec(4);

% inputs to SOM neurons
varies(end+1).conxn = '(ROn->X)';
varies(end).param = 'gSYN';
varies(end).range = gsynvec(5);

% cross-channel inhibition
varies(end+1).conxn = '(X->ROn)';
varies(end).param = '(gSYN,tauR,tauD,fF,tauF)';
varies(end).range = [gsynvec(6) ; 6; 12 ; 0.1 ; 120];

% noise at intermediate neurons
varies(end+1).conxn = 'ROn->ROn';
varies(end).param = 'FR';
varies(end).range = [8/options.nCells];

end
