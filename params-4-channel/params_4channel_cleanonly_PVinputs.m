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
%options.locNum = [];
options.locNum = 1:24;
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

NetconHandler;
GsynHandler;
% XRnetcon = netcons.XRnetcon;
% RCnetcon = netcons.RCnetcon;
% PEnetcon = netcons.PEnetcon;
% 
% XRgsyncon = gsyncons.XRgsyncon;
% RCgsyncon = gsyncons.RCgsyncon;
% OnRgsyncon = gsyncons.OnRgsyncon;
% PEgsyncon = gsyncons.PEgsyncon;

%Just going to use this to do our average of 3.
% Input strength
% varies(end+1).conxn = '(On->On,Off->Off)';
% varies(end).param = 'g_postIC';
% varies(end).range = 0.165*ones(1,3);



% (On->ROn) connections
% varies(end+1).conxn = '(On->ROn)';
% varies(end).param = 'gSYN';
% varies(end).range = varied_params;

% % PV->E strength
% varies(end+1).conxn = 'SOnOff->ROn';
% varies(end).param = 'gSYN';
% varies(end).range = {PEgsyncon,PEgsyncon};
% 
% 
% % ROn -> C strength
% varies(end+1).conxn = 'ROn->C';
% varies(end).param = 'gSYN';
% varies(end).range = {RCgsyncon,RCgsyncon};


% noise at intermediate neurons
varies(end+1).conxn = 'ROn->ROn';
varies(end).param = 'FR';
varies(end).range = 8/options.nCells;
