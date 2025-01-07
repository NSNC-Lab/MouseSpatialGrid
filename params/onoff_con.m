%% Options struct

options = struct;
options.nCells = 1;
options.opto = 0;

if options.opto, nSims = 5; else, nSims = 1; end

options.mex_flag = 0;
options.parfor_flag = 0;
options.plotRasters = 0;

% locNum should be empty for full grids
options.locNum = 15;
options.SpatialAttention = 0;

%% define network parameters
clear varies

trialInds = repmat(1:20,nSims,1);

% % % DO NOT CHANGE THIS % % %
varies(1).conxn = '(On->On,Off->Off)';
varies(1).param = 'trial';
varies(1).range =  trialInds(:)';
% % % DO NOT CHANGE THIS % % %

% Input strength
varies(end+1).conxn = '(On->On,Off->Off)';
varies(end).param = 'g_postIC';
varies(end).range = 0.27;

%Switch back Off->R1On to On->R1On
    
% E->E connections
varies(end+1).conxn = '(On->RIOn)';
varies(end).param = 'gSYN';
varies(end).range = 0.02;

varies(end+1).conxn = '(Off->RIOn)';
varies(end).param = 'gSYN';
varies(end).range = 0.005;

varies(end+1).conxn = '(Off->RIOff)';
varies(end).param = 'gSYN';
varies(end).range = 0.02;

varies(end+1).conxn = '(RIOn->ROn)';
varies(end).param = 'gSYN';
varies(end).range = 0.02;

%PVs
varies(end+1).conxn = '(SIOnOff->RIOn)';
varies(end).param = 'gSYN';
varies(end).range = 0.025;

varies(end+1).conxn = '(SIOnOff->RIOff)';
varies(end).param = 'gSYN';
varies(end).range = 0.025;

varies(end+1).conxn = '(SOnOff->ROn)';
varies(end).param = 'gSYN';
varies(end).range = 0.025;

% On -> PV
varies(end+1).conxn = '(On->SIOnOff)';
varies(end).param = 'gSYN';
varies(end).range = 0.02;

varies(end+1).conxn = '(RIOn->SOnOff)';
varies(end).param = 'gSYN';
varies(end).range = 0.02;

% Off-> PV
varies(end+1).conxn = '(Off->SIOnOff)';
varies(end).param = 'gSYN';
varies(end).range = 0.005;

varies(end+1).conxn = '(RIOff->SOnOff)';
varies(end).param = 'gSYN';
varies(end).range = 0.005;



% control and opto conditions 
varies(end+1).conxn = '(SOnOff,SIOnOff)';
varies(end).param = 'Itonic';
varies(end).range = 0; 

varies(end+1).conxn = '(ROn->ROn)';
varies(end).param = 'FR';
varies(end).range = 8;



%Figure 4 (Sweep over fp values for E->E)
% varies(end+1).conxn = '(On->ROn,Off->ROn)';
% varies(end).param = 'fP';
% varies(end).range = [0:0.1:1];

% %Sweep over fp values for E->PV)
% varies(end+1).conxn = '(On->SOnOff,Off->SOnOff)';
% varies(end).param = 'fP';
% varies(end).range = [0:0.1:1];

%Sweep over fp values for PV->E)
varies(end+1).conxn = '(SOnOff->ROn)';
varies(end).param = 'fP';
varies(end).range = [0:0.1:1];

% find varied parameter, excluding trials
varied_param = find( (cellfun(@length,{varies.range}) > 1 & ~cellfun(@iscolumn,{varies.range})));
if numel(varies(1).range) == 20, varied_param(1) = []; % delete trial varies
    end

if isempty(varied_param) % if no varied params, settle on 2nd entry in varies
    varied_param = 2;
end

% netcons can't be put in the varies struct (Dynasim doesn't recognize it?);
% it needs to be put in a different struct

netcons = struct; % row = source, column = target
netcons.XRnetcon = zeros(options.nCells,options.nCells);

% for simplification: use 1st varied param for 2d searches
expVar = [varies(varied_param(1)).conxn '-' varies(varied_param(1)).param];
expVar = strrep(expVar,'->','_');