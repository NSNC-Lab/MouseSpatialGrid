% NETCONS enable connectivity between different channels of the network
% row = source, column = target 
netcons = struct; 

% XRnetcon: SOM (X) -> E (Ron)
netcons.XRnetcon = [1,1,0,0;
                    0,1,0,0;
                    0,0,1,0;
                    0,0,0,1];

% PEnetcon: PV (SOnOff) -> E (Ron)
netcons.PEnetcon = eye(options.nCells);

% RCnetcon: E (Ron) -> C
netcons.RCnetcon = [1;1;1;1];

