% gsyn networks enable pathway specific synaptic strength modification
% row = source, column = target 
gsyncons = struct; 

% all other gsyn networks defined as identity in columnNetwork_simpler_onoff

% OnRgsyncon: Onset (On) -> ROn
% default was 0.02
gsyncons.OnRgsyncon = [0.020,0.020,0.020,0.020];

% XRgsyncon: SOM (X) -> E (Ron)
gsyncons.XRgsyncon =   [0.000,0.000,0.000,0.000;
                        0.005,0.000,0.000,0.000;
                        0.005,0.000,0.000,0.000;
                        0.005,0.000,0.000,0.000];

% RCgsyncon: E (Ron) -> C
% needs to be transposed for PSC3 mech
% gsyncons.RCgsyncon = [varied_struct.RtoC, 0.003, 0.003, 0.003].';

gsyncons.RCgsyncon = [0.006, 0.006, 0.006, 0.006].';


% gsyncons.PEgsyncon = [varied_struct.IntraPV,0.000,0.000,0.000;
                                    % 0.000,0.020,0.000,0.000;
                                    % 0.000,0.000,0.020,0.000;
                                    % 0.000,0.000,0.000,0.020];

gsyncons.PEgsyncon = [0.025,0.000,0.000,0.000;
                      0.000,0.025,0.000,0.000;
                      0.000,0.000,0.025,0.000;
                      0.000,0.000,0.000,0.025];

