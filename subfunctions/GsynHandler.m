% gsyn networks enable pathway specific synaptic strength modification
% row = source, column = target 
gsyncons = struct; 

% all other gsyn networks defined as identity in columnNetwork_simpler_onoff

% OnRgsyncon: Onset (On) -> ROn
% default was 0.02
gsyncons.OnRgsyncon = [0.015,0.015,0.015,0.015];

% XRgsyncon: SOM (X) -> E (Ron)
gsyncons.XRgsyncon =   [0.000,0.002,0.002,0.002;
                        0.002,0.000,0.000,0.000;
                        0.002,0.000,0.000,0.000;
                        0.002,0.000,0.000,0.000];

% RCgsyncon: E (Ron) -> C
% needs to be transposed for PSC3 mech
gsyncons.RCgsyncon = [0.01, 0.000, 0.000, 0.000].';


gsyncons.PEgsyncon = [0.025,0.0,0.00,0.0;
                      0.0,0.025,0.000,0.000;
                      0.0,0.000,0.025,0.000;
                      0.0,0.000,0.000,0.025];
