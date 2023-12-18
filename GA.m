%Load in target data
load('90_deg_target.mat')

%Initilize the GA
%Run SpikingNetowkr_withOffset with random inputs.
gsynvec = rand(8,1);
SpikingNetwork_withOffset%(gsynvec)