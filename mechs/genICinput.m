function sigIn = genICinput(trial,tauR,tauD,dt)
% input:
% output:
%   sigIn = a vector with 1s at spike times waveforms, with dimensions time x (nfreqs x nlocs)
%
% @ Erik Roberts, Kenny Chou
% Boston Univeristy, 2019
%
% TODO: allow resampling of spk_IC to different dt

if exist('IC_spks.mat','file')
    fileData = load('IC_spks.mat','spks');
else
    fileData = load('..\IC_spks.mat','spks');
end

% IC data: trial x location x time
% Desired: time x location x trial
% Permute data 
spk_IC = fileData.spks;
spk_IC = permute(spk_IC,[3,2,1]);
sigIn = squeeze(spk_IC(:,:,trial)); % time x location x cells

% ========================= create EPSP waveform =========================
% ============= convolve each time series with epsc waveform ==============
t_a = max(tauR,tauD)*7; % Max duration of syn conductance
t_vec = 0:dt:t_a;
tau2 = tauR;
tau1 = tauD;
tau_rise = tau1*tau2/(tau1-tau2);
b = ((tau2/tau1)^(tau_rise/tau1) - (tau2/tau1)^(tau_rise/tau2)^-1); % /tau2?
epsc =  - b * ( exp(-t_vec/tau1) - exp(-t_vec/tau2) ); % - to make positive
sigIn = conv2(sigIn,epsc','same');
% 
% t_a = dur; % Max duration of syn conductance
% t_vec = 0:dt:t_a;
% current =  ones(1,length(t_vec));
% sigIn = conv2(sigIn,current','same');

end