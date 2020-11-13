function sigIn = genICinput_V3(trial,dt,label)
% input:
% output:
%   sigIn = a vector with 1s at spike times waveforms, with dimensions time x (nfreqs x nlocs)
%
% @ Erik Roberts, Kenny Chou
% Boston Univeristy, 2019
%
% TODO: allow resampling of spk_IC to different dt

if exist(['IC_spks_' label(2:end-1) '.mat'],'file')
    fileData = load(['IC_spks_' label(2:end-1) '.mat'],'spks');
else
    fileData = load(['..\IC_spks_' label(2:end-1) '.mat'],'spks');
end

% IC data: trial x location x time
% Desired: time x location x trial
% Permute data 
spk_IC = fileData.spks;
spk_IC = permute(spk_IC,[3,2,1]);
sigIn = squeeze(spk_IC(:,:,trial)); % time x location x cells

end