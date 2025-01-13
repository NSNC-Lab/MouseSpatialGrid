function [s, to, fo, pg] = GaussianSpectrum(input, increment, winLength, samprate)
%
% Gaussian spectrum
% 	s = GuassianSpectrum(input, increment, winLength, samprate)
% 	Compute the guassian spectrogram of an input signal with a guassian
% 	window that is given by winLength. The standard deviation of that
% 	guassian is 1/6th of winLength.
%	Each time frame is [winLength]-long and
%	starts [increment] samples after previous frame's start.
%	Only zero and the positive frequencies are returned.
%   to and fo are the time and frequency for each bin in s and Hz
%   pg is a rumming rms.

%%%%%%%%%%%%%%%%%%%%%%%
% Massage the input
%%%%%%%%%%%%%%%%%%%%%%%

% Enforce even winLength to have a symmetric window
if rem(winLength, 2) == 1
    winLength = winLength +1;
end

% Make input it into a row vector if it isn't
if size(input, 1) > 1,
	input = input';
end;

% Padd the input with zeros

%IB 1/7/2025
%Here it seems that you would padd with the winLength instead of
%winlength/2. If you only pad with half of the window length then you risk
%alterin the signal's latency. For instance the window will have some
%non-padded values for the first value calculated by the stft. This means
%that the first value that you output in your stft is your n = winlength/2
%value at the output.

%Using matlabs fft function puts your lower freqncies first. I think that
%you might be able to use fftshift in order to correct for this, but by
%padding the front you can fix the un-causal relationships seen in the
%model.

%You musta also account for the increment which will change the window in
%steps of size 19. Therefore in order to stay causal you must also push the
%window by 19.

%pinput = zeros(1,length(input)+winLength);
%pinput(winLength/2+1:winLength/2+length(input)) = input;
%inputLength = length(pinput);

%<><><><><><><><>%
pinput = zeros(1,length(input)+increment+2*winLength);
pinput(winLength+1+increment:winLength+increment+length(input)) = input;
inputLength = length(pinput);
%<><><><><><><><>%

% The number of time points in the spectrogram
frameCount = floor((inputLength-winLength)/increment)+1;

% The window of the fft
fftLen = winLength;


%%%%%%%%%%%%%%%%%%%%%%%%
% Guassian window 
%%%%%%%%%%%%%%%%%%%%%%%%
nstd = 6;                   % Number of standard deviations in one window.
wx2 = ((1:winLength)-((winLength+1)/2)).^2;
wvar = (winLength/nstd)^2;
ws = exp(-0.5*(wx2./wvar));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize output "s" 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if rem(fftLen, 2)   
    % winLength is odd
    s = zeros((fftLen+1)/2+1, frameCount);
else
    % winLength is even 
    s = zeros(fftLen/2+1, frameCount);
end

pg = zeros(1, frameCount);
for i=1:frameCount
    start = (i-1)*increment + 1;
    last = start + winLength - 1;
    f = zeros(fftLen, 1);
    f(1:winLength) = ws.*pinput(start:last);
    pg(i) = std(f(1:winLength));
    %\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\\/\/\/\/\/\/\/\/\/
    specslice = fft(f);
    if rem(fftLen, 2)   % winLength is odd
        s(:,i) = specslice(1:((fftLen+1)/2+1));
    else
        s(:,i) = specslice(1:(fftLen/2+1));
    end
    %/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\\/\/\/\/\/\/\/\/\/

    %s(:,i) = specslice(1:(fftLen/2+1));
end

% Assign frequency_label
if rem(fftLen, 2)   % winLength is odd
    select = 1:(fftLen+1)/2;
else
    select = 1:fftLen/2+1;
end
fo = (select-1)'*samprate/fftLen;

% assign time_label
to = ((1:size(s,2))-1)'.*(increment/samprate);
return

