function [PSTH,t_vec] = makeplotPSTHs(varargin)

% outputs:
% PSTH = [numTargets x time] matrix
% t_vec = [1 x time] vector

raster = varargin{1};
if length(varargin) > 1
    numTargets = varargin{2};
else
    numTargets = 1;
end

% calculate PSTH of model results

% bin by 20 ms, assuming dt = 1
t_vec = 1:20:size(raster,2);

numTrials = size(raster,1)/numTargets; % number of trials to add for PSTH

for n = 1:numTargets
    temp = sum(raster(([1:numTrials]+numTrials*(n-1)),:));
    
    % sum up spikes within each time bin
    
    for t = 1:length(t_vec)-1
        PSTH(n,t) = sum(temp(t_vec(t):t_vec(t+1)));
    end
    PSTH(n,length(t_vec)) = sum(temp(t_vec(end):end));
        
    plot((t_vec-1)/1000,50*PSTH(n,:)/numTrials + (150*(n-1)));
    hold on;
end

t_vec = (t_vec-1)/1000;

if numTargets == 2
    plot([t_vec(1) t_vec(end)],[150 150],'k');
end

set(gca,'ytick',[0:50:300]);
set(gca,'yticklabels',{'0','50','100','0','50','100','150'});
ylim([0 300]);
ylabel('FR (Hz)');
xlabel('Time (s)');

% output t_vec in terms of s

end