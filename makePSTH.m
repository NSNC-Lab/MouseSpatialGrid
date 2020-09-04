function [PSTH,t_vec] = makePSTH(raster)

% calculate PSTH of model results

if size(raster,2) < 3000    % w/o spontaneous activity
    t_vec = 1:20:size(raster,2);
else % w/ spontaneous activity
    t_vec = 1:20:5000;
end

% add trials together
temp = sum(raster);

% sum up spikes within each time bin
for t = 1:length(t_vec)-1
    PSTH(t) = sum(temp(t_vec(t):t_vec(t+1)));   
end
PSTH(end+1) = sum(temp(t_vec(end):end));

t_vec = (t_vec-1)/1000;

end