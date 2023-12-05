function [pc,fr] = plotParamvsPerf_1D(varargin)

results = varargin{1};
nVaries = varargin{2}; % how much parameter sets are there, excluding the trials and number of repeat trials for laser
dt = varargin{3}; 

if nargin == 4
    pop = varargin{4};
else
    fnames = fieldnames(results);
    pop = fnames(contains(fnames,'_V_spikes'));
    pop = erase(pop{1},'_V_spikes');
end

nSims = length(results)/nVaries/20;

for ns = 1:nSims
    for n = 1:nVaries

        % fetch 20 trials to calculate performance at specific spot
        subData = results( ((ns-1)*nVaries + n) : nVaries*nSims : end);

        nCh = size(subData(1).([pop '_V_spikes']),2);

        for ch = 1:nCh

            for t = 1:20
                raster(t,:) = subData(t).([pop '_V_spikes'])(:,ch);
            end

            [pc.SPIKE(ns,n,ch),pc.ISI(ns,n,ch),pc.RISPIKE(ns,n,ch),pc.spkct(ns,n),fr(ns,n,ch)] = calcPCandFR(raster,20,dt);

        end

    end
end

end


function [pc_SPIKE,pc_ISI,pc_RISPIKE,pc_spkct,fr] = calcPCandFR(raster,numTrials,dt)

start_time = 300; % [ms]
end_time = start_time + 3000; % [ms];

% spks to spiketimes in a cell array of 20x2
spkTime = cell(numTrials,1);
for ii = 1:numTrials
    spkTime{ii} = find(raster(ii,:))*dt;
end
spkTime = reshape(spkTime,numTrials/2,2);
input = reshape(spkTime,1,numTrials);
fr = round(mean(cellfun(@(x) sum(x >= start_time & x < end_time) / 3,input)));

STS = SpikeTrainSet(input,start_time,end_time);

distMat = STS.SPIKEdistanceMatrix(start_time,end_time);
pc_SPIKE = calcpcStatic(distMat, numTrials/2, 2, 0);

distMat = STS.ISIdistanceMatrix(start_time,end_time);
pc_ISI = calcpcStatic(distMat, numTrials/2, 2, 0);

distMat = STS.RateIndependentSPIKEdistanceMatrix(start_time,end_time);
pc_RISPIKE = calcpcStatic(distMat, numTrials/2, 2, 0);

distMat = calcSpkCtDist(spkTime,start_time,end_time);
pc_spkct = calcpcStatic(distMat, numTrials/2, 2, 0);

end
