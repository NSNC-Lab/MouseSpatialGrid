function plotPerformanceGrids_v3(data,s,annotTable,subPops,targetIdx,mixedIdx,simOptions,expName)
% plot performance grids for specified subpopulation of neurons
%
% v3 changes the script into a function

subz = simOptions.subz;
locationLabels = simOptions.locationLabels;

% setup populations
popNames = {s.populations.name};
% subPops = popNames(~contains(popNames,'TD') & ~contains(popNames,'X')); %remove TD neuron
popNamesT = strcat({s.populations.name},'T');
popNamesM = strcat({s.populations.name},'M');
numPops = numel(subPops);
popSizes = [s.populations.size];
popSizes = popSizes(matches(popNames,subPops));
onlyC = any(contains(subPops,'C')) & length(subPops) == 1;

% rearrange subPops if C is not the only thing in subPops
if ~onlyC
    C_ind = find(contains(subPops,'C'));
    subPops = subPops([setdiff(1:length(subPops),C_ind) , C_ind]);
end

chanLabels = simOptions.chanLabels;

% check if varied parameter is vector or matrix
numVars = length(annotTable);

for vv = 1:numVars
    annotStr = annotTable{vv};
    
    % figure setup
    figwidth = max(length(popSizes)*300,900); % 300 * number of neurons?
    figheight = min(length(popSizes)*300,800);
    h = figure('position',[200 50 figwidth figheight]);
    
    plotheight = 1/numPops*0.55; % percentage
    
    xstart = 0.1;
    ystart = 0.1;
    for pop = 1:numPops
        
        % vary grid width depending on the number of neurons in subpopulation
        if popSizes(pop) == 1
            plotwidth = 0.233; % percentage, 0.7*num neurons
        else
            plotwidth = 0.14;
        end
        
        for chan = 1:popSizes(pop)

            % subfigure positions
            col = chan-1;
            row = pop-1;
            xoffset = xstart+plotwidth*(col)*1.2;
            yoffset = ystart+plotheight*1.6*row;

            % mixed cases
            if sum(ismember(subz,mixedIdx)) > 0
                for i = 1:length(mixedIdx)
                    idx = subz(mixedIdx(i) == subz);
                    perf.(subPops{pop})(i) = data(idx).perf.(subPops{pop}).(['channel' num2str(chan)])(vv);
                    fr.(subPops{pop})(i) = data(idx).fr.(subPops{pop}).(['channel' num2str(chan)])(vv);
                end
                subplot('Position',[xoffset yoffset plotwidth plotheight])
                plotPerfGrid(perf.(subPops{pop})',fr.(subPops{pop})',[]);
                
                % add axes
                if onlyC || (chan==1 && pop==1)
                    xticks(1:4); xticklabels(locationLabels)
                    yticks(1:4); yticklabels(fliplr(locationLabels))
                    xlabel('target location')
                    ylabel('masker location')
                end
            end

            % target-only configs
            if sum(ismember(subz,targetIdx)) > 0
                perf.(popNamesT{pop}) = zeros(1,4);
                fr.(popNamesT{pop}) = zeros(1,4);
                
                if ~isempty(targetIdx)
                    for i = 1:length(targetIdx)
                        idx = subz(targetIdx(i) == subz);
                        perf.(popNamesT{pop})(i) = data(idx).perf.(subPops{pop}).(['channel' num2str(chan)])(vv);
                        fr.(popNamesT{pop})(i) = data(idx).fr.(subPops{pop}).(['channel' num2str(chan)])(vv);
                    end
                end

                subplot('Position',[xoffset yoffset+plotheight+0.01 plotwidth plotheight*0.2])
                plotPerfGrid([perf.(popNamesT{pop})],[fr.(popNamesT{pop})],'');

                if chan == 1, ylabel(subPops{pop}); end
            end

            % if we're plotting the population before the C cell (R2On),
            % add channel label titles
            if popSizes(pop) > 1, title(chanLabels{chan},'fontweight','normal'); end
        end
    end
    
%     % simulation info
%     annotation('textbox',[xoffset+plotwidth*1.2 yoffset+plotheight*0.5 plotwidth plotheight],...
%            'string',annotStr,...
%            'FitBoxToText','on',...
%            'LineStyle','none')
       
    saveas(gcf,fullfile('simData',expName,['C_grid_vary' num2str(vv) '.png']));
end