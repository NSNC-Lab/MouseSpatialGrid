function makeGrids_bestIteration(data,varies,DirPart,data_perf,data_FR,best_iterations,loss)

set(0,'defaultfigurevisible','on');
    
% performance vector has dimensions [numSpatialChan,nvaried]
neurons = {'Ipsi. sigmoid','Gaussian','U','Cont. sigmoid'};

temp = {data.name};
temp(cellfun('isempty',temp)) = {'empty'}; %label empty content

targetIdx = find(contains(temp,'m0') & ~strcmp(temp,'s0m0.mat'));
%maskerIdx = find(contains(temp,'s0') & ~strcmp(temp,'s0m0.mat'));
%mixedIdx = find(~contains(temp,'m0') & ~contains(temp,'s0') & ~contains(temp,'empty'));
textColorThresh = 70;
numSpatialChan = 4;

figstr = 'CleanGrid vary ';

% if more than one parameter is varied
if sum(cellfun(@length,{varies(2:end).range}) > 1) > 1
    
    multiFlag = 1;
    
    temp = find((cellfun(@length,{varies.range}) > 1) == 1);
    temp(1) = [];
    variedConxns = {varies(temp).conxn};
    variedParams = {varies(temp).param};
    variedRanges = {varies(temp).range};
    
    for i = 1:length(temp)
        figstr = cat(2,figstr,[variedConxns{i},variedParams{i},'%0.3f, ']);
    end
    figstr(end-1:end) = []; % delete extra ', '
    
    [w1,w2] = meshgrid(variedRanges{1},variedRanges{2});
    c = cat(2,w1,w2);
    paramPairs = reshape(c,[],2);
else
    multiFlag = 0;
    ind = find(cellfun(@length,{varies(2:end).range}) ~= 1);
    variedParam = [varies(ind+1).conxn,varies(ind+1).param];
    figstr = cat(2,figstr,variedParam,'%0.3f');
end

width=8; hwratio=1;
x0=.08; y0=.08;
dx=.04; dy=.04;
lx=.125; ly=.125/hwratio;

x=-108:108;
tuningcurve = zeros(4,length(x));
ono = load('ono_curves_V2.mat','sigmoid','gauss','ushaped');

gsyn_str = data(targetIdx(1)).annot(contains(data(targetIdx(1)).annot,'RC_{gSYN} = '));

h = figure('visible','on');
figuresize(width, width*hwratio,h, 'inches')

for i = 1:length(best_iterations)
    
    vv = best_iterations(i);
    
    gSYNs = extractAfter(gsyn_str{vv},'RC_{gSYN} = ');
    gSYNs = str2num(gSYNs);
    
    tuningcurve(1,:) = ono.sigmoid * gSYNs(1)/0.21;
    tuningcurve(2,:) = ono.gauss * gSYNs(2)/0.21;
    tuningcurve(3,:) = ono.ushaped * gSYNs(3)/0.21;
    tuningcurve(4,:) = fliplr(ono.sigmoid) * gSYNs(4)/0.21;
    
    subplot(2,2,1); 
    plot(x,tuningcurve','b','linewidth',1);
    hold on;
    plot(x,sum(tuningcurve',2),'k','linewidth',2);
    ylim([0 max(sum(tuningcurve',2))+0.2]);
    xlim([x(1) x(end)]);
    set(gca,'xdir','reverse');
    xticks([-90,0:45:90]);
    xlabel('Azimuth');
    
    % C neuron; target or masker only cases
    perf.CT = zeros(1,4);
    fr.CT = zeros(1,4);
    if ~isempty(targetIdx)
        for i = 1:length(targetIdx)
            perf.CT(5-i) = data(targetIdx(i)).perf.C(vv);
            fr.CT(5-i) = data(targetIdx(i)).fr.C(vv);
            fr.R(:,5-i) = data(targetIdx(i)).fr.R(:,vv);
        end
    end    
    
    subplot('Position',[0.6 0.6 0.2 0.2/4])
    plotPerfGrid(perf.CT,[],[],textColorThresh);
    title('Model');       
    
    % show data grid next to model grid    
    subplot('Position',[0.6 0.5 0.2 0.2/4])
    plotPerfGrid(data_perf(1:4)',[],[],textColorThresh);
    title('Data');
    
    % calculate error and correlation with data
    [cc_clean,~] = calcModelLoss(perf.CT,data_perf(1:4)');

    str = {sprintf('Clean loss = %0.1f',loss(vv)),...
        sprintf('Clean perf C.C. = %0.3f',cc_clean),...
        sprintf('Clean avg. deviation = %0.1f',mean(abs(perf.CT-data_perf(1:4)'))),...
        data(targetIdx(1)).annot{vv,1:end}};
        
    annotation('textbox',[0.6 .35 0.2 0.1],...
           'string',str,...
           'FitBoxToText','on',...
           'LineStyle','none')
       
    % Show FR vs azimuth for clean data
    subplot(2,2,3)
    plot([-90 0 45 90],data_FR,'-b',...
        [-90 0 45 90],fliplr(fr.CT),'-r','linewidth',2);
    hold on
    plot([-90 0 45 90],ones(1,4)*mean(data_FR),'--b',...
        [-90 0 45 90],ones(1,4)*mean(fr.CT),'--r','linewidth',2);
    legend('Data','Model');
    xlabel('Azimuth');
    ylabel('Clean FR (Hz)')
    set(gca,'xdir','reverse');
    ylim([min([data_FR,fr.CT])-10 max([data_FR,fr.CT])+10]);
    xticks([-90,0:45:90]);

    % save grid
    saveas(gca,[filesep DirPart filesep 'best_iteration_V4_' num2str(vv) '.tiff'])
    clf
        
end
close;
end