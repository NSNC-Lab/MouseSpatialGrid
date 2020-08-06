function makeGrids_bestIteration(data,varies,DirPart,data_perf,data_FR,best_iterations,loss)

set(0,'defaultfigurevisible','on');
    
% performance vector has dimensions [numSpatialChan,nvaried]
neurons = {'Ipsi. sigmoid','Gaussian','U','Cont. sigmoid'};

temp = {data.name};
temp(cellfun('isempty',temp)) = {'empty'}; %label empty content

targetIdx = find(contains(temp,'m0') & ~strcmp(temp,'s0m0.mat'));
%maskerIdx = find(contains(temp,'s0') & ~strcmp(temp,'s0m0.mat'));
mixedIdx = find(~contains(temp,'m0') & ~contains(temp,'s0') & ~contains(temp,'empty'));
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

width=11; hwratio=0.75;
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
    
        % CT; target only cases
    perf.CT = zeros(1,4);
    fr.CT = zeros(1,4);
    
    if ~isempty(targetIdx)
        for i = 1:length(targetIdx)
            perf.CT(i) = data(targetIdx(i)).perf.C(vv);
            fr.CT(i) = data(targetIdx(i)).fr.C(vv);
        end
    end
    
    if ~isempty(mixedIdx)
    for i = 1:length(mixedIdx)
        perf.C(i) = data(mixedIdx(i)).perf.C(vv);
        fr.C(i) = data(mixedIdx(i)).fr.C(vv);
    end
    end
    
    gSYNs = extractAfter(gsyn_str{vv},'RC_{gSYN} = ');
    gSYNs = str2num(gSYNs);
    
    tuningcurve(1,:) = fliplr(ono.sigmoid) * gSYNs(1)/0.21;
    tuningcurve(2,:) = ono.ushaped * gSYNs(2)/0.21;
    tuningcurve(3,:) = ono.gauss * gSYNs(3)/0.21;
    tuningcurve(4,:) = ono.sigmoid * gSYNs(4)/0.21;
    
    % Show FR vs azimuth for clean data
    subplot('Position',[x0 y0 0.5-x0 0.4-y0])
    plot([-90 0 45 90],fliplr(data_FR),'-b',...
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
    
    subplot('Position',[x0 0.6-y0 0.5-x0 0.4-y0]); 
    plot(x,tuningcurve','b','linewidth',1);
    hold on;
    plot(x,sum(tuningcurve',2),'k','linewidth',2);
    ylim([0 max(sum(tuningcurve',2))+0.2]);
    xlim([x(1) x(end)]);
    set(gca,'xdir','reverse');
    xticks([-90,0:45:90]);
    xlabel('Azimuth');
   
    
    % plot model grid
    subplot('Position',[0.57 0.75+dy/4 dx*4 dy])
    plotPerfGrid(perf.CT,[],textColorThresh);
    title('Model');       
    
    subplot('Position',[0.57 0.75-dy*4 dx*4 dy*4]);
    plotPerfGrid(perf.C,[],textColorThresh);
    
    % show data grid next to model grid    
    subplot('Position',[0.8 0.75+dy/4 dx*4 dy])
    plotPerfGrid(data_perf(1:4)',[],textColorThresh);
    title('Data');
    
    subplot('Position',[0.8 0.75-dy*4 dx*4 dy*4]);
    plotPerfGrid(data_perf(5:end)',[],textColorThresh);
    
    % plot difference grid
    subplot('Position',[0.57 0.4+dy/4 dx*4 dy])
    plotPerfGrid(perf.CT-data_perf(1:4)',[],-5);
    title('Difference');
    
    subplot('Position',[0.57 0.4-dy*4 dx*4 dy*4]);
    plotPerfGrid(perf.C-data_perf(5:end)',[],-5);

    % calculate error and correlation with data
    [cc_clean,~] = calcModelLoss([perf.CT,perf.C],data_perf');

    str = {sprintf('Loss = %0.1f',loss(vv)),...
        sprintf('Perf C.C. = %0.3f',cc_clean),...
        sprintf('Clean perf deviation = %0.1f',mean(abs([perf.CT,perf.C]-data_perf'))),...
        data(targetIdx(1)).annot{vv,1:end}};
        
    annotation('textbox',[0.8 .35 0.2 0.1],...
           'string',str,...
           'FitBoxToText','on',...
           'LineStyle','none')

    % save grid
    saveas(gca,[filesep DirPart filesep 'best_iteration_' num2str(vv) '.tiff'])
    clf
        
end
close;
end