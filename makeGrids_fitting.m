function makeGrids_fitting(simdata,DirPart,data_perf,data_FR,data_FR_coloc,best_iterations)

set(0,'defaultfigurevisible','on');
    
temp = {simdata.name};
temp(cellfun('isempty',temp)) = {'empty'}; %label empty content

targetIdx = find(contains(temp,'m0') & ~strcmp(temp,'s0m0.mat'));
%maskerIdx = find(contains(temp,'s0') & ~strcmp(temp,'s0m0.mat'));
mixedIdx = find(~contains(temp,'m0') & ~contains(temp,'s0') & ~contains(temp,'empty'));
textColorThresh = 70;

width=11; hwratio=0.75;
x0=.08; y0=.08;
dx=.04; dy=.04;
lx=.16; ly=.16/hwratio;

x=-108:108;
tuningcurve = zeros(4,length(x));
ono = load('ono_curves_V2.mat','sigmoid','gauss','ushaped');

gsyn_str = simdata(targetIdx(1)).annot(contains(simdata(targetIdx(1)).annot,'RC_{gSYN} = '));

h = figure('visible','on');
figuresize(width, width*hwratio,h, 'inches')

if ~exist('best_iterations','var')
    inds = 1:size(simdata(1).annot,1);
else
    inds = best_iterations;  
end

for vv = inds
            
    % CT; target only cases
    perf.CT = zeros(1,4);
    fr.CT = zeros(1,4);
    
    if ~isempty(targetIdx)
        for i = 1:length(targetIdx)
            perf.CT(i) = simdata(targetIdx(i)).perf.C(vv);
            fr.CT(i) = simdata(targetIdx(i)).fr.C(vv);
        end
    end
    
    if ~isempty(mixedIdx)
    for i = 1:length(mixedIdx)
        perf.C(i) = simdata(mixedIdx(i)).perf.C(vv);
        fr.C(i) = simdata(mixedIdx(i)).fr.C(vv);
    end
    end
    
    gSYNs = extractAfter(gsyn_str{vv},'RC_{gSYN} = ');
    gSYNs = str2num(gSYNs); %#ok<ST2NM>
    
    tuningcurve(1,:) = fliplr(ono.sigmoid) * gSYNs(1)/0.18;
    tuningcurve(2,:) = ono.ushaped * gSYNs(2)/0.18;
    tuningcurve(3,:) = ono.gauss * gSYNs(3)/0.18;
    tuningcurve(4,:) = ono.sigmoid * gSYNs(4)/0.18;
    
    % Show FR vs azimuth for clean data
    
    % flip firing rates since 90° is first index
        subplot('Position',[x0 0.35 0.5-x0 0.3-y0])
    h1=plot([-90 0 45 90],fliplr(data_FR),'-b',...
        [-90 0 45 90],fliplr(fr.CT),'-r','linewidth',2);
    hold on
    plot([-90 0 45 90],ones(1,4)*mean(data_FR),'--b',...
        [-90 0 45 90],ones(1,4)*mean(fr.CT),'--r','linewidth',2);
    ylabel('Clean FR (Hz)')
    set(gca,'xdir','reverse');
    ylim([min([data_FR,fr.CT])-10 max([data_FR,fr.CT])+10]);
    xticks([-90,0:45:90]);
    legend([h1(1),h1(2)],'Data','Model');
    
    % for co-located data
        subplot('Position',[x0 y0 0.5-x0 0.3-y0])
    h2=plot([-90 0 45 90],fliplr(data_FR_coloc),'-b',...
        [-90 0 45 90],fliplr(fr.C),'-r','linewidth',2);
    hold on; 
    plot([-90 0 45 90],ones(1,4)*mean(data_FR_coloc),'--b',...
        [-90 0 45 90],ones(1,4)*mean(fr.C),'--r','linewidth',2);
    ylabel('Co-located FR (Hz)')
    set(gca,'xdir','reverse');
    ylim([min([data_FR_coloc,fr.C([1:5:end])])-10 max([data_FR_coloc,fr.C([1:5:end])])+10]);
    xticks([-90,0:45:90]);
    xlabel('Azimuth');
    legend([h2(1),h2(2)],'Data','Model');
    
    % make subplot of tuning curves
    subplot('Position',[x0 0.7-y0 0.5-x0 0.3-y0]);
    plot(x,tuningcurve','b','linewidth',1);
    hold on;
    plot(x,sum(tuningcurve),'k','linewidth',2);
    ylim([0 max(sum(tuningcurve))+0.2]);
    xlim([x(1) x(end)]);
    set(gca,'xdir','reverse');
    xticks([-90,0:45:90]);
    
    % plot model grid clean + co-located
    subplot('Position',[0.57 0.8+dy/4 lx ly/2])
    plotPerfGrid([perf.CT;perf.C],[],textColorThresh);
    title('Model');
    
    % show data grid next to model grid clean
    subplot('Position',[0.8 0.8+dy/4 lx ly/2])
    plotPerfGrid([data_perf(1:4)';data_perf(5:5:end)'],[],textColorThresh);
    title('Data');
    
    % plot difference grid clean
    subplot('Position',[0.57 0.4+dy/4 lx ly/2])
    plotPerfGrid([data_perf(1:4)'-perf.CT;data_perf([5:5:20])'-perf.C],[],-5);
    title('Difference');
    
    % calculate error and correlation with data
    [cc,loss] = calcModelLoss([perf.CT,perf.C],[data_perf(1:4)',data_perf([5:5:20])']);
    perf_str = sprintf('Clean+coloc perf C.C. = %0.3f',cc);
    dev_str = sprintf('Clean+coloc perf deviation = %0.1f',mean(abs([perf.CT,perf.C]-[data_perf(1:4)',data_perf([5:5:20])'])));
    
    % round gsyns for plotting purposes
    gSYNs_rounded = round(gSYNs*10000)/10000;
    gsyn_annot = ['RC_{gSYN} = ' mat2str(gSYNs_rounded)];
    
    str = {sprintf('Loss = %0.1f',loss(1)),...
        perf_str,...
        dev_str,...
        simdata(targetIdx(1)).annot{vv,1},...
        simdata(mixedIdx(1)).annot{vv,1},...
        gsyn_annot,...
        simdata(targetIdx(1)).annot{vv,end}};
    
    annotation('textbox',[0.75 0.35 0.2 0.1],...
        'string',str,...
        'FitBoxToText','on',...
        'LineStyle','none')

    % save grid
    saveas(gca,[filesep DirPart filesep 'best_iteration_' num2str(vv) '.tiff'])
     clf;  
end
%close;

end