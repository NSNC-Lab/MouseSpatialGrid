function plotFiringGrid(neuronPerf,neuronFiring,tit,textColorThresh)
% plots performance grids; can handle a data array of 16 or 8 elements
% works with mouseNetwork_main

if numel(neuronPerf) == 16
    % 4x4 grid
    [X,Y] = meshgrid(1:4,4:-1:1);
    str = cellstr(num2str(round(neuronFiring(:))));
    neuronPerf = reshape(neuronPerf,4,4);
    for i = 1:length(str)
       str{i} = ['(',str{i},')']; 
    end
    imagesc(flipud(neuronPerf));
    xticks([1:4]); xticklabels({'90','45','0','-90'})
    yticks([1:4]); yticklabels(fliplr({'90','45','0','-90'}))
    xlabel('Song Location')
    ylabel('Masker Location')
elseif numel(neuronPerf) == 8
    % 2x4 grid
    [X,Y] = meshgrid(1:4,1:2);
    str = cellstr(num2str(round(neuronFiring(:))));
    for i = 1:length(str)
       str{i} = ['(',str{i},')']; 
    end
    imagesc(neuronPerf);
    xticks([]);
    yticks(1:2); yticklabels({'target only','masker only'})
end

title(tit)
colormap('parula');
t = text(X(:),Y(:),str,'Fontsize',13,'HorizontalAlignment', 'Center');
for i= 1:numel(neuronPerf)
    if neuronPerf(i)>textColorThresh
        t(i).Color = 'k';
    else
        t(i).Color= 'w';
    end
end
caxis([50,100])
set(gca,'fontsize',12)

