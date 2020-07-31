function plotModelResults(ranges,DirPart,loss,model_FR,within_thresh)

% plot fr and loss functions

names = {'RC_{ipsi}','RC_{G}','RC_{U}','RC_{contra}'};
inputChans = cellfun(@(x) length(x) ~= 1,ranges);

% check how many dimensions in parameter space is used
dims = sum(inputChans);

if dims < 2
    temp = ranges{inputChans};
else
    inds = find(inputChans);
    sz1 = length(ranges{inds(1)}); sz2 = length(ranges{inds(2)});
    [X,Y] = meshgrid(ranges{inds(1)},ranges{inds(2)});
end

% plot loss vs param
figure;
if dims < 2
    title('Weights vs. Loss and Mean FR');    
    plot(temp,loss,'linewidth',1);
    hold on;
    scatter(temp(within_thresh),loss(within_thresh),80,[1 0 0],'filled');

    xlabel(names{inputChans}); ylabel('Loss');
    
    yyaxis right
    
    plot(temp,model_FR,'linewidth',1);
    hold on;
    scatter(temp(within_thresh),model_FR(within_thresh),80,[1 0 0],'filled');
    
    ylabel('Mean clean FR');
    
    legend('Loss','Loss within FR threshold','FR','FR within threshold','location','best');
    xlim([temp(1) temp(end)]);
    saveas(gcf,[filesep DirPart filesep 'MSE and FR grid search.png'])
    
elseif dims == 2
    Z = reshape(loss,sz2,sz1);
    surf(X,Y,Z,'FaceAlpha',0.3);
    
    hold on;
    scatter3(X(within_thresh),Y(within_thresh),Z(within_thresh),80,[1 0 0],'filled')
    
    xlabel(names{inds(1)}); ylabel(names{inds(2)}); zlabel('Loss');
    %zlim([0 200]);
    legend('Model','Models within threshold','location','best');
    title('Weights vs. Loss');
    saveas(gcf,[filesep DirPart filesep 'MSE grid search 2d space.fig'])
end
close;

% plot FR vs param for 2d parameter space
if dims == 2
    figure;
    Z = reshape(model_FR,sz2,sz1);
    
    surf(X,Y,Z,'FaceAlpha',0.3);
    title('Weights vs. Mean FR'); 
    hold on;
    scatter3(X(within_thresh),Y(within_thresh),Z(within_thresh),80,[1 0 0],'filled')
    xlabel(names{inds(1)}); ylabel(names{inds(2)}); zlabel('Mean clean FR');
    legend('Model','Models within threshold','location','best')
    
    saveas(gcf,[filesep DirPart filesep 'MSE FR 2d space.fig'])
    close;
end

end