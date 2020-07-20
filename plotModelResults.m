function plotModelResults(weights,DirPart,loss,frHistory,within_thresh,VRdiff)

% plot fr and loss functions

names = {'RC_{ipsi}','RC_{G}','RC_{U}','RC_{contra}'};
inputChans = find(sum(weights,2) ~= 0);

temp = weights(inputChans,:);

% check how many dimensions in parameter space is used
dims = 0;
for r = 1:length(inputChans)
    dims = dims + double(length(unique(temp(r,:))) ~= 1);
    if length(unique(temp(r,:))) == 1
        inputChans(r) = [];
        temp(r,:) = [];
    end
end

% plot loss vs param
figure;
if dims < 2
    title('Weights vs. Loss and Mean FR');
    [temp,I] = sort(temp);
    
    plot(temp,loss(I),'linewidth',1);
    hold on;
    scatter(temp(within_thresh),loss(within_thresh),80,[1 0 0],'filled');

    xlabel(names{inputChans}); ylabel('Loss');
    
    yyaxis right
    
    plot(temp,frHistory(I),'--k','linewidth',1);
    hold on;
    scatter(temp(within_thresh),frHistory(within_thresh),80,[1 0 0],'filled');
    
    ylabel('Mean clean FR');
    
    legend('Loss','Loss within FR threshold','FR','FR within threshold','location','best');
    xlim([temp(1) temp(end)]);
    saveas(gcf,[filesep DirPart filesep 'MSE and FR grid search.png'])
    
elseif dims == 2
    
    sz1 = length(unique(temp(1,:))); sz2 = length(unique(temp(2,:)));
    
    X = reshape(temp(1,:)',sz2,sz1); Y = reshape(temp(2,:)',sz2,sz1);
    Z = reshape(loss,sz2,sz1);
    surf(X,Y,Z);
    
    hold on;
    scatter3(X(within_thresh),Y(within_thresh),Z(within_thresh),80,[1 0 0],'filled')
    
    xlabel(names{inputChans(1)}); ylabel(names{inputChans(2)}); zlabel('Loss');
    zlim([0 200]);
    legend('Model','Models within threshold','location','best');
    title('Weights vs. Loss');
    saveas(gcf,[filesep DirPart filesep 'MSE grid search 2d space.fig'])
end

close;

% plot FR vs param for 2d parameter space
if dims == 2
    figure;
       
    Z = reshape(frHistory,sz2,sz1);
    
    surf(X,Y,Z);
    title('Weights vs. Mean FR'); 
    hold on;
    scatter3(X(within_thresh),Y(within_thresh),Z(within_thresh),80,[1 0 0],'filled')
    xlabel(names{inputChans(1)}); ylabel(names{inputChans(2)}); zlabel('Mean clean FR');
    legend('Model','Models within threshold','location','best')
    
    saveas(gcf,[filesep DirPart filesep 'MSE FR 2d space.fig'])
end


%% plot weights vs VR MSE
figure;
if dims < 2
    plot(temp,VRdiff(I),'linewidth',1);
    xlabel(names{inputChans}); ylabel('Sq. VR difference');
    hold on;
    scatter(temp(within_thresh),VRdiff(within_thresh),80,[1 0 0],'filled');
    title('Weights vs. Mean VR^{2} Difference');
    legend('Model','Models within threshold','location','best');
    xlim([temp(1) temp(end)]);
    saveas(gcf,[filesep DirPart filesep 'MSE and VR grid search.png'])
elseif dims == 2    
    Z = reshape(VRdiff,sz2,sz1);
    surf(X,Y,Z);
    
    hold on;
    scatter3(X(within_thresh),Y(within_thresh),Z(within_thresh),80,[1 0 0],'filled')
    xlabel(names{inputChans(1)}); ylabel(names{inputChans(2)}); zlabel('Sq. VR difference');
    title('Weights vs. Mean VR^{2} Difference');
    legend('Model','Models within threshold','location','best')
    
    saveas(gcf,[filesep DirPart filesep 'MSE and VR grid search 2d space.fig'])
end

end