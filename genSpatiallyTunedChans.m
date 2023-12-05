function [azi,spatialCurves,chanLabels] = genSpatiallyTunedChans(nCells)

% generate spatially-tuned channels based on # channels
 
% +90deg is contra in mouse experiments (left), while -90 is ipsi (right)

azi = -108:108;

% excitatory tuning curves
if nCells == 1
    spatialCurves(1,:) = 0.8*gaussmf(azi,[24 0]) + 0.2; % center

    chanLabels = {'center'};
elseif nCells == 3
    spatialCurves(1,:) = sigmf(azi,[0.1 22.5]); % contra-tuned sigmoid
    spatialCurves(2,:) = 0.8*gaussmf(azi,[24 0]) + 0.2; % center
    spatialCurves(3,:) = sigmf(azi,[-0.1 -22.5]); % ipsi-tuned sigmoid

    chanLabels = {'contra','center','ipsi'};
end

figure;
plot(azi,spatialCurves','linewidth',1);
chanNums = cellstr(string((1:3)'));
set(gca,'xdir','reverse');
xlabel('stimulus location'); xtickformat('degrees');
title('channel tuning curves','fontweight','normal');
xlim([-108 108]);

hold on;
locs = [90 45 0 -90];
for i = 1:4
    line([1 1 ]*locs(i),[0 1],'color','k','linestyle','--','displayname','')
end
legend(strcat(chanNums,{':'},chanLabels'),'location','best');

end