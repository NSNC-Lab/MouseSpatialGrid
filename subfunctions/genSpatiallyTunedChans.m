function [azi,spatialCurves,chanLabels,bestLocs] = genSpatiallyTunedChans(nCells)

% generate spatially-tuned channels based on # channels
 
% +90deg is contra in mouse experiments (left), while -90 is ipsi (right)

azi = -108:108;

% fits to panniello and walker 2018, already flipped so that positive is to
% the left
load('Panniello_fits.mat','spatial_fits');
chanLabels = {spatial_fits.label};

% excitatory tuning curves
if nCells == 1
    chanInds = find(contains({spatial_fits.label},'center'));
    bestLocs = 0;
    %spatialCurves(1,:) = 0.8*gaussmf(azi,[24 0]) + 0.2; % center
elseif nCells == 2
    chanInds = find(contains(chanLabels,{'contra','center'}));
    bestLocs = [90 0];
elseif nCells == 3
    chanInds = find(contains(chanLabels,{'contra','center','ipsi'}));
    bestLocs = [90 0 -90];
    %spatialCurves(1,:) = sigmf(azi,[0.1 22.5]); % contra-tuned sigmoid
    %spatialCurves(2,:) = 0.8*gaussmf(azi,[24 0]) + 0.2; % center
    %spatialCurves(3,:) = sigmf(azi,[-0.1 -22.5]); % ipsi-tuned sigmoid
elseif nCells == 4
    chanInds = find(contains(chanLabels,{'contra','45pref','center','ipsi'}));
    bestLocs = [90 45 0 -90];
end
chanLabels = chanLabels(chanInds);

for i = 1:nCells
    c_ind = find(strcmp({spatial_fits.label},chanLabels(i)));
    spatialCurves(i,:) = spatial_fits(c_ind).curve / max(spatial_fits(c_ind).curve);
end

figure;
plot(azi,spatialCurves','linewidth',1);
chanNums = cellstr(string((1:nCells)'));
set(gca,'xdir','reverse');
xlabel('stimulus location'); xtickformat('degrees');
title('channel tuning curves','fontweight','normal');
xlim([-108 108]);

hold on;
locs = [90 45 0 -90];
for i = 1:4
    line([1 1]*locs(i),[0 1],'color','k','linestyle','--','displayname','')
end
legend(strcat(chanNums,{':'},chanLabels'),'location','best');

end