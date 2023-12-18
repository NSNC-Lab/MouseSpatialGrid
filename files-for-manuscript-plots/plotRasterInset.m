%-- 7/13/22, 4:06 PM --%

%set(gcf,'position',[5 5 2.75 0.75])

% rate-based, control 3 laser 10

% timing-based, control 7 laser 2

inset = [2.45 2.75];
target_id = 2;

%% psths

plotsimPSTH('R2On',results,7)
plotsimPSTH('R2On',results,2)

% retrieve PSTHs from plotsimPSTH
figure(1);
subplot(2,1,target_id); a = findobj(gca,'type','bar');
ctrl = a;

figure(2);
subplot(2,1,target_id); a = findobj(gca,'type','bar');
laser = a;

% make inset figure, starting with the PSTHs
figure('unit','inches','position',[4 4 1.5 4.5]); subplot(3,1,3);
plot(ctrl.XData,smooth(ctrl.YData,3),'k','linewidth',2)
hold on;
plot(laser.XData,smooth(laser.YData,3),'r','linewidth',2)
xlim(inset); ylim([0 30]);
set(gca,'ytick',0:10:30);

%% rasters
plotsimRaster('R2On',results,7);
plotsimRaster('R2On',results,2);

% convert raster data from ms to s
figure(4); subplot(2,1,target_id); a = findobj(gca,'type','line');
ctrl = a;
ctrl.XData = ctrl.XData/1000-0.3;

% convert raster data from ms to s
figure(5); 
subplot(2,1,target_id); a = findobj(gca,'type','line');
laser = a;
laser.XData = laser.XData/1000-0.3;

%% back to inset
figure(3);
subplot(3,1,1); plot(ctrl.XData,ctrl.YData,'k'); ylim([0.5 10.5]);
xlim(inset)
subplot(3,1,2); plot(laser.XData,laser.YData,'r'); ylim([0.5 10.5]);
xlim(inset)

subplot(3,1,1); ylabel('Control');
box off
set(gca,'xtick',[],'ytick',[],'fontsize',8);
subplot(3,1,2); ylabel('Laser');
box off
set(gca,'xtick',[],'ytick',[],'fontsize',8);

subplot(3,1,3); box off; xlabel('Time (s)'); set(gca,'fontsize',8); ylabel('Spike count');
set(gcf,'position',[4 4 1.6 4.5]);
