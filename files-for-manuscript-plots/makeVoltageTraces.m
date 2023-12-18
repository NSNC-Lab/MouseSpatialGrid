figure; 

t = 0.1:0.1:3500;

fields = {'S2On','R2On','S2Off'};
titles = {'L2/3 PV_{On}','L2/3 E','L2/3 PV_{Off}'};
V_reset = [-52 -54 -52];
colors = [[200 0 0]/255;...
    0 0 0;...
    [200 0 0]/255];
styles = {'-','-','-'};

xlims = 1700;

xlims(2) = xlims(1)+250;

for s = 1:3
subplot(3,1,s); hold on;

V = example_sim(1).([fields{s} '_V']);
spk_inds = find(V(2:end) == V_reset(s) & V(1:end-1) > V_reset(s)) + 1;
V(spk_inds) = 0;

plot(t,V,'color',colors(s,:),'linewidth',1,'linestyle',styles{s});

xlim(xlims); ylim([-80 0]);
end

% scalebar

plot([xlims(1)+25 xlims(1)+25],[-50 -30],'k','linewidth',3);
plot([xlims(1)+25 xlims(1)+45],[-50 -50],'k','linewidth',3);

% paper size
for s = 1:3
subplot(3,1,s); set(gca,'fontsize',8,'xtick',[]);
end
set(gcf,'unit','inches','position',[5 5 1.5 1.75]);
savefig(gcf,'Voltage trace figure, paper');

% % presentation size
% for s = 1:3
% subplot(3,1,s); set(gca,'fontsize',18);
% title(titles{s},'fontweight','normal')
% end
% set(gcf,'unit','inches','position',[5 5 5 5.5]);
% savefig(gcf,'Voltage trace figure, presentation');

