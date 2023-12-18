labels = {'S2On','R2On','S2Off'};

a = findobj(gcf,'type','axes');

inds = [];
for n = 1:length(a)
    if any(strcmp(labels,a(n).Title.String{1}))
    inds = cat(1,inds,n);
    end
end

k = 1; PSTHs = [];
for i = inds'
    aa = findobj(a(i),'type','line');
    PSTHs = cat(1,PSTHs,aa(1).YData);
    x = aa(1).XData;
end

figure('unit','inches','position',[4 4 1.6 1.3]);
plot(x/10,PSTHs,'linewidth',1); ylim([0 30]);
ylabel('L2/3 spike count'); xlabel('Time (ms)'); box off
xlim([400 900]);set(gca,'fontsize',8,'xtick',[400 650 900]); 