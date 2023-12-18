s = findobj(gcf,'type','axes'); % subplots

zoom_in = [];

for n = 1:length(s)

    axes(s(n));

    b = get(s(n),'children');

    b.YData = smooth(b.YData,3);
    if ismember(n,[3 5 6])
        zoom_in = cat(1,zoom_in,b.YData);
    end

end

% make zoomed in PSTH for timing
figure('unit','inches','position',[ 4 4 1.6 1.3]);
t_inds = find([0:200:35000] >= 3600 & [0:200:35000] <= 8600);
colors = {'g','r','k'};

for n = 1:3
hold on; plot(360:20:860,zoom_in(n,t_inds),colors{n},'linewidth',1);
end
xlim([360 860]); xlabel('Time (ms)');
ylim([0 40]); ylabel('L4 spike count');
legend('PV_{Off}','PV_{On}','E'); set(gca,'fontsize',8)