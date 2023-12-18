

figure; bar(1:2,[0.025 0.035]); ylim([0 0.04]);
xlim([0.5 2.5])
set(gca,'fontsize',8);
set(gca,'xticklabel',{'w_{On}','w_{Off}'})
ylabel('g_{SYN} (\muS)')
title('Spike timing-based hotspots')
set(gcf,'unit','inches','position',[4 4 1.3 2])

figure; bar(1:2,[0.025 0.01]); ylim([0 0.04]);
xlim([0.5 2.5])
set(gca,'fontsize',8);
set(gca,'xticklabel',{'w_{On}','w_{Off}'})
ylabel('g_{SYN} (\muS)')
title('Firing rate-based hotspots')
set(gcf,'unit','inches','position',[4 4 1.3 2])