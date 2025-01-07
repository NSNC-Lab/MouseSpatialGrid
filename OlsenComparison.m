figure;
subplot(3,1,1);
load('on_dom.mat')
plot(avg_data); hold on;
% plot([200,200],[0,40],'k--'); hold on
% plot([350,350],[0,40],'k--'); hold on
% plot([500,500],[0,40],'k--'); hold on
% plot([650,650],[0,40],'k--'); hold on
% plot([800,800],[0,40],'k--'); hold on
% plot([950,950],[0,40],'k--'); hold on
title('ON')
xlim([0 1819])
subplot(3,1,2);
load('off_dom.mat')
plot(avg_data)
title('OFF')
xlim([0 1819])
subplot(3,1,3);
load('both_dom.mat')
plot(avg_data)
title('BOTH')
xlim([0 1819])


figure;
load('on_dom.mat')
plot(avg_data,'LineWidth',2); hold on
load('off_dom.mat')
plot(avg_data,'LineWidth',2)

xlim([1200 1819])