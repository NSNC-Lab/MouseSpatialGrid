close all;
figure();

subplot(2,1,1);
plot(0:0.1:1,track_Spike(1:11)); hold on
plot(0:0.1:1,track_ISI(1:11)); hold on
plot(0:0.1:1,track_RISpike(1:11)); hold on

%Find FR=30 crossover


x_value = 0.1;
y1_value = track_fr(2).channel1;


%Find performance crossover
y2_value = track_Spike(2);
%y_coord = track_Spike(j) + ((track_Spike(j+1)-track_Spike(j))/0.1)*x_value;


plot(0:0.1:1,y2_value*ones(1,11),'k--'); hold on
plot([x_value,x_value],[40, 100],'k--')



ylim([40 100]);
ylabel('Performance')
ytickformat('percentage');
xticks([0:0.2:1]);
set(gca, 'FontSize', 12);

subplot(2,1,2);
plot(0:0.1:1,[track_fr(1:11).channel1],'k'); hold on
plot([x_value,x_value],[0,80],'k--'); hold on
plot([0,1],[y1_value,y1_value],'k--'); hold on


ylim([0 80]);
ylabel('FR (HZ)')

sgtitle('E->E')

set(gcf, 'Position', [100, 100, 350, 700]);
set(gca, 'FontSize', 12);
xticks([0:0.2:1]);

print(gcf,'-vector','-dsvg',['C:\Users\ipboy\Documents\Modeling Paper\Figures\Figure4\Figure4_EE_ONDOM_ONOnly.svg']) % svg


figure();

subplot(2,1,1);
plot(0:0.1:1,track_Spike(12:22)); hold on
plot(0:0.1:1,track_ISI(12:22)); hold on
plot(0:0.1:1,track_RISpike(12:22)); hold on
x_value = 0.2;
y1_value = track_fr(14).channel1;


%Find performance crossover
y2_value = track_Spike(14);
%y_coord = track_Spike(j) + ((track_Spike(j+1)-track_Spike(j))/0.1)*x_value;


plot(0:0.1:1,y2_value*ones(1,11),'k--'); hold on
plot([x_value,x_value],[40, 100],'k--')

ylim([40 100]);
ylabel('Performance')
ytickformat('percentage');
set(gca, 'FontSize', 12);
xticks([0:0.2:1]);

subplot(2,1,2);
plot(0:0.1:1,[track_fr(12:22).channel1],'k'); hold on
plot([x_value,x_value],[0,80],'k--'); hold on
plot([0,1],[y1_value,y1_value],'k--'); hold on

ylim([0 80]);
ylabel('FR (HZ)')

sgtitle('E->PV')
set(gcf, 'Position', [100, 100, 350, 700]);
set(gca, 'FontSize', 12);
xticks([0:0.2:1]);

print(gcf,'-vector','-dsvg',['C:\Users\ipboy\Documents\Modeling Paper\Figures\Figure4\Figure4_EPV_ONDOM_ONOnly.svg']) % svg


figure();
subplot(2,1,1);
plot(0:0.1:1,track_Spike(23:33)); hold on
plot(0:0.1:1,track_ISI(23:33)); hold on
plot(0:0.1:1,track_RISpike(23:33)); 


x_value = 0.5;
y1_value = track_fr(28).channel1;


%Find performance crossover
y2_value = track_Spike(28);
%y_coord = track_Spike(j) + ((track_Spike(j+1)-track_Spike(j))/0.1)*x_value;


plot(0:0.1:1,y2_value*ones(1,11),'k--'); hold on
plot([x_value,x_value],[40, 100],'k--')


ylim([40 100]);
ylabel('Performance')
ytickformat('percentage');

set(gca, 'FontSize', 12);
xticks([0:0.2:1]);




subplot(2,1,2);
plot(0:0.1:1,[track_fr(23:33).channel1],'k'); hold on
plot([x_value,x_value],[0,80],'k--'); hold on
plot([0,1],[y1_value,y1_value],'k--'); hold on

ylim([0 80]);
ylabel('FR (HZ)')

sgtitle('PV->E')

set(gcf, 'Position', [100, 100, 350, 700]);
set(gca, 'FontSize', 12);
xticks([0:0.2:1]);

print(gcf,'-vector','-dsvg',['C:\Users\ipboy\Documents\Modeling Paper\Figures\Figure4\Figure4_PVE_ONDOM_ONOnly.svg']) % svg