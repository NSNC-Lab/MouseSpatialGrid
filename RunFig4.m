%Loop through variables

track_Spike = [];

track_ISI = [];

track_RISpike = [];

track_fr = [];

for var_change1 = [0:0.1:1]

    var_change2 = 0.2;
    var_change3 = 0.5;
    
    SpikingNetwork_paper;
    track_Spike = [track_Spike, perf.SPIKE];
    track_ISI = [track_ISI, perf.ISI];
    track_RISpike = [track_RISpike, perf.RISPIKE];
    track_fr = [track_fr, data(15).fr.R2On];

    
end

for var_change2 = [0:0.1:1]

    var_change1 = 0.1;
    var_change3 = 0.5;
    
    SpikingNetwork_paper;
    track_Spike = [track_Spike, perf.SPIKE];
    track_ISI = [track_ISI, perf.ISI];
    track_RISpike = [track_RISpike, perf.RISPIKE];
    track_fr = [track_fr, data(15).fr.R2On];

    

    
end

for var_change3 = [0:0.1:1]

    var_change1 = 0.1;
    var_change2 = 0.2;
    
    SpikingNetwork_paper;
    track_Spike = [track_Spike, perf.SPIKE];
    track_ISI = [track_ISI, perf.ISI];
    track_RISpike = [track_RISpike, perf.RISPIKE];
    track_fr = [track_fr, data(15).fr.R2On];
    
    
end



figure();

subplot(2,1,1);
plot(0:0.1:1,track_Spike(1:11)); hold on
plot(0:0.1:1,track_ISI(1:11)); hold on
plot(0:0.1:1,track_RISpike(1:11)); hold on

%Find FR=30 crossover
for j = 1:11
    if track_fr(j+1).channel1 <= 30
        if track_fr(j).channel1 > 30
            x_value = (30 - track_fr(j).channel1) * 1/((-track_fr(j).channel1 - track_fr(j+1).channel1)/0.1);
            break;
        end
    end

    if track_fr(j+1).channel1 >= 30
        if track_fr(j).channel1 < 30
            x_value = (30 - track_fr(j).channel1) * 1/((track_fr(j).channel1 - track_fr(j+1).channel1)/0.1);
            break;
        end
    end
end

cross_over = (j-1)*0.1 + x_value;

%Find performance crossover
y_coord = track_Spike(j) + ((track_Spike(j+1)-track_Spike(j))/0.1)*x_value;


plot(0:0.1:1,y_coord*ones(1,11),'k--')



ylim([40 100]);
ylabel('Performance')
ytickformat('percentage');
xticks([0:0.2:1]);
set(gca, 'FontSize', 12);

subplot(2,1,2);
plot(0:0.1:1,[track_fr(1:11).channel1],'k'); hold on
plot([cross_over,cross_over],[0,80],'k--')

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
plot(0:0.1:1,track_RISpike(12:22)); 


%Find FR=30 crossover
for j = 12:22
    if track_fr(j+1).channel1 <= 30
        if track_fr(j).channel1 > 30
            x_value = (30 - track_fr(j).channel1) * 1/((-track_fr(j).channel1 - track_fr(j+1).channel1)/0.1);
            break;
        end
    end

    if track_fr(j+1).channel1 >= 30
        if track_fr(j).channel1 < 30
            x_value = (30 - track_fr(j).channel1) * 1/((track_fr(j).channel1 - track_fr(j+1).channel1)/0.1);
            break;
        end
    end
end

cross_over = (j-11-1)*0.1 + x_value;

%Find performance crossover
y_coord = track_Spike(j) + ((track_Spike(j+1)-track_Spike(j))/0.1)*x_value;


plot(0:0.1:1,y_coord*ones(1,11),'k--')

ylim([40 100]);
ylabel('Performance')
ytickformat('percentage');
set(gca, 'FontSize', 12);
xticks([0:0.2:1]);

subplot(2,1,2);
plot(0:0.1:1,[track_fr(12:22).channel1],'k')
plot([cross_over,cross_over],[0,80],'k--')

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


%Find FR=30 crossover
for j = 23:33
    if track_fr(j+1).channel1 <= 30
        if track_fr(j).channel1 > 30
            x_value = (30 - track_fr(j).channel1) * 1/((-track_fr(j).channel1 - track_fr(j+1).channel1)/0.1);
            break;
        end
    end

    if track_fr(j+1).channel1 >= 30
        if track_fr(j).channel1 < 30
            x_value = (30 - track_fr(j).channel1) * 1/((track_fr(j).channel1 - track_fr(j+1).channel1)/0.1);
            break;
        end
    end
end

cross_over = (j-22-1)*0.1 + x_value;

%Find performance crossover
y_coord = track_Spike(j) + ((track_Spike(j+1)-track_Spike(j))/0.1)*x_value;


plot(0:0.1:1,y_coord*ones(1,11),'k--')

ylim([40 100]);
ylabel('Performance')
ytickformat('percentage');

set(gca, 'FontSize', 12);
xticks([0:0.2:1]);




subplot(2,1,2);
plot(0:0.1:1,[track_fr(23:33).channel1],'k')
plot([cross_over,cross_over],[0,80],'k--')

ylim([0 80]);
ylabel('FR (HZ)')

sgtitle('PV->E')

set(gcf, 'Position', [100, 100, 350, 700]);
set(gca, 'FontSize', 12);
xticks([0:0.2:1]);

print(gcf,'-vector','-dsvg',['C:\Users\ipboy\Documents\Modeling Paper\Figures\Figure4\Figure4_PVE_ONDOM_ONOnly.svg']) % svg













%imagesc(reshape(track_Spike_perf,[51,51]))
%kernal = [0.0030,0.0133,0.0219,0.0133,0.0030;0.0133,0.0596,0.0983,0.0596,0.0133;0.0219,0.0983,0.1621,0.0983,0.0219;0.0133,0.0596,0.0983,0.0596,0.0133;0.0030,0.0133,0.0219,0.0133,0.0030];

%figure();
%imagesc(conv2(reshape(track_Spike_perf,[51,51]),kernal,'same'))