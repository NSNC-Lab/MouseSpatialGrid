%Select which plot
plot_num = 10;

%Plot the performance


load("on_dom_results_LOWFR.mat")

figure(Position=[600,200,800,400]);
subplot('Position',[0.1,0.3,0.5,0.65])
plot(mean(song1_holder,1),'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5); hold on
fill([1:length(song1_holder), fliplr(1:length(song1_holder))], [mean(song1_holder,1)+std(song1_holder,1), fliplr(mean(song1_holder,1)-std(song1_holder,1))],[0.8500 0.3250 0.0980], ...
    'FaceAlpha', 0.1, 'EdgeColor', 'none');

% plot(mean(song2_holder,1),'m','LineWidth',1.5); hold on
% fill([1:length(song1_holder), fliplr(1:length(song1_holder))], [mean(song2_holder,1)+std(song2_holder,1), fliplr(mean(song2_holder,1)-std(song2_holder,1))], 'red', ...
%     'FaceAlpha', 0.1, 'EdgeColor', 'none');





load("on_only_results_LOWFR.mat")

song1_holder = song1_holder_only;
% song2_holder = song2_holder_only;
plot(mean(song1_holder,1),'Color',[0.4660 0.6740 0.1880],'LineWidth',1.5); hold on
fill([1:length(song1_holder), fliplr(1:length(song1_holder))], [mean(song1_holder,1)+std(song1_holder,1), fliplr(mean(song1_holder,1)-std(song1_holder,1))], [0.4660 0.6740 0.1880], ...
    'FaceAlpha', 0.1, 'EdgeColor', 'none');
% 
% plot(mean(song2_holder,1),'g','LineWidth',1.5); hold on
% fill([1:length(song1_holder), fliplr(1:length(song1_holder))], [mean(song2_holder,1)+std(song2_holder,1), fliplr(mean(song2_holder,1)-std(song2_holder,1))], 'red', ...
%     'FaceAlpha', 0.1, 'EdgeColor', 'none');
% yticks(linspace(0,4,5))

%plot(mean(song1_holder_only,1)-mean(song1_holder,1),'b','LineWidth',1.5); hold on

% mymap = zeros(length(song1_holder_only),3);
% 
% idx = mean(song1_holder_only,1) - mean(song1_holder,1) > 0;
% 
% for j = 1:length(mymap)
%     if idx(j)>0
%         mymap(j,:) = [1,0,0];
%     else
%         mymap(j,:) = [0,0,1];
%     end
% end
% 
% colormap(mymap)
% 
% diff_data = mean(song1_holder_only,1) - mean(song1_holder,1);
% for j = 1:length(diff_data)-1
%     plot([j j+1], [diff_data(j) diff_data(j+1)], 'Color', mymap(j,:), 'LineWidth', 1.5); hold on
% end


%yticklabels(linspace(0,100,5))
ylabel('Spikes/s')
xticks('')
xlim([1 length(song1_holder)])
ax = gca;
ax.YAxis.FontSize = 14;
ax.XAxis.FontSize = 14;
ylim([0 5])


subplot('Position',[0.1,0.1,0.5,0.15])
plot(song1(starting_sample*20:ending_sample*20),'k'); hold on
xlim([0 ending_sample*20-starting_sample*20])
ylim([-1 1])
yticks('')
xticks('')
ax = gca;
ax.YAxis.FontSize = 14;
ax.XAxis.FontSize = 14;

text(-0.02, 0.5, 'Target 1', 'Units', 'normalized', ...
    'FontSize', 12, 'HorizontalAlignment', 'right');



subplot('Position',[0.7,0.1,0.25,0.85])
bar([1.37,2.17,2.96],[mean(perf_only(1,1:end)),mean(perf_only(2,1:end)),mean(perf_only(3,1:2:end))],0.4,'FaceColor','k','LineWidth',2); hold on
errorbar([1.37,2.17,2.96],[mean(perf_only(1,1:end)),mean(perf_only(2,1:end)),mean(perf_only(3,1:end))],[std(perf_only(1,1:end))/sqrt(length(perf_only(1,1:end))),std(perf_only(2,1:end))/sqrt(length(perf_only(2,1:end))),std(perf_only(3,1:end))/sqrt(length(perf_only(3,1:end)))],"LineStyle","none",'Color',[0,0,0],'LineWidth',2);

bar([1.05,1.85,2.64],[mean(perf(1,1:end)),mean(perf(2,1:end)),mean(perf(3,1:end))],0.4,'FaceColor','none','LineWidth',2); hold on
errorbar([1.05,1.85,2.64],[mean(perf(1,1:end)),mean(perf(2,1:end)),mean(perf(3,1:end))],[std(perf(1,1:end))/sqrt(length(perf(1,1:end))),std(perf(2,1:end))/sqrt(length(perf(2,1:end))),std(perf(3,1:end))/sqrt(length(perf(3,1:end)))],"LineStyle","none",'Color',[0,0,0],'LineWidth',2);
    
xlim([0.7,3.3])
ylim([50 100])
yticks([50:10:100])
ax = gca;
ax.YAxis.FontSize = 14;
ax.XAxis.FontSize = 14;
xticks([1.2,2,2.8])
xticklabels({'SPIKE','ISI','RI-SPIKE'})
ylabel('Performance')
xlabel(sprintf('Spike distance\n measure'))
box off;
print(gcf,'-vector','-dsvg',['C:\Users\ipboy\Documents\Modeling Paper\Figures\Figure',num2str(plot_num),'\Resubmission2025_2\','On_Conv_PV_COMP_HIGH','.svg']) % svg


figure(Position=[600,200,200,400]);
load("on_dom_results_LOWFR.mat")
song1_holder = song1_holder(:,16:35);


plot(mean(song1_holder,1),'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5); hold on
fill([1:size(song1_holder,2), fliplr(1:size(song1_holder,2))], [mean(song1_holder,1)+std(song1_holder,1), fliplr(mean(song1_holder,1)-std(song1_holder,1))],[0.8500 0.3250 0.0980], ...
    'FaceAlpha', 0.1, 'EdgeColor', 'none');

load("on_only_results_LOWFR.mat")
song1_holder = song1_holder_only(:,16:35);

plot(mean(song1_holder,1),'Color',[0.4660 0.6740 0.1880],'LineWidth',1.5); hold on
fill([1:size(song1_holder,2), fliplr(1:size(song1_holder,2))], [mean(song1_holder,1)+std(song1_holder,1), fliplr(mean(song1_holder,1)-std(song1_holder,1))], [0.4660 0.6740 0.1880], ...
    'FaceAlpha', 0.1, 'EdgeColor', 'none');

box off;
xticks('')
yticks('')
axis off;

print(gcf,'-vector','-dsvg',['C:\Users\ipboy\Documents\Modeling Paper\Figures\Figure',num2str(plot_num),'\Resubmission2025_2\','Advancing','.svg'])
