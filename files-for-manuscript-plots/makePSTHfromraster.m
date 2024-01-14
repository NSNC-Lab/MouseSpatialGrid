% Use this function to overlay PSTHs on top of rasters (see Figures 4 and 5
% in the bioRxiv manuscript)

% You'll need to have one of the figures make by plotRasterTree open

s = findobj(gcf,'type','axes'); % subplots

zoom_in = [];

for n = 1:length(s)
    
    axes(s(n));
    
    b = get(s(n),'children');
    for m = 1:length(b)
        
        if numel(b(m).XData) == 2
            delete(b(m));
            continue;
        end
        
        % use only target 1 responses (10.5 - 20.5 in YData)
        trialinds = b(m).YData(2:3:end);
        t1times = b(m).XData(2:3:end); %
        t1times = t1times(trialinds >= 11.5);
        
        temp = find(trialinds >= 11.5);
        
        b(m).YData = b(m).YData(sort([(temp-1)*3+1,(temp-1)*3+2,(temp-1)*3+3]))-10;
        b(m).XData = b(m).XData(sort([(temp-1)*3+1,(temp-1)*3+2,(temp-1)*3+3]));
        b(m).LineWidth = 0.5; b(m).Color = [0.2 0.2 0.2];
        set(gca,'ydir','normal'); ylim([0.5 10.5]);
        psth = smooth(histcounts(t1times,0:200:35000),3);
        psth(end+1) = 0; hold on;

        if ismember(n,[5 6 8])
        zoom_in = cat(1,zoom_in,psth');
        end
        yyaxis right;
        plot([0:200:35000],psth,'linewidth',1,'color','k'); ylim([0 30]);
        
        xlim([3600 13600]);
        
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




