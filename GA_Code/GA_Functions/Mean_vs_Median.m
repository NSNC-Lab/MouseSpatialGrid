mean_vals = [];
median_vals = [];


for j = 1:length(data_tracker)
    for k = 1:length(data_tracker{j})
            mean_vals = [mean_vals,data_tracker{j}(k).perf.C.channel1];
            median_vals = [median_vals,data_tracker{j}(k).perfmed.C.channel1];
            
    end
end

histogram(mean_vals,'FaceColor','r','FaceAlpha',0.5,'BinWidth',5); hold on;
histogram(median_vals,'FaceColor','b','FaceAlpha',0.5,'BinWidth',5); hold on;



