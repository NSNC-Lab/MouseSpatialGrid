hotspot = [];
for g = 1:length(grid_tracker)

    hotspot = [hotspot grid_tracker{g}(2,2)];

end

histogram(hotspot)
mean(hotspot); 
std(hotspot);


[max_val,idx] = max(hotspot);

k = 0;

%Plot maximum

figure('Position', [100, 100, 425, 700]);
plot_performance(grid_tracker{idx},fr_tracker{idx},"",k)

saveas(gcf, 'Tar-90Best.svg');

%Plot mean performance

mean_perf = zeros(5,4);
mean_fr = zeros(5,4);


for f = 1:length(grid_tracker)
    mean_perf = mean_perf + grid_tracker{f};
    mean_fr = mean_fr + fr_tracker{f};
end

mean_perf = mean_perf/length(grid_tracker);
mean_fr = mean_fr/length(grid_tracker);

figure('Position', [100, 100, 425, 700]);
plot_performance(mean_perf,mean_fr,"",k)


saveas(gcf, 'Tar-90AVG.svg');




