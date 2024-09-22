%Scatter plot
x_points = [];
y_points = [];

for k = 27:length(state.curgrid)
    %Perf grid for all points
    x_points = [x_points; reshape(state.curgrid{k},[3000,1,1])];
    %Fr frid for all poitns
    y_points = [y_points; reshape(state.curfrgrid{k},[3000,1,1])];
end

figure;
scatter(x_points,y_points)

%Take the average for all of the performance values, for what the average
%firing rate is

all_avgs = [];
points = [];
for j = 40:100
    cur_perf = [];
    for m = 1:length(x_points)
        if x_points(m) == j
            cur_perf = [cur_perf,y_points(m)];
        end
    end
    
    if length(cur_perf) > 0
        all_avgs = [all_avgs,mean(cur_perf)];
        points = [points,j];
    end

    
end

figure;
plot(points,all_avgs,'LineWidth',2)