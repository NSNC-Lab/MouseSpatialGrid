addpath("subfunctions")
addpath("plotting")

close all

%Insert name of graph
graph_name = "tar90 mask-90 No XC";
sgtitle(graph_name);

%Insert Subgraph names
subgraphs_names = {"Best Result","Average Last Generation"};

%Insert Grid data to graph
datas = {state.bestApproximateGridHistory{end},mean(state.curgrid{end},3)};

%Insert FR data to graph
fr_datas = {state.frGridHistory{end},mean(state.curfrgrid{end},3)};

%Add 2. 1 for the target grids and 1 for the firing rate plot
num_figs = length(datas) + 2;

for k = 1:num_figs
    if k<=(num_figs-2)
        plot_performance(datas{k},fr_datas{k},subgraphs_names{k},k)
    elseif k==(num_figs-1)
        plot_performance(Target_grid,Target_fr_grid,"Target Grid",k)
    else
        plot_fr(fr_datas,Target_fr_grid,subgraphs_names,k)
    end
end
