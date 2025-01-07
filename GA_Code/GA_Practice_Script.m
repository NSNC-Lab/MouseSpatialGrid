%0. Set up paths, and things that will be the same across all sims.
Init_Handler;

%Changed upper left to 87 7/31 (Might allow us to get better firing rates)

%1. Load in target spatial grid (Practice grid for now)
%Hard
%Target_grid = [100,99,99,100;100,84,79,67;87,60,56,56;66,59,58,53;47,54,56,56];
%Target_fr_grid = [9.1,15.5,21.8,26.1;9.27,6.4,6.9,7.63;5.2,4.57,5.63,3.5;2.7,2.9,3.3,2.7;1.1,1.2,1.5,1.8];

%"Easier"   326 CH26
%Target_grid = [78,82,83,74;71,64,55,59;61,56,58,52;55,58,50,52;57,49,53,59];
%Target_fr_grid = [14.8,13.8,14.2,21.3;8.3,9.9,13.6,11.6;7.5,7.4,7.6,9.1;5.8,7.4,8.8,8.2;6,6.2,7.4,6];

%326 CH23
%Target_grid = [49,58,58,64;53,49,71,53;62,51,52,55;54,54,63,53;49,50,55,54];
%Target_fr_grid = [5.53,4.17,3.93,2.47;1.83,2.9,4.4,5.27;2.6,4.07,4.87,4.8;4.23,4.17,3.57,2.87];

%326 CH22
% Target_grid = [76,81,84,76;74,69,58,50;91,79,58,65;66,50,53,56;57,61,66,57];
% Target_fr_grid = [12.4,13.53,22.67,15.73;6.07,7.57,8.3,12.73;11.2,11.4,14.77,13.87;5.9,9.13,7.87,8.37;7.37,6.27,7.9,7.03];

%326 CH31
% Target_grid = [75,79,86,82;90,83,76,57;77,69,55,51;64,55,58,56;55,60,56,55];
% Target_fr_grid = [18.83,11.13,9.33,11.63;35.5,30.57,31.27,11.53;21.63,17.1,13.93,14.77;18.23,15.27,21.03,16.83;17.8,17.7,19.83,17.97];

%323 CH24
% Target_grid = [72,82,89,61;77,73,62,50;63,63,53,80;63,54,56,59;45,58,58,55];
% Target_fr_grid = [38.9,54.43,60.9,21.93;21.87,28.6,23,15.03;41.23,42.87,44.13,35.97;31.63,36.47,34.8,33.3;29.43,23,22.27,23.93];

%327 CH18
% Target_grid = [88,90,91,83;64,68,51,69;65,74,69,57;58,64,50,51;54,53,62,68];
% Target_fr_grid = [9.1,9.63,7.47,6.3;6.53,6.1,3.6,3.8;4.87,6.37,4.6,4.3;5.37,6.13,4.6,4.7;6.53,5.6,4.93,6.6];

%327 CH20
Target_grid = [77,68,80,82;83,74,73,54;67,60,56,56;55,52,53,58;50,56,51,57];
Target_fr_grid = [17.53,8.83,22.97,29.97;16.97,16.97,20.1,25.3;16.8,12.43,22.7,23.5;11.13,14,16.03,20.7;18.67,11.27,15.57,9.13];



%2. Define the number of variables
Calculate_Vars;

%3. Set GA options
optionsGA = optimoptions('ga', 'PopulationSize',50, 'MaxGenerations', 150, 'Display', 'iter', 'CreationFcn', @create_population,...
                         'MutationFcn', @mutate_populationV2, 'OutputFcn', @ga_output_function,'CrossoverFraction',0.5, ...
                         'SelectionFcn','selectiontournament','CrossoverFcn','crossovertwopoint','InitialPopulationMatrix', state.Population); 

%4. Define the fitness function
fitnessFunction = @(x) fitnessfunc_general(Target_grid,Target_fr_grid, reshape(x, [1, nVars]),optionsGA, plot_all,all_vars,XC_var_locs);

%5. Run the GA
[x, fval] = ga(fitnessFunction, nVars, [], [], [], [], [], [], [], optionsGA);

%Perhaps add functionality to adjust mutation rate and strength from here.