%Add a path to the GA functions
addpath('GA_Functions\')
clc;
%clear all;

%Move to the correct simulation path
cd('..');

%1. Load in target spatial grid (Practice grid for now)
Target_grid = [98,85,80,75;80,65,63,60;65,63,60,50;65,63,60,50;65,50,50,50];

%2. Define the number of variables
nVars = 3;

%3. Set GA options
%initialPopulation = state.Population;
optionsGA = optimoptions('ga', 'PopulationSize', 100, 'MaxGenerations', 30, ...
                         'Display', 'iter', 'CreationFcn', @create_population, ...
                         'MutationFcn', @mutate_population, 'OutputFcn', @ga_output_function);

                         %'InitialPopulationMatrix', initialPopulation, ...
                         % 'MutationFcn',{@mutationgaussian 1 1}, ...


%4. Define the fitness function
fitnessFunction = @(x) fitness_fun(Target_grid, reshape(x, [1, nVars]),optionsGA);

%5. Run the GA
[x, fval] = ga(fitnessFunction, nVars, [], [], [], [], zeros(1, nVars), 100*ones(1, nVars), [], optionsGA);

%Plot the evolution
%animationPlotter

