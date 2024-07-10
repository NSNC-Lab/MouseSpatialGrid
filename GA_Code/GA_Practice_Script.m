%Add the spikingNetwork to the path
cd('..');

%Debug to remove warnings
%dbstop if warning

%1. Load in target spatial grid (Practice grid for now)
Target_grid = [98,85,80,75;80,65,63,60;65,63,60,50;65,63,60,50;65,50,50,50];

%2. Define the number of variables
nVars = 2;

%2.1 Record everything
grid_rec = {};
param_rec = {};

%3. Define the fitness function
fitnessFunction = @(x) fitness_fun(Target_grid, reshape(x, [1, nVars]),grid_rec,param_rec);

%4. Set GA options
optionsGA = optimoptions('ga', ...
    'PopulationSize', 5, ...
    'MaxGenerations', 10, ...
    'Display', 'iter', ...
    'CreationFcn', @create_population, ...
    'MutationFcn', @mutate_population, ...
    'OutputFcn', @ga_output_function);

%5. Run the GA
[x, fval] = ga(fitnessFunction, nVars, [], [], [], [], zeros(1, nVars), 100*ones(1, nVars), [], optionsGA);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %Custom Functions%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Fitness function
% function fitness = fitness_fun(Target_grid, varied_params,grid_rec,param_rec)
% 
%     %varied_struct.OnRon = varied_params(1:4);
%     %varied_struct.RonC = varied_params(5:8);
% 
%     varied_struct.IntraPV = varied_params(1);
%     varied_struct.RtoC = varied_params(2);
% 
%     SpikingNetwork_withOffset;
% 
%     %fitness = sum(sum(abs(Target_grid - approximate_grid)));
%     fitness = abs(Target_grid(2,1) - approximate_grid(2,1));
% 
% 
%     %Look at correlation coeficint
%     %Look at MSE
%     %Look at regularization
% 
% 
% end

% Fitness function
function fitness = fitness_fun(Target_grid, varied_params, grid_rec, param_rec)
    % Access the variables from the main workspace
    persistent bestFitness bestApproximateGrid
    if isempty(bestFitness)
        bestFitness = inf;
        bestApproximateGrid = [];
    end
    
    varied_struct.IntraPV = varied_params(1);
    varied_struct.RtoC = varied_params(2);
    
    % Assume this function updates approximate_grid
    SpikingNetwork_withOffset;
        
    fitness = abs(Target_grid(2,1) - approximate_grid(2,1));

    % Update the best fitness and best approximate grid if current fitness is better
    if fitness < bestFitness
        bestFitness = fitness;
        bestApproximateGrid = approximate_grid;
    end

    % Save the variables back to the main workspace
    assignin('base', 'bestFitness', bestFitness);
    assignin('base', 'bestApproximateGrid', bestApproximateGrid);
end

% Custom population creation function
function population = create_population(GenomeLength, ~, optionsGA)
    n = optionsGA.PopulationSize;
    population = rand(n, GenomeLength)* (1/10);
end

% Custom mutation function
function mutant = mutate_population(parents, optionsGA, nvars, FitnessFcn, state, thisScore, thisPopulation)
    mutationRate = 0.2;
    mutationStrength = 0.005;
    mutant = thisPopulation(parents,:);
    for i = 1:size(mutant, 1)
        if rand < mutationRate
            mutant(i,:) = mutant(i,:) + (rand(1, nvars) - 0.5) * mutationStrength;
        end
    end
end

% Custom Output function
function [state, optionsGA, optchanged] = ga_output_function(optionsGA, state, flag)
    optchanged = false;

    % Retrieve the best fitness and best approximate grid from the workspace
    bestFitness = evalin('base', 'bestFitness');
    bestApproximateGrid = evalin('base', 'bestApproximateGrid');
    
    % Display or store the best fitness and best approximate grid as needed
    % fprintf('Generation: %d, Best Fitness: %.4f\n', state.Generation, bestFitness);
    % disp('Best Approximate Grid:');
    % disp(bestApproximateGrid);

    Store the best fitness and best approximate grid after each generation
    if ~isfield(state, 'bestFitnessHistory')
       state.bestFitnessHistory = [];
       state.bestApproximateGridHistory = {};
    end
    state.bestFitnessHistory(end+1) = bestFitness;
    state.bestApproximateGridHistory{end+1} = bestApproximateGrid;
    
    % Save the updated state to the base workspace
    assignin('base', 'state', state);
end

