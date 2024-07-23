
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
    persistent bestFitness bestApproximateGrid bestvars
    if isempty(bestFitness)
        bestFitness = inf;
        bestApproximateGrid = [];
        bestvars = [];
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
        bestvars = varied_params;
    end

    % Save the variables back to the main workspace
    assignin('base', 'bestFitness', bestFitness);
    assignin('base', 'bestApproximateGrid', bestApproximateGrid);
    assignin('base', 'bestvars', bestvars);
end

% Custom population creation function
function population = create_population(GenomeLength, ~, optionsGA)
    n = optionsGA.PopulationSize;
    population = rand(n, GenomeLength)* (1/10);
end

% Custom mutation function
function mutant = mutate_population(parents, optionsGA, nvars, FitnessFcn, state, thisScore, thisPopulation)
    mutationRate = 0.1;
    mutationStrength = 0.001;
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
    bestvars = evalin('base', 'bestvars');
    
    % Display or store the best fitness and best approximate grid as needed
    % fprintf('Generation: %d, Best Fitness: %.4f\n', state.Generation, bestFitness);
    % disp('Best Approximate Grid:');
    % disp(bestApproximateGrid);

    %Store the best fitness and best approximate grid after each generation
    if ~isfield(state, 'bestFitnessHistory')
       state.bestFitnessHistory = [];
       state.bestApproximateGridHistory = {};
       state.bestvars = {};
    end
    state.bestFitnessHistory(end+1) = bestFitness;
    state.bestApproximateGridHistory{end+1} = bestApproximateGrid;
    state.bestvars{end+1} = bestvars;
    
    % Save the updated state to the base workspace
    assignin('base', 'state', state);
end
