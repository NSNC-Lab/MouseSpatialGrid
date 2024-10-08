%1. Load in target spatial grid (Practice grid for now)
Target_grid = [98,85,80,75;88,85,79,50;85,79,65,50;79,65,50,50;65,50,50,50];

%2. Define the number of variables (5x4 matrix flattened into a 20-element vector)
nVars = 20;

%3. Define the fitness function
fitnessFunction = @(x) fitness_fun(Target_grid, reshape(x, [5, 4]));

%4. Set GA options
options = optimoptions('ga', ...
    'PopulationSize', 8, ...
    'MaxGenerations', 50, ...
    'Display', 'iter', ...
    'CreationFcn', @create_population, ...
    'MutationFcn', @mutate_population);

%5. Run the GA
[x, fval] = ga(fitnessFunction, nVars, [], [], [], [], zeros(1, nVars), 100*ones(1, nVars), [], options);

%6. Reshape the best solution back to a 5x4 matrix
bestSolution = reshape(x, [5, 4]);
disp('Best Solution:');
disp(bestSolution);
disp(['Best Fitness: ', num2str(fval)]);

% Fitness function
function fitness = fitness_fun(Target_grid, approximate_grid)
    fitness = sum(sum(abs(Target_grid - approximate_grid)));
end

% Custom population creation function
function population = create_population(GenomeLength, ~, options)
    n = options.PopulationSize;
    population = rand(n, GenomeLength) * 100;
end

% Custom mutation function
function mutant = mutate_population(parents, options, nvars, FitnessFcn, state, thisScore, thisPopulation)
    mutationRate = 0.1;
    mutationStrength = 10;
    mutant = thisPopulation(parents,:);
    for i = 1:size(mutant, 1)
        if rand < mutationRate
            mutant(i,:) = mutant(i,:) + (rand(1, nvars) - 0.5) * mutationStrength;
        end
    end
end