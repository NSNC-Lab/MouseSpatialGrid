%Get the mean fitness
mean_fitnesses = [];
for k = 1:length(state.curfitness)
    mean_fitnesses = [mean_fitnesses mean(state.curfitness{k})]
end

%Get the convergences
convergences = state.convergence;

%Plot them
figure;
plot(mean_fitnesses); hold on
plot(convergences);

legend({'mean fitness','convergence'})