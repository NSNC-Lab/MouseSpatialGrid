%Todos 
%Get error bars (1std) for rate modulation graph using current simulation
figure
plot_data = [];
errorbs = [];

for h = 1:5

   norm = [RM(:,h).R2On]/mean([RM(:,1).R2On]);
   plot_data = [plot_data mean(norm)];
   errorbs = [errorbs std(norm)];
end

plot(plot_data,'k','LineWidth',0.2); hold on
errorbar(plot_data,errorbs,'k','LineWidth',0.2)
xlim([0.9,5.1])

