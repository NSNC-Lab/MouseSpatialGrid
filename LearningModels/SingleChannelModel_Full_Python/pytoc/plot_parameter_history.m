figure;

colors = {'b','g','r','k','c'};

for k = 1:5
    plot(mean(params(:,k,:),3),colors{k}); hold on
end
legend({'On->ROn','Off->ROn','On->PV','Off->PV','PV->ROn'})
%%
figure;
spy(best_output)

%%
figure;
spy(squeeze(output(2,:,:,:)))

