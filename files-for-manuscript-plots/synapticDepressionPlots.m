x = 0:0.1:1;

figure('unit','inches','position',[ 4 4 2.25 3.5]);
subplot(2,1,1);
plot(x,pc.SPIKE); hold on; plot(x,pc.ISI,x,pc.RISPIKE); ylim([50 100])
legend('SPIKE','ISI','RI-SPIKE'); set(gca,'fontsize',8,'xticklabel',[],'xtick',0:0.2:1)
ylabel('Performance'); ytickformat('percentage')

subplot(2,1,2);
plot(x,fr,'k');
xlabel('f_{P}, E->E'); 
sgtitle('E->E synaptic depression strength');
ylabel('FR (Hz)'); ylim([0 80]);
set(gca,'fontsize',8,'xtick',0:0.2:1);
