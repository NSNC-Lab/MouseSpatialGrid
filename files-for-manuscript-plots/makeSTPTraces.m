% load('example_traces_fig2.mat')

figure('Units','inches','Position',[2 2 1.3 1]);
plot(0.1:0.1:3500,example_sim.R1On_On_PSC_syn,'k','linewidth',1); xlim([415 515]);

figure('Units','inches','Position',[2 2 1.3 1]);
plot(0.1:0.1:3500,example_sim.S1On_On_PSC_syn,'k','linewidth',1); xlim([415 515]);

figure('Units','inches','Position',[2 2 1.3 1]);
plot(0.1:0.1:3500,snn_out.R1On_S1On_PSC_syn./(snn_out.R1On_V+80),'r','linewidth',1); xlim([415 515]);

% figure('Units','inches','Position',[2 2 1.3 1]);
% plot(0.1:0.1:3500,snn_out.R1_X1_PSC_syn(:,1),'g','linewidth',1); xlim([400 500]);