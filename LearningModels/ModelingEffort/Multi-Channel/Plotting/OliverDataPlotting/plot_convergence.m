figure;
subplot(4,2,1)
plot([18663.95,18333.21])
title('Onset->PV')
subplot(4,2,2)
plot([13182.29,12031.04])
title('PV->Ron')
subplot(4,2,3)
plot([81863.39,72514.83])
title('E->E')
subplot(4,2,4)
plot([18205.70,18103.73])
title('Offset->PV')

%%
figure;
subplot(3,1,1)
plot([17775.41,17727.77,17659.06,17608.02,17560.75,17525.22,17462.14,17393.78,17357.47,17315.16])
title('bin size = PSTH bin size = 200 == 20ms')
subplot(3,1,2)
plot([18801.21,18912.23,19054.37,19146.52,19229.80,19313.11,19366.50,19417.19,19467.67,19508.67])
title('bin size = PSTH bin size = 2000 == 200ms')
subplot(3,1,3)
plot([19206.26,19207.36,19206.58,19208.99,19209.77,19209.85,19211.69,19211.93,19211.80,19211.54])
title('bin size = PSTH bin size = 20 == 2ms')

%%
figure;
plot(squeeze(params))