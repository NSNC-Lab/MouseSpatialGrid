%Performance
figure;
bar([mean(perf(1,:)),mean(perf(2,:)),mean(perf(3,:))]); hold on
errorbar([std(perf(1,:)),std(perf(1,:)),std(perf(1,:))]);