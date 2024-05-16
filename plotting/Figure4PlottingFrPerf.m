figure(201)

%pc.fr

plot(0:0.1:1,[pc.fr],'k')

figure(202)

plot(0:0.1:1,[pc(15).perf.SPIKE]); hold on
plot(0:0.1:1,[pc(15).perf.ISI]); hold on 
plot(0:0.1:1,[pc(15).perf.RISPIKE])
