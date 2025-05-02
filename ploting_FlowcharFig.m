close all
figure("Position",[500,100,750,750]);
subplot("Position",[0.15,0.85,0.7,0.1])
plot(song1)
title('Input Simulus','FontSize',12)
xlim([0 646000])
xticks([70000:192000:646000])
xticklabels([0:3])
yticks('')
xlabel('Time (s)')
ax = gca;
ax.FontSize = 9; 

subplot("Position",[0.58,0.47,0.33,0.27])
surf(strf.t,strf.f,strf.w1,'EdgeColor','none')
ylim([526 8100])
view(0,90)
title('STRF','FontSize',12)
ylabel('Hz')
xlabel('Time (s)')
ax = gca;
ax.FontSize = 9; 

subplot("Position",[0.12,0.47,0.33,0.27])
surf(specs.t,specs.f,transpose(specs.songs{1}),'EdgeColor','none')
ylim([526 8100])
xlim([0 3.23])
view(0,90)
title('Stimulus Spectrogram','FontSize',12)
xticks([0.35:1.92/2:3.23])
xticklabels([0:3])
ylabel('Hz')
xlabel('Time (s)')
ax = gca;
ax.FontSize = 9; 

subplot("Position",[0.12,0.08,0.33,0.27])
plot(fr_target_on{1})
xlim([0 32301])
view(0,90)
title('Onset FR','FontSize',12)
xticks([3500:(1.92/2)*10000:32301])
xticklabels([0:3])
ylabel('FR')
xlabel('Time (s)')
ax = gca;
ax.FontSize = 9; 

subplot("Position",[0.58,0.08,0.33,0.27])
plot(fr_target_off{1})
xlim([0 32301])
view(0,90)
title('Offset FR','FontSize',12)
xticks([3500:(1.92/2)*10000:32301])
xticklabels([0:3])
ylabel('FR')
xlabel('Time (s)')
ax = gca;
ax.FontSize = 9; 