timeCourse = 0.0055;
Eweight = 1;
Iweight = 1.2;
n1= 5;
n2 = 7;
x = linspace(0,0.09,18000);
y = exp(-x./timeCourse).*(((Eweight.*(x./timeCourse).^n1)./factorial(n1))-((Iweight.*(x./timeCourse).^n2)/factorial(n2)));
figure;
plot(x,y)

[song1,~] = audioread('200k_target1.wav');

song1_og = song1;

song1 = envelope(abs(song1), 200,'rms');

y2 = conv(song1,y);
y2 = y2-mean(y2);

y2_onset = [];
for k = 1:length(y2)
    if y2(k)>0
        y2_onset = [y2_onset,y2(k)];
    else
        y2_onset = [y2_onset,0];
    end
end

y2_onset(1:70000) = 0;

y2 = -y2;

y2_offset = [];
for k = 1:length(y2)
    if y2(k)>0
        y2_offset = [y2_offset,y2(k)];
    else
        y2_offset = [y2_offset,0];
    end
end

y2_offset(1:70000) = 0;

%length(y2_onset)
%y2_offset = -y2(-y2>=0);

figure;
subplot(5,1,1)
plot(song1_og)
title('original Stimulus')
subplot(5,1,2)
plot(y2_onset)
title('Temporal Convolution (onset)')
subplot(5,1,3)
plot(y2_offset)
title('Temporal Convolution (offset)')
subplot(5,1,4)
plot(fr_target_on{1,1})
title('Original (onset)')
subplot(5,1,5)
plot(fr_target_off{1,1})
title('Original (offset)')

