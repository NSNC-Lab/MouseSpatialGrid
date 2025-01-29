%Calculate RI-SPIKE/ISI
offset = [];
onset = [];

for j = 0.005:0.005:0.05
    for k = 0.005:0.005:0.05
        offset = [offset, j];
        onset = [onset, k];
    end
end

ratio_grid = perf(3,:)./perf(2,:);
ratio_grid_reshaped = reshape(ratio_grid,10,10);

offset_reshaped = reshape(offset,10,10);
onset_reshaped = reshape(onset,10,10);

fr_reshaped = reshape(FR,10,10);

RI_reshaped = reshape(perf(3,:),10,10);

figure;
imagesc(ratio_grid_reshaped)

%Convolve this image in order to smooth it

filter = [1,1,1;1,1,1;1,1,1];

ratio_grid_convolved = conv2(ratio_grid_reshaped,filter, 'valid')./sum(sum(filter));
%ratio_grid_convolved = ratio_grid_convolved(2:end-1,2:end-1); %Get rid of excess

figure;
imagesc(ratio_grid_convolved);
colorbar


yticklabels([0.01:0.005:0.045]')
xticklabels([0.05:0.005:0.085])

title('RI-SPIKE/ISI')
ylabel('Offset')
xlabel('Onset')


fr_convolved = conv2(fr_reshaped, filter,'valid')./sum(sum(filter));

figure;
imagesc(fr_convolved);
colorbar


yticklabels([0.01:0.005:0.045]')
xticklabels([0.05:0.005:0.085])

title('FR')
ylabel('Offset')
xlabel('Onset')


RI_convolved = conv2(RI_reshaped, filter,'valid')./sum(sum(filter));

figure;
imagesc(RI_convolved);
colorbar


yticklabels([0.01:0.005:0.045]')
xticklabels([0.05:0.005:0.085])

title('RI_{Performance}')
ylabel('Offset')
xlabel('Onset')


