function makeParallelPlot(data,within_thresh,loss)

temp = {data.name};
temp(cellfun('isempty',temp)) = {'empty'}; %label empty content

targetIdx = find(contains(temp,'m0') & ~strcmp(temp,'s0m0.mat'));

gsyn_strs = cellfun(@str2num,extractAfter({data(targetIdx(1)).annot{:,2}},'RC_{gSYN} = '),'UniformOutput',false);
all_gsyns = zeros(length(gsyn_strs),4);
for i = 1:length(gsyn_strs)
    all_gsyns(i,:) = gsyn_strs{i};
end

% group cases by 1) those outside of threshold and then 2) by loss value
grpvars = ones(size(gsyn_strs));    % 1: outside of FR threshold
grpvars(within_thresh & loss > 150) = 2;
grpvars(within_thresh & loss > 100 & loss < 150) = 3;
grpvars(within_thresh & loss > 80 & loss < 100) = 4;
grpvars(within_thresh & loss > 60 & loss < 80) = 5;
grpvars(within_thresh & loss < 60) = 6;

% different colors for top 5 parameter sets
topcolors = [1,0,0;...  % red
    0.13,0.55,0.13;...  % forest green
    0,0,1;...   % blue
    0.8,0.8,0;...   % gold
    1,0.5,0];   % orange

temp = cool(max(grpvars)-1);

colors = {[0,0,0]};
for i = 1:max(grpvars)-1
    colors{i+1} = temp(i,:);
end
linewidths = [0.01 0.1 0.5 1 2 4];
linetype = {':','-','-','-','-','-'};

jitter_scale = 0.004;
jitter = jitter_scale*rand(length(all_gsyns),4)-jitter_scale/2;

x = ones(length(all_gsyns),4).*[1:4];

ymin = min(loss(within_thresh));
ymax = max(loss(within_thresh));

loss_scale = 0.5;
loss_offset = zeros(length(all_gsyns),4);
loss_offset(within_thresh,:) = loss_scale*repmat((ymax-loss(within_thresh))/(ymax-ymin),1,4);

y = all_gsyns-mean(all_gsyns,2)+jitter+loss_offset;

% plot sets
figure('position',[200 200 400 800]);
all_lines = plot(x',y');

% replace sets by colors
for i = 1:max(grpvars)
    hold on;
    temp = all_lines(grpvars==i);
    numlines(i) = sum(grpvars == i);
    for li = 1:length(temp)
       temp(li).LineWidth = linewidths(i);
       temp(li).Color = colors{i};
       temp(li).LineStyle = linetype{i};
    end
    hlines(i) = temp(1);
end

% replace colors of best-performing sets
temp2 = sort(loss(within_thresh),'ascend');
for t = 1:5
    toplines(t) = all_lines(loss == temp2(t));
    toplines(t).Color = topcolors(t,:);
end

set(gca,'xtick',1:4);
set(gca,'xticklabels',{'RC_{ipsi}','RC_{Gauss}','RC_{U}','RC_{contra}'});
set(gca,'xdir','reverse');
xlim([0.5 4.5])
ylim([min(y,[],'all') max(y,[],'all')])
ylabel('RC_{gSYN}');

% set up legend

labels = {'outside FR','loss > 150','loss [100,150]','loss [80,100]',...
    'loss [60,80]','loss < 60'};

numvars = 0;
t = 0;
while numvars < 5
    numvars = numvars + numlines(end-t);
    t = t+1;
end

    hlines(end-t+1:end) = [];
    labels(end-t+1:end) = [];
    hlines(end+1:end+5) = toplines;
    labels(end+1:end+5) = {'best','2nd best','3rd best','4th best','5th best'};


legend(hlines,labels,'location','best')

end