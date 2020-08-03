function makeParallelPlot(data,within_thresh,loss)

temp = {data.name};
temp(cellfun('isempty',temp)) = {'empty'}; %label empty content
targetIdx = find(contains(temp,'m0') & ~strcmp(temp,'s0m0.mat'));

% get RC_gsyn values
gsyn_strs = cellfun(@str2num,extractAfter({data(targetIdx(1)).annot{:,2}},'RC_{gSYN} = '),'UniformOutput',false);
all_gsyns = zeros(length(gsyn_strs),4);
for i = 1:length(gsyn_strs)
    all_gsyns(i,:) = gsyn_strs{i};
end

% different colors for top 5 parameter sets
topcolors = [1,0,0;...  % red
    0.13,0.55,0.13;...  % forest green
    0,0,1;...   % blue
    0.8,0.8,0;...   % gold
    1,0.5,0];   % orange

x = ones(length(all_gsyns),4).*[1:4];

% scale sets by loss
ymin = min(loss(within_thresh));
ymax = max(loss(within_thresh));

loss_scale = 1;
loss_offset = zeros(length(all_gsyns),4);
loss_offset(within_thresh,:) = loss_scale*repmat((ymax-loss(within_thresh))/(ymax-ymin),1,4);

y = all_gsyns-mean(all_gsyns,2)+loss_offset;

% plot
figure('position',[200 200 600 600]);
all_lines = plot(x(within_thresh,:)',y(within_thresh,:)','k');

% replace colors of best-performing sets
[~,ind] = sort(loss(within_thresh),'ascend');
width = fliplr([1.5,1.5,2,3,4]);

for t = 1:5
    toplines(t) = all_lines(ind(t));
    toplines(t).Color = topcolors(t,:);
    toplines(t).LineWidth = width(t);
end

set(gca,'xtick',1:4);
set(gca,'xticklabels',{'RC_{ipsi}','RC_{Gauss}','RC_{U}','RC_{contra}'});
set(gca,'xdir','reverse');
xlim([0.5 4.5])
ylim([min(y,[],'all') max(y,[],'all')])
set(gca,'ytick',[0 1]);
set(gca,'yticklabels',{'Largest loss','Smallest loss'});
ytickangle(gca,45)

labels = {'Best','2^{nd}','3^{rd}','4^{th}','5^{th}'};
legend(toplines,labels,'location','best')

title('Parameter sets ordered by performance and FR loss');

end