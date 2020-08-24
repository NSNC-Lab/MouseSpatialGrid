function plotFRvsgSYN(simdata,within_thresh)

temp = {simdata.name};
temp(cellfun('isempty',temp)) = {'empty'}; %label empty content

targetIdx = find(contains(temp,'m0') & ~strcmp(temp,'s0m0.mat'));
%maskerIdx = find(contains(temp,'s0') & ~strcmp(temp,'s0m0.mat'));
% mixedIdx = find(~contains(temp,'m0') & ~contains(temp,'s0') & ~contains(temp,'empty'));
model_FR = [];

for i = 1:length(targetIdx)
    model_FR(:,i) = simdata(targetIdx(i)).fr.C;
end
model_FR = mean(model_FR,2);

% get RC_gsyn values
gsyn_strs = cellfun(@str2num,extractAfter({simdata(targetIdx(1)).annot{:,2}},'RC_{gSYN} = '),'UniformOutput',false);
all_gsyns = zeros(length(gsyn_strs),1);
for i = 1:length(gsyn_strs)
    all_gsyns(i) = sum(gsyn_strs{i});
end

figure('visible','on');
scatter(all_gsyns,model_FR,'filled','b'); hold on;
scatter(all_gsyns(within_thresh),model_FR(within_thresh),'filled','r');
xlabel('sum(RC_{gSYN})');
ylabel('Mean clean FR');

% fit conductance sum vs FR line to all non-zero sets

inds = all_gsyns ~= 0;

fit = [ones(size(all_gsyns(inds))) , all_gsyns(inds)] \ model_FR(inds);
plot(all_gsyns(inds),fit(1) + fit(2)*all_gsyns(inds),'k','linewidth',1);
title('Sum of synaptic conductances vs. model firing rate');

legend('All sets','Fits within FR threshold','Fit to non-zero sets','location','best');

end