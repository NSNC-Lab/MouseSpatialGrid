
load("example_cells_Isaac.mat")

figure(10)
for iTrial = 1:10

        spks = cell2mat(data(1).Spks_clean(iTrial,1,1))'
        xspikes = repmat(spks,3,1);
        yspikes = nan(size(xspikes));

        if ~isempty(yspikes)
            yspikes(1,:) = iTrial-1;
            yspikes(2,:) = iTrial;
        end
        

        plot(xspikes,yspikes,'Color','k'); hold on
end

all = [];
for iTrial = 1:10
    all = [all;cell2mat(data(1).Spks_clean(iTrial,1,1))];
end

%Resize to be between 0 and 3s. 

figure(50);
nbins = 250
h = histogram(all,nbins)

z = histcounts(all,nbins)

z2 = z(51:200)

save('90_deg_target.mat',"z2")

%[~,spk_inds] = find(raster);




    % for tid = 1:2
    %     raster = data(subz(i)).spks.C.channel1((1:10) + 10*(tid-1),:);
    %     [~,spk_inds] = find(raster);
    %     data(subz(i)).output_PSTH(tid,:) = histcounts(spk_inds,psth_vec);
    % end 