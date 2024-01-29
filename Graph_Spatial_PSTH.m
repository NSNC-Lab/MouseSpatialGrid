close all;

%In data grab the things we want
names = {'On','Off','SOnOff','ROn'};
chans = {'channel1','channel2','channel3','channel4'};
orientation = [2,5,7,9];

for m = 1:length(locs)
    for k = 1:length(locs)
        figure(k+(m-1)*4)
        for j = 1:length(names)
            
            subplot(3,3,orientation(j))
        
            all = 0;
            nbins = 100;


            for i = 1:20
                %Important (There will be no response on other channels because this is just the )
                spks = find(data(subz(m)).spks.(names{j}).(chans{k})(i,:)==1)*3/7000/7000;
                xspikes = repmat(spks,3,1);
                yspikes = nan(size(xspikes));
                yspikes(1,:) = i-1;
                yspikes(2,:) = i;
            
                plot(xspikes,yspikes,'Color','k'); hold on

                all = all + data(subz(m)).spks.(names{j}).(chans{k})(i,:);
                ylim([0 4])
            end

            h = plot(movmean(all,100));
        end
    end
end




%Add C for each target location
%Convert to PSTHs