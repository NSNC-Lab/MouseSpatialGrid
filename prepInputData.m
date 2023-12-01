
padToTime = 3500; %ms

% z : 1-4 , masker only
% z : 5,10,15,20 , target only
% z : 6-9 and etc. , mixed trials

% locs are based on the locations for stimuli in mouse experiments
% +90, 45, and 0 are to the left (contra)
% -90 to the right (ipsi)
locs = [90 45 0 -90];

temp = cellfun(@numel,fr_target_on);
tmax = max(temp);

x = -108:108;

% excitatory tuning curves
if options.nCells == 1
    tuningcurve(1,:) = 0.8*gaussmf(x,[24 0]) + 0.2; % center
elseif options.nCells == 3
    tuningcurve(1,:) = sigmf(x,[0.1 22.5]); % contra-tuned sigmoid
    tuningcurve(2,:) = 0.8*gaussmf(x,[24 0]) + 0.2; % center
    tuningcurve(3,:) = sigmf(x,[-0.1 -22.5]); % ipsi-tuned sigmoid
end

if ~isempty(options.locNum)
    subz = 1:24;
else
    subz = options.locNum;
end

% only create spks file if not done yet

% if ~exist(fullfile(study_dir,'solve','IC_spks_on.mat'),'file') %|| options.regenSpks
    labels = {'on','off'};
    for ICtype = [1 2]
        % divide all times by dt to upsample the time axis
        spks = [];

        % load fr traces, scale by weights
        for z = subz

            % target and masker weights
            if z <= 4  % masker only
                tloc(z) = nan;
                mloc(z) = locs(z);
            elseif mod(z,5) == 0 % target only
                tloc(z) = locs(floor(z/5));
                mloc(z) = nan;
            else % mixed
                tloc(z) = locs(floor(z/5));
                mloc(z) = locs(mod(z,5));
            end

            singleConfigSpks = zeros(20,1,tmax);

            for t = 1:20 % trials [1:10]
                for ch = 1:options.nCells % channels [(1:4),(1:4)]

                    t_wt = tuningcurve(ch,x == tloc(z));
                    m_wt = tuningcurve(ch,x == mloc(z));

                    if isempty(t_wt), t_wt = 0; end
                    if isempty(m_wt), m_wt = 0; end

                    if t <= 10 %song 1
                        singleConfigSpks(t,ch,:) = t_wt.*eval(['fr_target_' labels{ICtype} '{1}']) + m_wt.*fr_masker{t};
                    else
                        singleConfigSpks(t,ch,:) = t_wt.*eval(['fr_target_' labels{ICtype} '{2}']) + m_wt.*fr_masker{t-10};
                    end

                    if t_wt + m_wt >= 1
                        singleConfigSpks(t,ch,:) = singleConfigSpks(t,ch,:) / (t_wt + m_wt);
                    end
                end
            end

            % format of spks is : [trial x channel x time]
            % pad each trial to have duration of timePerTrial
            if size(singleConfigSpks,3) < padToTime/dt
                padSize = padToTime/dt-size(singleConfigSpks,3);
                singleConfigSpks = cat(3,singleConfigSpks,zeros(20,options.nCells,padSize));
            end

            % increase or decrease gain factor of strf
            % spks = cat(3,spks,singleConfigSpks(:,1,:)) * newStrfGain / strfGain;
            spks = cat(3,spks,singleConfigSpks) * newStrfGain / strfGain;
            % strfGain is from default_STRF_with_offset.mat == 0.1

        end
        % format of spks should be : [time x channel x trial]
        spks = permute(spks,[3 2 1]);
        save(fullfile(study_dir, 'solve',['IC_spks_' labels{ICtype} '.mat']),'spks','dt');
    end
% end
