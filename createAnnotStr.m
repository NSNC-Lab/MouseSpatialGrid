function annot = createAnnotStr(data)

paramstr = {data(1).varied{2:end}};
gSYNs = []; gs = 1;
i = 1;
for aa = 1:length(data.varied)-1
    if contains(['data.' paramstr{aa}],'R_C_gSYN')
        gSYNs = cat(2,gSYNs,eval(['data.' paramstr{aa}]));
        gs = gs + 1;
    else
        annot{:,i} = sprintf('%s = %.5f',paramstr{aa},...
            eval(['data.' paramstr{aa}]));
        i = i + 1;
    end
end
% round-up gSYNs so that annotation strings don't have very long numbers
annot{:,end+1} = ['RC_{gSYN} = ' mat2str(round(10000*gSYNs)/10000)];

end