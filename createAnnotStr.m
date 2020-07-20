function annotstr = createAnnotStr(data,varies,perf)

paramstr = {data(1).varied{2:end}};
annotstr{:,1} = ['Disc = ' num2str(mean(max(perf.C(vv))))];

RCs = []; rc = 1;
gSYNs = []; gs = 1;
i = 1;
for aa = 1:length(varies)-1
    if contains(['data.' paramstr{aa}],'R_C_inputChan')
        RCs = cat(2,RCs,eval(['data.' paramstr{aa}]));
        rc = rc + 1;
    elseif contains(['data.' paramstr{aa}],'R_C_gSYN')
        gSYNs = cat(2,gSYNstr,eval(['data.' paramstr{aa}]));
        gs = gs + 1;
    else
        annotstr{:,i+1} = sprintf('%s = %.3f',paramstr{aa},...
            eval(['data.' paramstr{aa}]));
        i = i + 1;
    end
end
annotstr{:,end+1} = ['RC_{gSYN} = ' mat2str(gSYNs)];
annotstr{:,end+1} = ['RC_{netcon} = ' mat2str(RCs)];

end