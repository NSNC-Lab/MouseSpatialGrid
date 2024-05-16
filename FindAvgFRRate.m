for h = 1:100
    holder(h)=length(find(results(h).R2On_V_spikes == 1))/3.5;
    %holder2(h)=length(find(results(h).R1On_V_spikes == 1))/3.5;
end

mean(holder)
std(holder)
%mean(holder2)

%0.4 fp
%FR S2Onoff 72
%FR S1Onoff 108

%Very low varience between trials

%0 fp

%Might be worth rerunning offset only at lower e->e