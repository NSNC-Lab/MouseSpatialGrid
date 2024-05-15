%Number of nodes (x data)
%heirarchy of nodes (y data)
channels = [1,2,3,4];
le = length(s.populations)-1;

Nodesx = [];
Nodesy = [];
for b = 1:length(channels)
    for a = 1:length(s.populations)
    
    %b = 1;
        if strcmp(s.populations(a).name,'On')
            Nodesx = [Nodesx channels(b)];
            Nodesy = [Nodesy 0];
        
        elseif strcmp(s.populations(a).name,'Off')
            Nodesx = [Nodesx channels(b)+0.5];
            Nodesy = [Nodesy 0];
        
        elseif strcmp(s.populations(a).name,'ROn')
            Nodesx = [Nodesx channels(b)];
            Nodesy = [Nodesy 2];
        
        elseif strcmp(s.populations(a).name,'ROff')
            Nodesx = [Nodesx channels(b)+0.5];
            Nodesy = [Nodesy 2];
        
        elseif strcmp(s.populations(a).name,'SOnOff')
            Nodesx = [Nodesx channels(b)+0.25];
            Nodesy = [Nodesy 1];
        
        elseif strcmp(s.populations(a).name,'TD')
            Nodesx = [Nodesx channels(b)+0.75];
            Nodesy = [Nodesy 3];
        
        elseif strcmp(s.populations(a).name,'X')
            Nodesx = [Nodesx channels(b)-0.25];
            Nodesy = [Nodesy 1];
        
        end

    end
end

for t = 1:length(s.populations)
    if strcmp(s.populations(t).name,'C')
        c_pos = t;
        c_num = s.populations(t).size;
        
        for h = 1:c_num 
            Nodesx = [Nodesx (length(channels)+1)*1*h/(1+c_num)];
            Nodesy = [Nodesy 3.5];
        end
        
    end
end




%Should be able to somehow automattically populate this with s.connections
s1 = [];
t = [];

inhibs_s = [];
inhibs_t = [];


for c = 1:length(channels)
    for x = 1:length(s.connections)
        str = extractBefore(s.connections(x).direction,"-");
        stra = extractAfter(s.connections(x).direction,">");
        if   strcmp(str,'On') 
            s1 = [s1 (channels(c)-1)*le+1];
        elseif strcmp(str,'Off')
            s1 = [s1 (channels(c)-1)*le+2];
        
        elseif strcmp(str,'ROn')

            s1 = [s1 (channels(c)-1)*le+3];

        elseif strcmp(str,'ROff')
            s1 = [s1 (channels(c)-1)*le+4];
        
        elseif strcmp(str,'SOnOff')
            s1 = [s1 (channels(c)-1)*le+5];
            inhibs_s = [inhibs_s (channels(c)-1)*le+5];
            flag = 1;
    
        elseif strcmp(str,'TD')
            s1 = [s1 (channels(c)-1)*le+6];
        
        elseif strcmp(str,'X')
            if strcmp(stra,'ROn')
                holder = 1;
            else
                s1 = [s1 (channels(c)-1)*le+7];
            end
    
        elseif strcmp(str,'C')
            s1 = [s1 (channels(c)-1)*le+8];
    
        end
    
    
        %%%%%%%%%%
    
        if   strcmp(stra,'On') 
            t = [t (channels(c)-1)*le+1];
            if flag == 1
                flag = 0;
                inhibs_t = [inhibs_t (channels(c)-1)*le+1];
            end


        elseif strcmp(stra,'Off')
            t = [t (channels(c)-1)*le+2];
            if flag == 1
                flag = 0;
                inhibs_t = [inhibs_t (channels(c)-1)*le+2];
            end
        
        elseif strcmp(stra,'ROn')
           
            if flag == 1
                flag = 0;
                inhibs_t = [inhibs_t (channels(c)-1)*le+3];
            end
            if strcmp(str,'X')
                holder = 1;
            else
                t = [t (channels(c)-1)*le+3];
            end
        
        elseif strcmp(stra,'ROff')
            t = [t (channels(c)-1)*le+4];
            if flag == 1
                flag = 0;
                inhibs_t = [inhibs_t (channels(c)-1)*le+4];
            end
        
        elseif strcmp(stra,'SOnOff')
            t = [t (channels(c)-1)*le+5];
            if flag == 1
                flag = 0;
                inhibs_t = [inhibs_t (channels(c)-1)*le+5];
            end
    
        elseif strcmp(stra,'TD')
            t = [t (channels(c)-1)*le+6];
            if flag == 1
                flag = 0;
                inhibs_t = [inhibs_t (channels(c)-1)*le+6];
            end
        
        elseif strcmp(stra,'X')

            t = [t (channels(c)-1)*le+7];

            
            if flag == 1
                flag = 0;
                inhibs_t = [inhibs_t (channels(c)-1)*le+7];
            end
        
        %We may need to find a way to make this dynamic
        elseif strcmp(stra,'C')
            t = [t length(Nodesx)];
            if flag == 1
                flag = 0;
                inhibs_t = [inhibs_t length(Nodesx)];
            end
    
        end
       


    end
end

%Just for example (BE SURE TO REMOVE THIS!!!!!)
% netcons.XRnetcon = [1,1,0,1;
%                     0,1,1,1;
%                     1,1,1,0;
%                     1,0,1,1];


inhibs_t2 = [];
inhibs_s2 = [];
%Now lets handle cross channel connections
%Probably should handle RC netcon here actually eventually.
for q = 1:length(netcons.XRnetcon)
    for p = 1:length(netcons.XRnetcon)
        if netcons.XRnetcon(p,q) == 1
            t = [t (q-1)*le+3];
            inhibs_t2 = [inhibs_t2 (q-1)*le+3];
            s1 = [s1 (p-1)*le+7];
            inhibs_s2 = [inhibs_s2 (p-1)*le+7];
        end
    end
end



%Be careful! Do not override the s struct
%s1 = [1 1 2 2 5 5 3 6 6 3];
%t = [3 5 4 5 3 4 8 3 4 7];
weights = ones(1,length(t));

figure(20)

G = digraph(s1,t,weights);
%G = digraph(s,t)

%IMPORTANT:::: Look at
%https://www.mathworks.com/help/matlab/ref/matlab.graphics.chart.primitive.graphplot.highlight.html
%for graph colors

p = plot(G,'b','Xdata',Nodesx,'Ydata',Nodesy); hold on

%Find the inhibitory nodes
%X
%SonOff
%TD?

highlight(p, inhibs_s, inhibs_t,'EdgeColor','red')
highlight(p, inhibs_s,'NodeColor','red')

highlight(p, inhibs_s2,'NodeColor','red')
highlight(p, inhibs_s2, inhibs_t2,'EdgeColor',[1 0.6 0.6])
highlight(p, inhibs_s2, inhibs_t2,'LineStyle','--')


%Now we have to go through and "highlight" which means to change the size
%of the line widths that connect neurons.
syns = [];

%Might have to do this for inhibs_2/s2
for u = 1:length(s1)
    start_type = mod(s1(u),7);
    end_type = mod(t(u),7);
    
    if start_type == 0
        start_type = 7;
    end
    if end_type == 0
        end_type = 7;
    end

    tstring = s.populations(start_type).name + "->" + s.populations(end_type).name;
    for j = 1:length(s.connections)
        if strcmp(tstring,s.connections(j).direction)
            target = j;
        end
    end

    gsyn_m = s.connections(target).parameters(2);

    if gsyn_m{1} > 0.2
        end_gsyn = 0.5;
    else
        end_gsyn = (gsyn_m{1}^1.5)*750;
    end

    syns = [syns gsyn_m];
    
    highlight(p, s1(u), t(u), 'LineWidth',end_gsyn,'ArrowSize',end_gsyn*5)

    
end



%What are the things that go into the magnitude of the connection?
%Gsyn
%Netcon
    %PEnetcon
    %RCnetcon

%annotation('arrow', [.56 .58], [.221 .223],'color',[0 0 1]);

%ToDos
%Make all of the connections dynamic    %Still needs some work
      %Should make it so that if a new node  is added to the system or
      %taken away it still works.
      %Should work with different cortical neuron setups.

%Extend to 4 channels                   %Done
%Change colors                          %Mostly done
%Change line/arrow weights              %Mostly done
%Put in arrows                          %Done
%Find a way to handle x inhibition without it being too crowded   %Done

%Dales law implement   %Mostly Done

%Things to touch on in today's meeting
%1. Netcon is not implemented yet for line width multiplier
%2. Labels on lines
%3. Further flexibility

%Todos (After meeting)
%Make TDs red
%Replicate Aim network model