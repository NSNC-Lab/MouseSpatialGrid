%Generate 10 spike trains per song

%Song 1 -- Randomly generated but the same accross all trains for 1





%Song 2


%Add Noise for each pair of songs -- (bit flip per song pair i.e. song1a and song2a)




%1. Random choise classification

%Choose a template pair

%Calculate discrimitability for other pairs

%Compare these values


%2. Matrix classification 


%Choose # of trials

n_trials = 2;


one_mat = ones([n_trials,n_trials]);
zero_mat = zeros([n_trials,n_trials]);


%Calcualte a distMat
%distMat = [zero_mat,zero_mat;one_mat,zero_mat];

%Run a few times just to make sure it converges on 50%

corrects = [];

%for epochs = 1:10000

%distMat = rand([n_trials*2,n_trials*2]);
%Create n by n symetric matrix

A = rand(n_trials*2); 
distMat = (A + A') / 2;
distMat(1:n_trials*2+1:end) = 0; % Set the diagonal elements to 0
%This is a silly trick since you can index a matrix at point k by doing
%row*cloumn length +remainder. n+1 will always be the diagnol.




figure();
imagesc(distMat)

%Do matrix comparison

total = 0;
correct = 0;

for song = 1:2
    
    for m = 1:n_trials
        for k = 1:n_trials
            for j = 1:n_trials
                if j ~= k
                    
                    if song == 1
                        
                        if distMat(k,j) == distMat(m+n_trials,j)
                            correct = correct + 0.5;
                            comparison = 0;
                        else
                            comparison = distMat(k,j) < distMat(m+n_trials,j);
                        end
                    else
                        
                        if distMat(m,j+n_trials) == distMat(k+n_trials,j+n_trials)
                            correct = correct + 0.5;
                            comparison = 0;
                        else
                            comparison = distMat(m,j+n_trials) > distMat(k+n_trials,j+n_trials);
                        end
                    end
    
                    if comparison == 1
                        correct = correct + 1;
                    end
        
                    total = total + 1;
       
                    
                end
            end
        end
    end
end

total;

correct = correct/total;

corrects = [corrects,correct];


%end

mean(corrects)

