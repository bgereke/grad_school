function [Trains] = glm(Stim1,k,h,b,dt,numtrials)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [Trains] = glm(Stim1,k,h,b,numtrials) simulates the response of
% a GLM neuron to a 1D stimulus.
% Stim1 - the 1D stimulus
% k - stimulus filter
% h - spike-history filter
% b - constant
% dt - time step (sets rate)
% numtrials - number of trials
% Trains - numtrials x length(Stim1)-length(k) matrix of spike trains
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Trains = zeros(numtrials,length(Stim1));

for i = 1:numtrials
    
    y_hist = zeros(size(h));
    
    for j = 1:length(Stim1)-length(k)
        
        lambda = exp(dot(k,Stim1(j:j+length(k)-1))+dot(h,y_hist)+b);
        
        if dt*lambda > 1
            lambda = 1/dt;
        end
        
        Trains(i,j+length(k)-1) = binornd(1,dt*lambda);
        
        if Trains(i,j+length(k)-1) == 1
            y_temp = y_hist;
            y_hist = zeros(size(h));
            y_hist(2:end) = y_temp(1:end-1);
            y_hist(1) = 1;
        else
            y_temp = y_hist;
            y_hist = zeros(size(h));
            y_hist(2:end) = y_temp(1:end-1);
        end
    end
end
 