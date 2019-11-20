function [Trains] = gammatrain(shape,rate,dt,T,numtrials)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [Trains] = gammatrain(shape,rate,dt,T,numtrials) returns spike 
% trains with gamma distributed isi's
% shape - shape parameter for the gamma distribution
% rate - spike rate of the process
% dt - time step
% T - trial length
% numtrials - number of spike trains
% Trains - numtrials x T/dt matrix of spike trains
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Trains = zeros(numtrials,T/dt);
ISIs = gamrnd(shape,1/(rate*shape),numtrials,T*rate*3);
Times = cumsum(ISIs,2);
Times(Times>T) = 0;
Times = ceil(Times/dt);

for i=1:numtrials
    Trains(i,Times(i,1:sum(Times(i,:)~=0))) = 1;
end
