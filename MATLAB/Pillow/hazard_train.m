function [Trains] = hazard_train(H,dt,T,numtrials)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [Trains] = hazard_train(H,dt,T,numtrials) returns spike trains 
% by an inhomogenous Poisson process taking a hazard function as its rate.
% H - hazard function
% dt - time step
% T - trial length
% numtrials - number of trials
% Trains - numtrials x T/dt matrix of spike trains
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Trains = zeros(numtrials,T/dt);

for i = 1:numtrials
    spidx = 0;
    while spidx < T/dt
        sp = rand(size(H)) <= dt*H;
        while isempty(sp)
            sp = rand(size(H)) <= dt*H;
        end
        idx = find(sp,1);
        if idx < 1/dt
            spidx = spidx+idx;
            if spidx > T/dt
                break
            end
            Trains(i,spidx) = 1;
        end
    end      
end