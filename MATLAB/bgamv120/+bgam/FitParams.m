classdef FitParams
    %bgam.FitParams - A dummy class to enter parameters related to a bgam fit    
    %
    %bgam.FitParams Properties:
    %    dlink - Distribution/Link, of type bgam.Dlink
    %    niters - Number of boosting iterations
    %    beta - Greediness of boosting procedure
    %    b - start offset
    %    optimizeStepSize - if true, optimize stepsize every iteration
    %    fitFraction - Percentage of observations used every boosting
    %                  iteration
    %    displayFreq - Number of boosting iterations to wait until display
    %
    %See Also bgam.Dlink
    properties
         %Specify distribution/link combo
         %See Also bgam.Dlink
         dlink  = bgam.dlink.NormalIdentity; 
         
         % Number of boosting iterations
         niters = 500; 
            
         %Greediness of boosting procedure
         beta   = 0.1; 
            
         %A vector the size of y specifying a start offset at each
         %observation
         b      = 0;
         
         %Whether to explicitly optimize the step size on each boosting 
         %iteration or rely on an estimate (ie as many Netown iterations as 
         %required or just one). false is faster, slightly less accurate
         optimizeStepSize = true;
            
         %The percentage of y values used at every point of the boosting
         %procedure. Smaller than 1 changes the boosting
         %procedure to stochastic gradient boosting.
         fitFraction = 1;
         
         %The number of boosting iterations to wait until displaying the current
         %likelihood
         displayFreq = 50;
            
         %A variable for internal use specifying an old fit to continue
         %boosting
         restartInfo = []; 
         
    end
    
    methods      
        function [] = displayFun(this,iter,d,maxd)
            if mod(iter,this.displayFreq) == 0
                fprintf('% 6d      %8.1f      %.3f\n', iter,d,1-d/maxd);
            end
        end
        
        function [] = displayChFun(this)
            if this.displayFreq < Inf
                fprintf('   Iter      Deviance      D^2\n');
            end
        end
    end
end