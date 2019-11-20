classdef PoissonExponential < bgam.Dlink
    %An implementation of the Poisson/exponential distro/inverse link
    %ie Poisson regression
    properties
    end
    

    methods
        function [ll] = computeLikelihoodAfterNl(~,y,x)
            ll = y'*log(x+1e-12) - sum(x);
        end
        
        function [ll] = computeLikelihood(~,y,eta)
            ll = y'*eta-sum(exp(eta));
        end
        
        function [g,H,L] = computeDlikelihoodX(~,y,eta,X)
            r = exp(eta);
            g = X'*(y-r);
            H = -X'*bsxfun(@times,X,r);
            if nargout > 2
                L = y'*eta-sum(r);
            end
        end
        
        function [g] = computeDlikelihoodEta(~,y,eta)
            g = y - exp(eta);
        end
    end    
end

