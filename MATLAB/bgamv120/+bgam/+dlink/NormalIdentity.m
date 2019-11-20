classdef NormalIdentity < bgam.Dlink
    %An implementation of the identity link/normal distribution
    %ie least-squares
    methods
        function [ll] = computeLikelihoodAfterNl(~,y,x)
            ll = -1/2*sum((y-x).^2);
        end
        
        function [ll] = computeLikelihood(~,y,eta)
            ll = -1/2*sum(bsxfun(@minus,y,eta).^2);
        end
        
        function [g,H,L] = computeDlikelihoodX(~,y,eta,X)
            g = X'*(y-eta);
            H = -X'*X;
            L = -1/2*sum((y-eta).^2);
        end
        
        function [g] = computeDlikelihoodEta(~,y,eta)
            g = y - eta;
        end
    end
end

