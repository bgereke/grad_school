classdef BinomialLogistic < bgam.Dlink
    %An implementation of the Binomial/logistic distro/inverse link
    %ie Logistic regression
    %Set ntrials to more than 1 for repeated observations
    properties
        ntrials = 1;
    end

    methods
        
        function [ll] = computeLikelihoodAfterNl(this,y,x)
            ll =  y'*log(x/this.ntrials+1e-12) + (this.ntrials - y)'*log(1-x/this.ntrials+eps);
        end
        function [ll] = computeLikelihood(this,y,eta)
            ll = -(this.ntrials-y)'*eta-this.ntrials*sum(safelogexp(eta));
        end
        
        function [g,H,L] = computeDlikelihoodX(this,y,eta,X)
            leta = logistic(eta);
            g = X'*(y-this.ntrials*leta);
            H = -this.ntrials*(X'*bsxfun(@times,X,leta.*(1-leta)));
            if nargout > 2
                L = -(this.ntrials-y)'*eta+this.ntrials*sum(log(leta));
            end
        end
        
        %Compute first derivative of the likelihood with
        %respect to eta
        function [g] = computeDlikelihoodEta(this,y,eta)
            g = y - this.ntrials*logistic(eta);
        end
    end
end

function [y] = safelogexp(y)
    y = -y;
    idx1 = y < -6;
    idx2 = y > 5;
    y(idx1) = exp(y(idx1));
    y(~idx1 & ~idx2) = log(1+exp(y(~idx1 & ~idx2)));
end

function [y] = logistic(y)
    y = 1./(1+exp(-y));
end
