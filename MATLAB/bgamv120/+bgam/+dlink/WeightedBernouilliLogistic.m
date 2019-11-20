classdef WeightedBernouilliLogistic < bgam.Dlink
    %Just like BinomialLogistic for binary outputs only, with different
    %weights for 0 and 1 classes (special use). 
    properties
        weights = [1,1];
    end

    methods
        
        function [ll] = computeLikelihoodAfterNl(this,y,x)
            ll =  this.weights(1)*sum(log(1-x(y==0)+eps)) + this.weights(2)*sum(log(x(y==1)+1e-12));
        end
        function [ll] = computeLikelihood(this,y,eta)
            ll = -this.weights(1)*sum(safelogexp(-eta(y==0,:))) - this.weights(2)*sum(safelogexp(eta(y==1,:)));
        end
        
        function [g,H,L] = computeDlikelihoodX(this,y,eta,X)
            %eta = single(eta);
            leta = logistic(eta);
            rb = zeros(size(y));
            rb(y==0) = -this.weights(1)*leta(y==0);
            rb(y==1) =  this.weights(2)*(1-leta(y==1));
            g = X'*rb;
            H = -X'*bsxfun(@times,X,leta.*(1-leta).*(y*this.weights(2)+(1-y)*this.weights(1)));
            if nargout > 2
                L = this.computeLikelihood(y,eta);
            end
        end
        
        %Compute first derivative of the likelihood with
        %respect to eta
        function [g] = computeDlikelihoodEta(this,y,eta)
            p = logistic(eta);
            g = zeros(size(y));
            g(y==1) =  this.weights(2)*(1-p(y==1));
            g(y==0) = -this.weights(1)*p(y==0);
            %g = this.weights(2)*(1-p(y==1))-this.weights(1)*p(y==0);
        end
    end
end

function [y] = safelogexp(y)
    y = -y;
    %Avoids overflows
    idx1 = y < -20;
    idx2 = y > 20;
    y(idx1) = exp(y(idx1));
    y(~idx1 & ~idx2) = log(1+exp(y(~idx1 & ~idx2)));
end

function [y] = logistic(y)
    y = 1./(1+exp(-y));
end
