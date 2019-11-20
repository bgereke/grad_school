classdef Stub < bgam.Trainer
    %bgam.train.Stub - Stub trainer
    %A stub is a function f(x) such that f(x) = 1 when x <= thresh
    %                                         = 0 otherwise
    %Stubs can be trained efficiently because the inner product of f(x) and y
    %is the sum of all y that are smaller than thresh, which can be
    %computed for all thresh by shuffling y and taking a cumulative sum..
    properties(Hidden)
        fudgeFactor = 1e-12;
    end
    
    methods(Hidden=true)
        function initData(this,X)
            [sorted,ranks] = sort(X);
            mask = [diff(sorted) ~= 0;false(1,size(X,2))];
            ips = 1./sqrt((1:size(X,1))'.*(size(X,1)-1:-1:0)'+this.fudgeFactor);
            %save column ranks, mask because of potential repeated values
            cache.ips = ips;
            cache.ranks = uint32(ranks); %1/2 the memory
            cache.mask = mask;
            this.cache = cache;
        end
        
        function [lrn] = findBestLearner(this,r,X)
            %Reshuffle y according to ranks of columns of X, then take
            %cumulative sum to obtain inner product between step function
            %and y
            if exist('getBestSpotStubLearner') == 3
                %Call a C function to do this
                [bestidxs,bestccs] = getBestSpotStubLearner(r,this.cache.ranks,this.cache.ips,this.cache.mask);
                [~,bestidx] = max(bestccs);
                bestaidx = bestidxs(bestidx) + 1;
            else
                %Do it in m code (much slower)
                cumsums = cumsum(r(this.cache.ranks),1);
                %Convert to correlation
                scaledccs = abs(bsxfun(@times,cumsums,this.cache.ips).*this.cache.mask);
                [~,bestspot] = max(scaledccs(:));
            
                %Find best column
                [bestaidx,bestidx] = ind2sub(size(X),bestspot);
            end
            %Find best threshold
            sortedc = sort(X(:,bestidx));
            
            lrn.gain = 1;
            lrn.index = bestidx;
            lrn.thresh = sortedc(bestaidx);
            lrn.offset = -bestaidx/size(X,1);
        end
    end
    
    methods
        function [y] = evaluate(~,lrn,X)
            y = lrn.gain*(lrn.offset + (X(:,lrn.index)<=lrn.thresh));
        end
        
        %Evaluate the function f(x) for a given covariate index
        function [y] = evaluateOverRange(this,learners,index,x)
            idx = find(cell2mat({learners(:).index})==index);
            y = zeros(length(x),1);
            x = x(:);
            for ii = idx
                y = y + learners(ii).gain*(learners(ii).offset + (x<=learners(ii).thresh));
            end
        end
    end
end