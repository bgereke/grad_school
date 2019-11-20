classdef Linear < bgam.Trainer
    %bgam.train.Linear - Linear trainer
    %
    %See also bgam.Trainer
    %
    %bgam.train.Linear Methods: 
    %    learnersToPath - Translates list of learners to a path of w's
    %    learnersToW    - Translates list of learners to a single w
    properties(Hidden)
        %Numerically stabilizes std computation
        fudgeFactor = 1e-12;
    end
    
    properties(Access=private)
        %Otherwise doc doesn't parse the class correctly
        nothingness;
    end
    
    methods(Hidden)
        %Computes std of columns
        function initData(this,X)
            this.cache = sqrt(sum(bsxfun(@minus,X,mean(X)).^2)+this.fudgeFactor);
        end
        
        %Finds column with highest correlation to current residual
        function [lrn] = findBestLearner(this,r,X)
            scaledcc = (X'*r)./this.cache';
            [~,bestidx] = max(abs(scaledcc));
            
            lrn.gain = 1;
            lrn.index = max(bestidx);
        end
    end
    
    methods(Access=public)
        %Trivial
        function [y] = evaluate(~,lrn,X)
            y = lrn.gain*X(:,lrn.index);
        end
        
        %this.learnersToPath(~,learners,numws)
        function [wpath] = learnersToPath(~,learners,nw)
            idx = sub2ind([nw,length(learners)],cell2mat({learners(:).index}),...
                    1:length(learners));
            wpath = zeros(nw,length(learners));    
            wpath(idx) = cell2mat({learners(:).gain});
            wpath = cumsum(wpath,2);
        end
        
        %this.learnersToW(learners,numws)
        function [wpath] = learnersToW(~,learners,nw)
            idx = sub2ind([nw,length(learners)],cell2mat({learners(:).index}),...
                    1:length(learners));
            wpath = zeros(nw,length(learners));    
            wpath(idx) = cell2mat({learners(:).gain});
            wpath = sum(wpath,2);
        end
    end

end