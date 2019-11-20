classdef Trainer < handle
    %bgam.Trainer - A class that:
    % - stores parameters relevant to learners of a given type
    % - finds a learner of a given type that best matches a residual
    % - evaluates a learner at a point
    %
    %Trainer is an abstract class only; only subclasses of Trainer should
    %be instantiated. Currently available trainers are:
    % - bgam.train.Linear - trainer for linear learners
    % - bgam.train.Stub   - trainer for stub learners (step functions)
    %
    %A trainer must implement the function evaluate:
    %   function [r] = evaluate(this,learner,X)
    %This evaluates a learner (for example one stored in the learners
    %property of an instance of bgam.Fit) for the design matrix X.
    %
    %A trainer may also implement convenience functions which are only
    %relevant to its learners, for example to visualize a batch of learners
    %of a particular type.
    %
    %The following information is of relevance only if you want to
    %implement your own trainer.
    %---------------------------
    %During learning (in fitbgam for example), initData(this,X) is called
    %once (at the start), while findBestLearner(this,residual,X) is called
    %at every boosting iteration. The purpose of initData(this,X) is to
    %compute and store information that is needed when finding the best 
    %learner; this info should only depends on the data (X). For example, for a
    %linear learner the standard deviation of the columns of X needs to be
    %computed, but it is wasteful to compute this every time
    %findBestLearner is called. The information is only relevant during
    %learning; it should be stored in .cache, and will be discarded after
    %learning (fitbgam calls cleanupData after learning everything).
    %
    %If, however, FitParams.fitFraction is set to less than 1, then the data
    %is a different subset of X on every boosting iteration. In this case
    %initData is called every boosting iteration right before
    %findBestLearner.
    %
    %The function 
    %   [learner] = findBestLearner(this,residual,X)
    %Should find the best learner of the appropriate class given data and a
    %residual. In the linear case for example this is the index of the
    %column of X for which abs(corr(residual, X)) is maximized. learner
    %should be a struct with a property 'gain', which should be set to 1. 
    %fitbgam will overwrite the gain property after determining the optimal
    %gain. Although it would be natural for learner to be an object, there
    %is a huge overhead in assigning object to a large vector as opposed to
    %assigning structs. Presumably a Matlab bug.
    %
    %residual always has zero mean and the implementor should take care
    %that the learning procedure is insensitive to the mean of a column, as
    %the implementation always reoptimizes c after every boosting
    %iteration.
    %
    %See bgam.train.Linear for an example of an implementation of
    %bgam.Trainer
    properties(Hidden=true,Access=protected)
        cache %See notes on initData
    end
        
    
    methods(Hidden=true)
        %See notes on initData
        function [] = cleanupData(this)
            this.cache = [];
        end
    end
    
    methods(Abstract,Hidden=true)
        %Stores information relevant to learning in .cache; see class 
        %comment for details.
        [] = initData(this,X)
        %Find best learner given data
        [learner] = findBestLearner(this,residual,X)
    end
    
    methods(Abstract)
        %Evaluates a learner at a given point. Should not rely on .cache
        [r] = evaluate(this,learner,X)
    end
end

