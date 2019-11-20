classdef LinearParams < bgam.LearnerParams
    %LinearLearnerParams - stores the parameters for a linear learner
    properties
        factory = bgam.lrn.LinearFactory;
        learner = bgam.lrn.Linear;
        %Used in preventing divide by zero errors in computing standard deviations
    end
end