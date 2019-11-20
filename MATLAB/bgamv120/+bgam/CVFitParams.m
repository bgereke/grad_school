classdef CVFitParams < bgam.FitParams
    %bgam.CVFitParams - A dummy class to enter parameters related to a bgam fit    
    %                   Inherits from bgam.FitParams
    %
    %bgam.CVFitParams Properties:
    %    blockSize - Number of iterations to do on a fold before display
    %    iterMargin - How many iterations without decrease in CV deviance
    %                 before giving up
    %    devianceMargin - Number of units of deviance relative to minimum
    %                     that are also considered minima
    %    parallel - if true, fit all folds in parallel using parfor
    %
    %See Also bgam.FitParams
    properties
        %How many boosting iterations to perform on a fold before displaying
        blockSize = 50;
        
        %Determines how many boosting iterations without a decrease in the
        %CV deviance before deciding to go forward with the final fit
        iterMargin = 50;        
        
        %If greater than 0, will select a final number of iterations not
        %corresponding to the actual minimum of CV deviance but to the
        %minimum number of iterations corresponding to CV deviance with
        %devianceMargin of the minimum
        devianceMargin = 0;
        
        %If true, fit all folds in parallel with parfor. Requires parallel
        %computing toolbox, matlabpool open must have been called
        %beforehand
        parallel = 0;
        

        displayStartOfBlockFun = @bgam.CVFitParams.defaultDisplayStartOfBlock
        displayEndOfBlockFun   = @bgam.CVFitParams.defaultDisplayEndOfBlock
        displayBlockFun        = @bgam.CVFitParams.defaultDisplayBlock
        displayHeaderFun       = @bgam.CVFitParams.defaultDisplayHeader
        displayEndOfCV         = @bgam.CVFitParams.defaultDisplayEndOfCV
    end
    
    methods(Static)
        function [] = defaultDisplayStartOfBlock(ii)
            fprintf('  %8.4g      ',ii);
        end
        function [] = defaultDisplayEndOfBlock(d)
            fprintf('    %8.2f\n',d);
        end
        function [] = defaultDisplayBlock(d)
            fprintf('%8.1f ',d)
        end
        function [] = defaultDisplayEndOfCV(d)
            fprintf('\nEnd of cross-validation.\nOptimal number of iterations: %d\nStarting main fit.\n\n',d)
        end
        
        function [] = defaultDisplayHeader(kfold)
            fprintf('\nStarting cross validation\n\n');
            fprintf('                                      CV deviance\n');
            fprintf('    #iters        ');

            for ii = 1:kfold
                fprintf('fold %d   ', ii);
            end

            fprintf('    total\n');
        end
    end
end