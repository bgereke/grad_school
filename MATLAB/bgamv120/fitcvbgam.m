function [thefit] = fitcvbgam(y,X,trainer,fitParams,folds)
    %[thefit] = fitcvbgam(y,X,trainer,fitParams,folds)
    %Computes a boosted GAM with number of boosting of iterations
    %determined by cross validation
    %
    %y,X,trainer are as in fitbgam
    %fitParams is an instance of bgam.CVFitParams which is like
    %bgam.FitParams but with a couple of extra parameters
    %folds is the same height as X and is composed of indices from 1 to N
    %for N-fold cross-validation. Determines which observation is included
    %in which fold. 
    %
    %Cross-validation goes roughly like this:
    %    for each fold -> ii
    %         Call fitbgam with X(~(folds==ii),:)
    %         Compute CV likelihood of models on X(fold==ii,:)
    %    Determine iteration with maximum likelihood
    %    Boost whole model for that number of iterations
    %
    %In reality the fits are performed in small batches of say 50
    %iterations, all folds are cycled through, 50 more iterations, etc.
    %until the likelihood has been stable for a bit, both to provide
    %quicker feedback and stop cross-validating when it seems like the
    %optimal number of iterations is long past.
    %
    %thefit is an instance of bgam.CVFit which is a subclass of bgam.Fit
    %with a couple of extra parameters
    %
    %See also TestBgam, fitbgam, bgam.CVFit, bgam.CVFitParams
    nfolds = max(folds);
    for ii = 1:nfolds
        thefits{ii} = bgam.Fit;
    end
    
    minDeviances = zeros(nfolds,1);
    maxDeviances = zeros(nfolds,1);
    refFitParams = fitParams;
    
    oldDisplayFreq = fitParams.displayFreq;
    refFitParams.displayFreq = Inf;
    
    if length(refFitParams.b) == 1
        refFitParams.b = zeros(length(y),1);
    end
    
    refFitParams.displayHeaderFun(nfolds);

    jj = 0;
    etas = cell(nfolds,1);
    cvdeviances = Inf*ones(refFitParams.niters,nfolds);
    
    if fitParams.parallel && matlabpool('size') == 0
        warning('bgam:parallel','Parallel was selected but matlabpool size is 0');
    end
    
    while jj < refFitParams.niters
        
        if jj == 0
            fitParams.displayStartOfBlockFun(0);
        else
            fitParams.displayStartOfBlockFun(jj+fitParams.blockSize);
        end
        
        rg = [1,fitParams.blockSize]+jj;
        rg = rg(1):rg(2);

        if jj == 0
            for ii = 1:nfolds
                etas{ii} = zeros(nnz(folds==ii),1);
            end
        end
        
        if fitParams.parallel

            parfor ii = 1:nfolds
                [thefits{ii},cvdeviancep,etas{ii},minDeviance,maxDeviance] = fitAFold(ii,jj,folds,refFitParams,thefits{ii},y,X,etas{ii},minDeviances(ii),maxDeviances(ii),trainer);
                minDeviances(ii) = minDeviance;
                maxDeviances(ii) = maxDeviance;
                cvdeviances(rg,ii) = cvdeviancep;
            end

            for ii = 1:nfolds
                if jj > 0
                    tgt = cvdeviances(jj+fitParams.blockSize,ii);
                else
                    tgt = maxDeviances(ii)-minDeviances(ii);
                end

                fitParams.displayBlockFun(tgt);
            end
        else
            for ii = 1:nfolds
                [thefits{ii},cvdeviancep,etas{ii},minDeviance,maxDeviance] = fitAFold(ii,jj,folds,refFitParams,thefits{ii},y,X,etas{ii},minDeviances(ii),maxDeviances(ii),trainer);
                minDeviances(ii) = minDeviance;
                maxDeviances(ii) = maxDeviance;
                cvdeviances(rg,ii) = cvdeviancep;

                if jj > 0
                    tgt = cvdeviances(jj+fitParams.blockSize,ii);
                else
                    tgt = maxDeviances(ii)-minDeviances(ii);
                end

                fitParams.displayBlockFun(tgt);
            end
        end
        
        if jj == 0
            fitParams.displayEndOfBlockFun(sum(maxDeviances)-sum(minDeviances));
            fitParams.displayStartOfBlockFun(jj+fitParams.blockSize);
            for ii = 1:nfolds
                fitParams.displayBlockFun(cvdeviances(fitParams.blockSize+jj,ii));
            end
        end
        fitParams.displayEndOfBlockFun(sum(cvdeviances(fitParams.blockSize+jj,:)));
        
        jj = jj + fitParams.blockSize;
        
        
        sumds = sum(cvdeviances,2);
        [themin minidx] = min(sumds);
        
        if jj - minidx > fitParams.iterMargin
            %Done
            break;
        end
    end
    
    
    cvdeviances = cvdeviances(1:jj,:);
    
    %Determine number of iterations to boost
    if minidx == refFitParams.niters
        niters = refFitParams.niters;
    else
        niters = find(sumds-themin <= refFitParams.devianceMargin,1);
    end
    
    fitParams.displayEndOfCV(niters);
    
    %Now fit the whole things
    refFitParams.niters = niters;
    refFitParams.displayFreq = oldDisplayFreq;
    thefit2 = fitbgam(y,X,trainer,refFitParams);
    
    %Can't figure out how to typecast a bgam.Fit into a bgam.CVFit, hence
    %this hack
    thefit = bgam.CVFit;
    fn = properties(thefit2);
    for ii = 1:length(fn)
        thefit.(fn{ii}) = thefit2.(fn{ii});
    end
    
    thefit.cvdeviances = cvdeviances;
    
    thefit.cvd2s = 1-sum(cvdeviances,1)'/sum(maxDeviances);
    
    thefit.fitParams = refFitParams;
    thefit.cvdeviance = sum(cvdeviances(niters,:));
    thefit.subfits = thefits;
end

function [thefit,cvdeviance,eta,minDeviance,maxDeviance] = fitAFold(ii,jj,folds,refFitParams,thefit,y,X,eta,minDeviance,maxDeviance,trainer)
    valset = folds==ii;
    fitset = ~valset;

    fitParams = refFitParams;
    fitParams.niters = fitParams.blockSize + jj;
    fitParams.b = fitParams.b(fitset);

    if jj ~= 0
        fitParams.restartInfo = thefit;
    end

    if iscell(X)
        pX = {X{1}(fitset,:),X{2:end}};
        thefit = fitbgam(y(fitset),pX,trainer,fitParams);
    else
        thefit = fitbgam(y(fitset),X(fitset,:),trainer,fitParams);
    end

    rg = [1,fitParams.blockSize]+jj;
    addedX = thefit.cumEvaluateRange(X(valset,:),rg);
    rg = rg(1):rg(2);
    
    ll = refFitParams.dlink.computeLikelihood(y(valset),...
            bsxfun(@plus,bsxfun(@plus,refFitParams.b(valset)+eta,addedX),thefit.cs(rg)'));
    eta = addedX(:,end) + eta;

    if jj == 0
        %compute max and min deviance on fold
        minDeviance = -2*refFitParams.dlink.computeLikelihoodAfterNl(...
            y(valset),y(valset));
        maxDeviance = -2*refFitParams.dlink.computeLikelihood(...
            y(valset),ones(nnz(valset),1)*thefit.cs(1)+refFitParams.b(valset));
    end

    cvdeviance = -2*ll - minDeviance;
end