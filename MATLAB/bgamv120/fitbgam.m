function [thefit] = fitbgam(y,X,trainer,fitParams)
    %function [thefit] = fitbgam(y,X,trainer,fitParams)
    %Fit boosted generalized additive model
    %A GAM is composed of an additive component:
    %
    % eta_k = sum_i learner_i(X_kj)
    %
    %The learner can be an arbitrarily nonlinear function, although
    %generally it is taken to be low-dimensional (acting on only a few of 
    %the columns of X simultaneously) to reduce the number of
    %parameters to fit. If the learner is linear then we obtain a
    %generalized linear model.
    %
    %The measured response y_k is assumed to be related to the internal
    %response by:
    %
    % E(y_k) = h(eta_k + eta0_k) = mu_k
    %
    %h^-1 is known as the link, and h as the inverse link. eta0_k is an offset.
    %Furthermore:
    %
    % y_k ~ SomeExponentialFamilyDistribution(mu_k)
    %
    %For example, if Distro = Normal, h = Identity -> Linear regression
    % if Distro = Binomial, h = Logistic function  -> Logistic regression
    % if Distro = Poisson,  h = Exponential        -> Poisson regression
    %
    %A model can be fit by maximizing the likelihood L(y,eta+eta0) of the data, and
    %this is a convex problem provided some minimal conditions on h.
    %Boosting approaches the problem a little differently by doing this in
    %a stepwise fashion. Suppose that the current eta has a given value.
    %Then d(L(y,eta+eta0))/d((eta)_k) is a derivative of the likelihood, so by
    %moving eta in the direction of this derivative we can increase the
    %likelihood. Boosting consists of finding, at every step,
    %the parameters of the learner that best matches the derivative 
    %d(L(y,eta+eta0))/d((eta)_k) and moving the model (eta) in this direction.
    %For any given learner there is always an optimal stepsize alpha such 
    %L(y,eta+eta0+alpha*learner) is maximized; the boosting procedure generally
    %works better by taking smaller steps, of the form 
    %eta + alpha*beta*learner, where beta < 1.
    %
    %After a great number of iterations the model starts to fit to noise;
    %regularization is done by truncating the model after a certain number
    %of iterations, and this is usually determined by (cross-)validation
    %
    %This implementation of GAM boosting has eta0_k = b_k + c, where b_k is
    %user-given (optional, defaults to 0) and c is an offset refit after
    %every boosting iteration.
    %
    %Parameters:
    %
    %y: observed values, n x 1
    %X: design matrix  , n x m (n observations, m covariates)
    %trainer: An instance of a subclass of bgam.Trainer containing 
    % information about the type of learners one wants to use. For example,
    % setting this to an instance of bgam.train.Linear will yield a model
    % with linear learners, equivalent to a sparse GLM. 
    %fitParams: An instance of bgam.FitParams containing parameters related to
    % the fit, for example, the number of boosting 
    % iterations, the distribution and link, b, beta, etc. 
    %
    %Returns:
    %fit: an instance of bgam.Fit containing the learned learners, the
    %deviance (goodness-of-fit) values, etc.
    %
    %References:
    %Friedman, Hastie and Tibshirani. Additive logistic regression: a 
    %statistical view of boosting. Ann. Statist. Volume 28, Number 2 (2000), 
    %337-407.
    %Bühlmann and Hothorn. Boosting Algorithms: Regularization, Prediction
    %and Model Fitting. Statist. Sci. Volume 22, Number 4 (2007), 477-505.
    %Wood. Generalized Additive Models: an introduction with R. CRC Press,
    %2006.
    %Hastie, T. J. and Tibshirani, R. J. (1990). Generalized Additive
    %Models. Chapman & Hall/CRC.
    %
    %See also TestBgam, bgam.Trainer, bgam.FitParams, bgam.Fit
    mindeviance = -2*fitParams.dlink.computeLikelihoodAfterNl(y,y);
    
    if isempty(fitParams.restartInfo)
        ii = 1;
        thefit = bgam.Fit();
        b = fitParams.b;
        if length(b) == 1
            b = zeros(size(y));
        end
        eta = b;
        thefit.cs = zeros(fitParams.niters,1);
        thefit.deviances = zeros(fitParams.niters,1);
        thefit.d2s = zeros(fitParams.niters,1);
        
        %Although it would be natural to have an array of objects here,
        %assignment becomes the bottlenecks. Hence learners is composed of
        %vanilla objects
        thefit.learners = [];
        
        [c,ll] = optimizeC(y,eta,fitParams.dlink,0);
        
        
        thefit.maxdeviance = -2*ll-mindeviance;
        
        fitParams.displayChFun();
        fitParams.displayFun(0,thefit.maxdeviance,thefit.maxdeviance);
    else
        thefit = fitParams.restartInfo;
        oldlen = length(thefit.learners);
        newlen = fitParams.niters;
        dlen = newlen - oldlen;
        
        ii = oldlen + 1;
        
        thefit.deviances = [thefit.deviances;zeros(dlen,1)];
        thefit.cs = [thefit.cs;zeros(dlen,1)];

        c = thefit.cs(oldlen);
        
        b = fitParams.b;
        eta = thefit.eta + b;
        fitParams.restartInfo = []; %To prevent from a fit containing its restart info for a fit that contains its restart info that ....
    end
    
    if fitParams.fitFraction == 1
        trainer.initData(X);
    end
    
    while ii <= fitParams.niters
        %Get current derivative of likelihood
        g = fitParams.dlink.computeDlikelihoodEta(y,eta+c);
        g = g - mean(g);
        
        if fitParams.fitFraction < 1
            rg = rand(size(g))<fitParams.fitFraction;
            g = g(rg);
            if iscell(X)
                pX = {X{1}(rg,:),X{2:end}};
            else
                pX = X(rg,:);
            end
                
            trainer.initData(pX);
            learner = trainer.findBestLearner(g-mean(g),pX);
        else
            %Get learner
            learner = trainer.findBestLearner(g,X);
        end
        
        %Eval learner
        deta = trainer.evaluate(learner,X);
                
        %Move model in direction of learner
        if fitParams.optimizeStepSize
            %Full optimization with backtracking
            [alpha ~] = optimizeA(y,eta+c,deta,fitParams.dlink);
        else
            %Partial optimization
            [g,H] = fitParams.dlink.computeDlikelihoodX(y,eta+c,deta-mean(deta));
            alpha = -g/H;
        end
        
        eta = eta + fitParams.beta*alpha*deta;
        learner.gain = double(alpha*fitParams.beta);
        
        if fitParams.optimizeStepSize
            [c,ll] = optimizeC(y,eta,fitParams.dlink,c);
        else
            c = c - fitParams.beta*alpha*mean(deta);
            [g,H] = fitParams.dlink.computeDlikelihoodX(y,eta+c,ones(size(eta)));
            c = c - g/H;
            ll = fitParams.dlink.computeLikelihood(y,eta+c);
        end
        
        thefit.deviances(ii) = double(-2*ll - mindeviance);
        thefit.d2s(ii) = double(1-(-2*ll - mindeviance)/thefit.maxdeviance);
        thefit.cs(ii) = double(c);
        
        if length(thefit.learners) < fitParams.niters
            if isempty(thefit.learners)
                thefit.learners = learner;
            end
            thefit.learners(fitParams.niters) = learner; %Preassign
            thefit.learners = reshape(thefit.learners,fitParams.niters,1);
        end
        
        thefit.learners(ii) = learner;
        fitParams.displayFun(ii,thefit.deviances(ii),thefit.maxdeviance)
        
        ii = ii + 1;
    end
    trainer.cleanupData();
    
    thefit.fitParams = fitParams;
    thefit.eta = eta - b;
    thefit.trainer = trainer;
    thefit.deviance = thefit.deviances(end);
end

function [c,L] = optimizeC(y,eta,dlink,c0)
    if dlink.isContinuous
        epsilon = 1e-3;
        c = c0;
        [g,H,L] = dlink.computeDlikelihoodX(y,eta+c,ones(size(y),class(y)));
        L0 = L;
        gain = 1;
        
        while 1 %Newton steps with backtracking if necessary
            dc = g/H;
            c = c - gain*dc;
            oldH = H;
            oldg = g;
            oldL = L;
            
            [g,H,L] = dlink.computeDlikelihoodX(y,eta+c,ones(size(y),class(y)));
            
            %Use finite difference instead
            %delta = 1e-8;
            %ga = (dlink.computeLikelihood(y,eta+c+delta)-dlink.computeLikelihood(y,eta+c))/delta
            
            %g
            %ga
            %ga-g
            
            if L < L0
                %backtrack
                c = c + gain*dc;
                gain = gain*.1;
                H = oldH;
                g = oldg;
                L = oldL;
                if gain < 1e-12
                    error('Could not find a direction in which to change alpha and c for the likelihood to increase');
                end
            elseif abs(L - L0) < epsilon
                break;
            else
                L0 = L;
            end
        end
    else
        %TODO
    end
end

function [alpha,c] = optimizeA(y,eta,deta,dlink)
    if dlink.isContinuous
        epsilon = 1e-3;
        ab = zeros(2,1,class(deta));
        X = [deta,ones(size(deta),class(deta))];
        [g,H,L] = dlink.computeDlikelihoodX(y,eta,X);
        L0 = L;
        gain = 1;
            
        while 1 %Newton steps with backtracking if necessary
            dab = H\g;
            ab = ab - gain*dab;
            oldH = H;
            oldg = g;
            oldL = L;
            [g,H,L] = dlink.computeDlikelihoodX(y,eta+X*ab,X);

            if L < L0
                %backtrack
                ab = ab + gain*dab;
                gain = gain*.1;
                H = oldH;
                g = oldg;
                L = oldL;
                if gain < 1e-12
                    warning('bgam:regression','Could not find a direction in which to change alpha and c for the likelihood to increase');
                    break;
                end
            elseif abs(L - L0) < epsilon
                break;
            else
                L0 = L;
            end
        end
        
        alpha = ab(1);
        c = ab(2);
    else
        %TODO
    end
end

