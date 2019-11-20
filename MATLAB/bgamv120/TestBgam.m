%Several examples on how to use bgam.*. Open in Matlab editor and read
%comments

%Create some simulated data for Poisson regression
X = randn(1000,500);
w = exp(-(1:500)'*.2).*randn(500,1);
y = poissrnd(exp(X*w));

%Now infer w given y, X through bgam
%Specify the trainer for linear learners
trainer = bgam.train.Linear(); %Doesn't have adjustable parameters

%Start specifying parameters of the fit
fparams = bgam.FitParams();

%Specifiy the use of Poisson distro and exponential inverse link
fparams.dlink = bgam.dlink.PoissonExponential();
%Type help bgam.Dlink for other distro-link combos

%total number of boosting iterations and greediness:
fparams.niters = 5000;
fparams.beta  = .02;

%And fit this model
tic;
thefit = fitbgam(y,X,trainer,fparams);
toc;

%The linear trainer includes a convenience function to translate the
%sequence of learners sum i to n learners_i(X) -> a path of w in X*w
wpath = trainer.learnersToPath(thefit.learners,500);

%Display wpath:
subplot(2,1,1);semilogx(wpath');
xlabel('# boosting iterations');
ylabel('weights');
%Note only a few ws ever become non-zero. Thus the model is sparse

% Now to decide how many iterations correspond to the best model, we might 
% use a secondary dataset on which the data is untrained and evaluate
% the quality of the model after each iteration
% Create this second set of simulated data
X2 = randn(1000,500);
y2 = poissrnd(exp(X2*w));
etapath = thefit.cumEvaluate(X2);
qualities = fparams.dlink.computeLikelihood(y2,bsxfun(@plus,etapath,thefit.cs'));

%And plot to view the validated log likelihood as a function of number
%of boosting iterations
subplot(2,1,2);cla;semilogx(qualities);
hold on;
[thebest,bestiter] = max(qualities);
plot(bestiter,qualities(bestiter)+150,'vr','MarkerSize',10);
hold off;
xlabel('#iterations');
ylabel('Validated goodness of fit');

%Usually this should be a little less than 5000, but boosting is
%relatively insensitive to total number of iterations (as long as it's the
%right order of magnitude)

%%

%Example 2: Stub learner, binomial/logistic
%Here there are two covariates, and the relationship between the covariates
%and the response is very nonlinear.
%eta = sin(pi*x_1) + sin(pi*2*x_2)

clear thefit;clear trainer;
%Generate some data
X = 2.5*(rand(500,2)-.5);
eta = 3*sum(sin(pi*X*diag([1,2])),2);
y = binornd(1,1./(1+exp(-eta)));

%Now infer w given y, X through bgam
trainer = bgam.train.Stub(); %Doesn't have adjustable parameters

%Start specifying parameters of the fit
fparams = bgam.FitParams();

%Specifiy the use of Binomial distro and logistic inverse link
fparams.dlink = bgam.dlink.BinomialLogistic();
%Type help bgam.Dlink for other distro-link combos

%total number of boosting iterations and greediness:
fparams.niters = 5000;
fparams.beta  = .02;

%fitfraction = 1 is faster (much faster if you compile the mex file, see
%getBestSpotStubLearner.readme.txt), but generally the results are less
%visually appealing than fitFraction << 1 (although the diff. in likelihood
%between the two is marginal)
fparams.fitFraction = 1;

%And fit this model
tic;
thefit = fitbgam(y,X,trainer,fparams);
toc;

%Plot the learned functions using evaluateOverRange convenience function
rg = -1.25:.01:1.25;
subplot(1,2,1);plot(rg,trainer.evaluateOverRange(thefit.learners,1,rg),rg,3*sin(pi*rg),X(:,1),-ones(size(X,1),1),'.');
subplot(1,2,2);plot(rg,trainer.evaluateOverRange(thefit.learners,2,rg),rg,3*sin(2*pi*rg),X(:,1),-ones(size(X,1),1),'.');

%Note that the gains are typically off because of a shrinkage effect

%%
%Example 3, cross-validation, Normal/identity (least squares), linear
%learner

%Random data
X = randn(500,500);
w = exp(-(1:500)'*.2).*randn(500,1);
y = X*w + randn(size(X,1),1);

%Note equal number of parameters and observations, plus some noise.
%Regularization (in the way of early stopping) is important in this case.

%Now infer w given y, X through bgam
%Specify the trainer for linear learners
trainer = bgam.train.Linear(); %Doesn't have adjustable parameters

%Start specifying parameters of the fit
fparams = bgam.CVFitParams();

%Specifiy the use of Poisson distro and exponential inverse link
fparams.dlink = bgam.dlink.NormalIdentity();
%Type help bgam.Dlink for other distro-link combos

%total number of boosting iterations and greediness:
fparams.niters = 5000;
fparams.beta  = .01;

%Assign which observation goes in which validation fold (here 10-fold
%validation is used)
folds = ceil(rand(size(y))*10);

%And fit this model with cross-validation
tic;
thefit = fitcvbgam(y,X,trainer,fparams,folds);
toc;

%And plot to view the validated log likelihood as a function of number
%of boosting iterations
semilogx(sum(thefit.cvdeviances'));
xlabel('#iterations');
ylabel('CV deviance (lower = better)');

%%
%Cross-validation can be parallelized through the use of parfor
%Requires parallel computing toolbox
%Do this once to initialize parallel environment:
matlabpool open;

%%
%Specify the use of parallel processing:
fparams.parallel = 1;
fparams.blockSize = 200; %Increase block size to minimize communication overhead

%This requires a large amount of RAM
%Fit
tic;
thefit = fitcvbgam(y,X,trainer,fparams,folds);
toc;

%On my 6-core i7 970 CPU, this is faster by a factor 2.7x