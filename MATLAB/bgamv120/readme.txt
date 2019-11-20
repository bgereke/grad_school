bgam - Boosted Generalized Additive Models package
---
Implements boosting for the Generalized Additive and Linear Models (GAM and GLM). 
Extensible, fully documented. Implements linear and stub learners, 
least-squares/logistic/Poisson regression.

The generalized linear model (GLM) is a flexible generalization of ordinary 
least squares regression. The GLM generalizes linear regression by allowing 
the linear model to be related to the response variable via a link function 
and by allowing the magnitude of the variance of each measurement to be a 
function of its predicted value. (Wikipedia)

A common example of a GLM is binomial-logistic distribution/inverse link 
GLM (aka logistic regression), where:

eta = X*w, y ~ Binomial( logistic (eta ))

This GLM allows one to tackle classification problems (where the output is 0 or 1)
in a quasi-linear way.

The generalized additive model (GAM) is a generalization of the GLM where the internal
dynamics are nonlinear, but nevertheless additive:

eta_i = f_1(X^(i,1)) + f_2(X^(i,2)) + ...

f_i are known as smoothers or (in the context of boosting) as learners. Boosting is 
a method of fitting GAMs and by extension GLMs by building up a model (eta) iteratively, 
by, at every iteration, adding to the model the learner most similar to the gradient 
of the likelihood with respect to eta. Regularization is usually done by early-stopping
where the optimal number of iterations is determined through validation.

bgam is a well-documented package that implements boosting with GAMs. 
It currently implements linear learners and stubs (depth-1 trees). Implemented distro-link 
combos include Gaussian/identity, Binomial/Logistic, Poisson/exponential. The package is 
object-oriented and new distro-link combos and learners can be implemented and used 
with ease. The package includes facilities for cross-validation, including a parallel implementation
supporting the parallel computing toolbox. It also allows a subset
of the data to be used at any boosting iteration (stochastic gradient boosting).

Open up TestBgam.m in the editor for several usage examples.

Contributions and requests for new features are welcome.
Author: Patrick Mineault (patrick DOT mineault AT gmail DOT com)

History:

02/07/2011 - 1.2.0 - Parallel cross-validation added
04/03/2011 - 1.1.0 - Added .mex file for StubLearner
                     Tweaked cross-validation code
24/01/2011 - 1.0.1 - Removed temporary .m~ files from .zip
23/01/2011 - 1.0.0 - Initial release

References:

Friedman, Hastie and Tibshirani. Additive logistic regression: a 
statistical view of boosting. Ann. Statist. Volume 28, Number 2 (2000), 
337-407.
Bühlmann and Hothorn. Boosting Algorithms: Regularization, Prediction
and Model Fitting. Statist. Sci. Volume 22, Number 4 (2007), 477-505.
Wood. Generalized Additive Models: an introduction with R. CRC Press,
2006.
Hastie, T. J. and Tibshirani, R. J. (1990). Generalized Additive
Models. Chapman & Hall/CRC.

