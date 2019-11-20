% demo4_1spacetimefilter.m
%
% Tutorial script illustrating maximum likelihood / maximally informative
% dimensions (MID) estimation for LNP model with a single space-time filter.
%
% This demo borrows heavily from demo1, and differs mostly due to the need
% to plot the space-time filter as a 2D image.

% initialize paths
initpaths;

datasetnum = 3;  % for spatial binary white noise stimulus
trainfrac = .8; % fraction of data to use for training (remainder is "test data")

% Load data divided into training and test sets
[Stim_tr,sps_tr,Stim_tst,sps_tst,RefreshRate,filts_true] = loadSimDataset(datasetnum,trainfrac);

slen_tr = size(Stim_tr,1);   % length of training stimulus / spike train
slen_tst = size(Stim_tst,1); % length of test stimulus / spike train
nx = size(Stim_tr,2); % number of spatial pixels
nsp_tr = sum(sps_tr);   % number of spikes in training set
nsp_tst = sum(sps_tst); % number of spikes in test set

%% == 1. Compute STA and estimate (piecewise constant) nonlinearity using histograms ====

nkt = 12; % number of time bins to use for filter 
% This is important: normally would want to vary this to find optimal filter length

% Compute STA
sta = simpleSTC(Stim_tr,sps_tr,nkt);  % compute STA
sta = sta./norm(sta(:));  % normalize sta to be a unit vector

% Estimate piecewise-constant nonlinearity ("reconstruction" method using histograms)
nhistbins = 15; % # histogram bins to use
[fnlhist,xbinedges] = fitNlin_hist1D(Stim_tr, sps_tr, sta, RefreshRate, nhistbins); % estimate 1D nonlinearity 

%% == 2. iSTAC (information-theoretic spike-triggered average and covariance) estimator [OPTIONAL] ===

% (Note that iSTAC is not using the temporal basis used for next
% two models, which would denoise filter estimate slightly; note that
% nevertheless it's nearly as good as estimator with RBF nonlinearity) 

nFilts = 1; % number of filters to compute

% Compute iSTAC estimator
fprintf('\nComputing iSTAC estimate\n');
[sta,stc,rawmu,rawcov] = simpleSTC(Stim_tr,sps_tr,nkt);  % compute STA and STC
istacFilt = compiSTAC(sta(:),stc,rawmu,rawcov,nFilts); % find iSTAC filters

% Fit iSTAC exponentiated-quadratic nonlinearity using maximum likelihood
pp_istac = fitNlin_expquad_ML(Stim_tr,sps_tr,istacFilt,RefreshRate); 

%% == 3.  Set up struct and basis for LNP model   ========

% Set up fitting structure and compute initial logli
mask = [];  % time range to use for fitting (set to [] if not needed).
pp0 = makeFittingStruct_LNP(sta,RefreshRate,mask); % initialize param struct

% == Set up temporal basis for representing filters  ====
% (try changing these params until basis can accurately represent STA).
ktbasprs.neye = 0; % number of "identity"-like basis vectors
ktbasprs.ncos = 6; % number of raised cosine basis vectors (DETERMINES BASIS DIMENSIONALITY)
ktbasprs.kpeaks = [0 3*nkt/4]; % location of 1st and last basis vector bump
ktbasprs.b = 50; % determines how nonlinearly to stretch basis (higher => more linear)
[ktbas, ktbasis] = makeBasis_StimKernel(ktbasprs, nkt); % make basis
filtprs_basis = (ktbas'*ktbas)\(ktbas'*sta);  % filter represented in new basis
sta_basis = ktbas*filtprs_basis;

% Plot STA vs. best reconstruction in temporal basis
tt = (-nkt+1:0);  % time points in Stim_tr filter
subplot(221); % ----
plot(tt,ktbasis); xlabel('time bin'); title('temporal basis'); axis tight;
subplot(223); % ----
plot(tt,sta,tt,sta_basis,'r--', 'linewidth', 2); 
axis tight; title('STA and basis fit');
xlabel('time bin'); legend('sta', 'basis fit');
subplot(222); imagesc(sta); axis image; title('sta');
subplot(224); imagesc(sta_basis); axis image; title('basis fit');

% Compute R^2 for the STA in this basis
fprintf('R^2 = %.3f\n', 1-sum((sta(:)-sta_basis(:)).^2)/sum((sta(:)-mean(sta(:))).^2));

% Insert filter basis into fitting struct
pp0.k = sta_basis; % insert sta filter
pp0.kt = filtprs_basis; % filter coefficients (in temporal basis)
pp0.ktbas = ktbas; % temporal basis
pp0.ktbasprs = ktbasprs;  % parameters that define the temporal basis

%% == 4. Maximum likelihood estimation of filter under exponential nonlinearity [OPTIONAL]

negL0 = -logli_LNP(pp0,Stim_tr,sps_tr);  % negative log-likelihood at initial point
fprintf('\nFitting LNP model w/ exp nonlinearity\n');

% Do ML estimation of model params (with temporal basis defined in pp0)
opts = {'display', 'off', 'maxiter', 100};
[pp_exp,negLexp,Cexp] = fitLNP_1filt_ML(pp0,Stim_tr,sps_tr,opts); % find MLE by gradient ascent
eb1 = sqrt(diag(Cexp(1:nkt,1:nkt))); % 1SD error bars on filter (if desired)


%% == 5.  ML / MID estimation: filter and non-parametric nonlinearity with rbf basis ===

% Set parameters for radial basis functions (RBFs), for parametrizing nonlinearity
fstruct.nfuncs = 5; % number of RBFs (experiment with this)
fstruct.epprob = [.01, 0.99]; % cumulative probability outside outermost basis function peaks (endpoints)
fstruct.nloutfun = @logexp1;  % log(1+exp(x))  % nonlinear output function
fstruct.nlfuntype = 'cbf'; % for 1D nonlinearity, cbf is the same as rbf

% Initialize nonlinearity with filter fixed
fprintf('\nInitializing RBF nonlinearity\n');
[pp_cbf,negLrbf0] = fitNlin_CBFs(pp_exp,Stim_tr,sps_tr,fstruct);  % initialize nonlinearity while holding filter fixed

% Do maximum likelihood fit of filter and nonlinearity
fprintf('Jointly optimizing filter and RBF nonlinearity\n');
opts = {'display', 'off'}; % optimization parameters
[pp_cbf,negLcbf] = fitLNP_multifilts_cbfNlin(pp_cbf,Stim_tr,sps_tr,opts); % jointly fit filter and nonlinearity

% % In practice, can get both of the steps above with the single function call:
% [ppcbf,negLcbf] = fitLNP_multifiltsCBF_ML(pp0,Stim_tr,sps_tr,1,fstruct,pp0.k);


%% 6. ====== Make plots & report performance ============

% true 1st filter (from simulation)
uvec = @(x)(x(:)./norm(x(:))); % anonymous function to convert vector to unit vector
rs = @(x)(reshape(x,nkt,nx)); % reshape as image
trueK = rs(uvec(filts_true(:,1))); % true filter (as unit vector)

% -- Plot filter and filter estimates (as unit vectors) ---------
subplot(251);  imagesc(1:nx,tt,trueK); title('first model filter');
subplot(252);  imagesc(1:nx,tt,sta); title('sta');
subplot(253);  imagesc(1:nx,tt,rs(istacFilt)); title('istac');
subplot(254);  imagesc(1:nx,tt,pp_exp.k); title('ML-exp');
subplot(255);  imagesc(1:nx,tt,pp_cbf.k); title('ML-rbf');

% Note the filter labelled "true" is just the first filter of the true
% model (which had 5 filters) -- not clear that this filter is really most

% ---- Compute nonlinearities for plotting ---------
xnl = (xbinedges(1)-.1):.1:(xbinedges(end)+0.1); % x points for evaluating nonlinearity
ynl_hist = fnlhist(xnl); % histogram-based (piecewise constant) nonlinearity 
ynl_istac = pp_istac.nlfun(xnl*norm(pp_istac.k)); % istac exponentiated-quadratic nonlinearity
ynl_exp = exp(xnl*norm(pp_exp.k)+pp_exp.dc);  % exponential nonlinearity
ynl_cbf = pp_cbf.nlfun(xnl*norm(pp_cbf.k));   % rbf nonlinearity

% ---- Plot nonlinearities --------------------------
subplot(234); 
plot(xnl, ynl_hist,xnl,ynl_istac,xnl,ynl_exp,xnl,ynl_cbf, 'linewidth',2);
axis tight; set(gca,'ylim',[0 200]);
ylabel('rate (sps/s)'); xlabel('filter output');
legend('sta-hist','istac-expquad','ML-exptl','ML-cbf','location', 'northwest');
title('estimated nonlinearities');


%% 7. ==== Compute training and test performance in bits/spike =====

% Compute the log-likelihood under constant rate (homogeneous Poisson) model
muspike_tr = nsp_tr/slen_tr;       % mean number of spikes / bin, training set
muspike_tst = nsp_tst/slen_tst; % mean number of spikes / bin, test set
LL0_tr =   nsp_tr*log(muspike_tr) - slen_tr*muspike_tr; % log-likelihood, training data
LL0_tst = nsp_tst*log(muspike_tst) - slen_tst*muspike_tst; % log-likelihood test data

% A. Compute logli for lnp with histogram nonlinearity
pp_sta = pp0; % make struct for the sta+histogram-nonlinearity model
pp_sta.k = sta; % insert STA as filter
pp_sta.dc = 0; % remove DC component (if necessary)
pp_sta.kt = []; pp_sta.ktbas = []; % remove basis stuff (just to make sure it isn't used accidentally)
pp_sta.nlfun =  fnlhist;
LLsta_tr = logli_LNP(pp_sta,Stim_tr,sps_tr); % training log-likelihood
[LLsta_tst,rrsta_tst] = logli_LNP(pp_sta,Stim_tst,sps_tst); % test log-likelihood

% B. Compute logli for lnp with iSTAC (exponentiated-quadratic) nonlinearity
LListac_tr = logli_LNP(pp_istac,Stim_tr,sps_tr); % training log-likelihood
[LListac_tst,rristac_tst] = logli_LNP(pp_istac,Stim_tst,sps_tst); % test log-likelihood

% C. Compute logli for lnp with exponential nonlinearity
%ratepred_pGLM = exp(pGLMconst + Xdsgn*pGLMfilt); % rate under exp nonlinearity
LLexp_tr = logli_LNP(pp_exp,Stim_tr,sps_tr); % train
[LLexp_tst,rrexp_tst] = logli_LNP(pp_exp,Stim_tst,sps_tst); % test

% D. Compute logli for lnp with cbf/rbf nonlinearity
LLcbf_tr = logli_LNP(pp_cbf,Stim_tr,sps_tr); % train
[LLcbf_tst,rrcbf_tst] = logli_LNP(pp_cbf,Stim_tst,sps_tst); % test

% Single-spike information:
% ------------------------
% The difference of the loglikelihood and homogeneous-Poisson
% loglikelihood, normalized by the number of spikes, gives us an intuitive
% way to compare log-likelihoods in units of bits / spike.  This is a
% quantity known as the ("empirical" or "sample") single-spike information.
% [See Brenner et al, Neural Comp 2000; Williamson et al, PLoS CB 2015].

f1 = @(x)((x-LL0_tr)/nsp_tr/log(2)); % compute training single-spike info 
f2 = @(x)((x-LL0_tst)/nsp_tst/log(2)); % compute test single-spike info
% (if we don't divide by log 2 we get it in nats)

SSinfo_tr = [f1(LLsta_tr), f1(LListac_tr), f1(LLexp_tr), f1(LLcbf_tr)];
SSinfo_tst = [f2(LLsta_tst),f2(LListac_tst), f2(LLexp_tst), f2(LLcbf_tst)];

fprintf('\nSingle-spike information (bits/spike):\n');
fprintf('------------------------------------- \n');
fprintf('Train: sta-hist:%.2f  istac: %.2f  exp:%.2f  cbf:%.2f\n', SSinfo_tr);
fprintf('Test:  sta-hist:%.2f  istac: %.2f  exp:%.2f  cbf:%.2f\n', SSinfo_tst);

% Plot test single-spike information
subplot(235);
axlabels = {'sta','istac','exp','cbf'};
bar(SSinfo_tst); ylabel('bits / sp'); title('test single spike info');
set(gca,'xticklabel',axlabels, 'ylim', [min(SSinfo_tst*.9), max(SSinfo_tst)*1.05]);

% ==== Last: plot the rate predictions for the two models =========
subplot(236); 
iiplot = 1:100; % time bins to plot
stem(iiplot,sps_tst(iiplot), 'k'); hold on;
plot(iiplot,rrsta_tst(iiplot)/RefreshRate, ...
    iiplot,rristac_tst(iiplot)/RefreshRate, ...
    iiplot,rrexp_tst(iiplot)/RefreshRate, ...
    iiplot,rrcbf_tst(iiplot)/RefreshRate,'linewidth',2); 
 hold off; title('rate predictions on test data');
ylabel('spikes / bin'); xlabel('time (bins)');
set(gca,'xlim', iiplot([1 end]));
legend('spike count', 'sta-hist', 'istac','ML-exp', 'ML-cbf');
