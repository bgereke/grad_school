% Set parameters
dtStim = 1; % Bin size for simulating model & computing likelihood (in units of stimulus frames)
dtSp = 1;  % Bin size for simulating model & computing likelihood (must evenly divide dtStim);
nkt = 500;    % Number of time bins in stimulus filter k
Stim = zscore(vel');
Stim = zeros(size(Stim));

% Compute Spike-triggered averages
sps1 = ridges(:,37);
sta1 = simpleSTC(Stim,sps1,nkt); % Compute STA 1
sta1 = reshape(sta1,nkt,[]); 

% Initialize param struct for fitting 
gg0 = makeFittingStruct_GLM(dtStim,dtSp);  % Initialize params for fitting struct 

% Initialize fields (using h and k bases computed above)
ggsim2 = makeSimStruct_GLM(nkt,dtStim,dtSp);  % Create GLM struct with default params
ggsim2.ktbasprs.kpeaks = [0 nkt];
ggsim2.ktbasprs.neye = 0;
ggsim2.ktbasprs.ncos = 12;
ggsim2.ktbasprs.b = 1000;
[ktbas,ktbasis] = makeBasis_StimKernel(ggsim2.ktbasprs,nkt);
knots = linspace(1,nkt,ggsim2.ktbasprs.ncos);
ktbas = fnval(csaps(knots,eye(length(knots)),1),1:nkt)';
gg0.ktbas = ktbas; % k basis

%cross-history bases
ihbasprs.ncols = 12;  % Number of basis vectors for post-spike kernel
ihbasprs.hpeaks = [dtSp 0.75*nkt*dtSp];  % Peak location for first and last vectors
ihbasprs.b = dtSp*5;  % Determines how nonlinear to make spacings
ihbasprs.absref = []; % absolute refractory period (optional)
% [iht,ihbas,ihbasis] = makeBasis_PostSpike(ihbasprs,dtSp);
knots = linspace(1,nkt,ihbasprs.ncols);
ihbas = fnval(csaps(knots,eye(length(knots)),1),1:nkt)';
gg0.ihbas = ihbas; % h self-coupling basis
gg0.ihbas2 = ihbas; % h coupling-filter basis
nktbasis = size(ktbas,2); % number of basis vectors in k basis
nhbasis = size(ihbas,2); % number of basis vectors in h basis
gg0.kt = (ktbas\sta1); % initial params from scaled-down sta 
gg0.k = gg0.ktbas*gg0.kt;  % initial setting of k filter
% gg0.k = gg0.k/norm(gg0.k);
gg0.ihw = zeros(nhbasis,1); % params for self-coupling filter
% gg0.ihw2 = zeros(nhbasis,97); % params for cross-coupling filter
% gg0.ih = [gg0.ihbas*gg0.ihw gg0.ihbas2*gg0.ihw2];
gg0.ih = [gg0.ihbas*gg0.ihw];
gg0.iht = iht;
gg0.dc = 0; % Initialize dc term to zero
gg0.couplednums = []; % number of cell coupled to this one (for clarity)

% Set spike responses for cell 1 and coupled cell
gg0.sps = ridges(:,37);  
% gg0.sps2 = ridges;gg0.sps2(:,[1 37 100]) = [];

% Compute initial value of negative log-likelihood (just to inspect)
[neglogli0,rr] = neglogli_GLM(gg0,Stim);

% Do ML fitting
fprintf('Fitting neuron 1:  initial neglogli0 = %.3f\n', neglogli0);
opts = {'display', 'iter', 'maxiter', 100};
[gg1, neglogli1] = MLfit_GLM(gg0,Stim,opts); % do ML (requires optimization toolbox)