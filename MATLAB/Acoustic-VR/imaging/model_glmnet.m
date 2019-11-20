%calium imaging model using glmnet
nd = 6; %number of days
nc = 208; %number of cells
nl = 0; %number of lags
opts = glmnetSet();
opts.intr = 1; opts.standardize = 1; %don't include intercept or standardize (I do that)
opts.alpha = 0;
F = cell(1,nd); %flourescence signals for each day
S = cell(1,nd); %spikes detected for each day
lags = cell(1,nd); %lagged design matrices
C = zeros(nc,nc,nl+1,nd); %cross-coefficient matrix for each lag and day
R = zeros(nc,nd,nd); % r-squared for each cell/day combo
expgrid = linspace(0,0.1*40,40);

%construct design matrix for each day
for d = 1:nd    
   %z-score flourescence traces 
   F{d} = zscore(datasets{d}.dF);   
   S{d} = zeros(size(F{d}));
   %compute and concatenate all lags
   lags{d} = zeros(size(F{d},1)-nl,(nl+1)*nc);
   for c = 1:nc      
      sm = adsmo(F{d}(:,c),1:size(F{d},1),1,100); 
      sm = adsmo(sm,1:size(F{d},1),1,50); %smooth flourescence traces
      [pks, locs] = findpeaks(-sm,'Threshold',0.5); %remove negative spikes
      for l=1:length(locs), sm(locs(l))=mean([sm(locs(l)-1) sm(locs(l)+1)]);end
%       [dconv, ~] = deconv(sm,exp(-expgrid)); %exponential deconvolution
%       s = std(dconv(dconv<0));
%       S{d}(dconv>2*s,c) = 1; %spike detection
      S{d}(:,c) = sm;
      l = phasespace(S{d}(:,c),nl+1,1); 
      lags{d}(:,1+(nl+1)*(c-1):(nl+1)*c) = fliplr(l);      
   end    
end

%fit glmnet on each day
for d = 1:nd   
   %run glmnet for each cell
   for c = 1:nc
      
       notc = 1:size(lags{d},2);
       notc(1+(nl+1)*(c-1):c*(nl+1)) = []; %cell indices excluding current cell 
       X = [lags{d}(:,notc) sum(lags{d}(:,notc),2)];
       cvfit = cvglmnet(X,lags{d}(:,1+(nl+1)*(c-1)),'gaussian',opts);  %fit model     

       %put coefficients into cross-coefficient matrix
       coef = cvglmnetCoef(cvfit,'lambda_min'); %model coefficients
       coef(1) = [];
       for l = 1:nl+1     
           notc = 1:nc+1;
           notc(c) = [];
           C(c,notc,l,d) = coef(l:(nl+1):end);
       end
       
       %get r-squared for each day
       for dd = 1:nd
           notc = 1:size(lags{d},2);
           notc(1+(nl+1)*(c-1):c*(nl+1)) = []; %cell indices excluding current cell
           X = [lags{dd}(:,notc) sum(lags{dd}(:,notc),2)];
           pred = cvglmnetPredict(cvfit,X,'lambda_min');
%            dev = S{dd}(nl+1:end,c).*log(S{dd}(nl+1:end,c)./exp(pred))-S{dd}(nl+1:end,c)+exp(pred);
%            dev(isnan(dev)) = exp(pred(isnan(dev)));
%            dev = 2*sum(dev);
%            mu = mean(S{dd}(nl+1:end,c))*ones(size(S{dd}(nl+1:end,c)));
%            ndev = S{dd}(nl+1:end,c).*log(S{dd}(nl+1:end,c)./mu)-S{dd}(nl+1:end,c)+mu;      
%            ndev(isnan(ndev)) = mu;
%            ndev = 2*sum(ndev);
           R(c,d,dd) = 1 - sum((S{dd}(nl+1:end,c)-pred).^2)/sum((S{dd}(nl+1:end,c)).^2);
%            R(c,d,dd) = (ndev-dev)/ndev;
       end
       
   end
   disp(sprintf('%s%s',num2str(d/nd*100),'% Complete'));
end

