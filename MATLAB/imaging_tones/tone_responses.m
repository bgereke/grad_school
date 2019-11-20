% tone responses
expgrid = linspace(0,0.1*40,40);
smraw = zeros(size(raw));
dconv = zeros(size(raw));dconv(end-38:end,:) = [];
spks = zeros(size(dconv));
nc = size(raw,2);
nt = size(raw,1);

for c = 1:nc     
      sm = adsmo(dF(:,c),1:nt,1,50); 
      sm = adsmo(sm,1:nt,1,25); %smooth flourescence traces
      [pks, locs] = findpeaks(-sm,'Threshold',0.5); %remove negative spikes
      for l=1:length(locs), sm(locs(l))=mean([sm(locs(l)-1) sm(locs(l)+1)]);end
      smraw(:,c) = sm;
      [dconv(:,c), ~] = deconv(sm,exp(-expgrid)); %exponential deconvolution
      s = std(smraw(:,c));
      spks(smraw(:,c)>2*s,c) = 1; %spike detection     
end 

ns = 10000;
sample = zeros(ns,15);

for s = 1:ns
    
    sam = randsample(modpos(:,1),797,true);
    sample(s,:) = hist(sam,15);
    
end