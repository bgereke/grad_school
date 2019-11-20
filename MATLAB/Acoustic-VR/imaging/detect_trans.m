function [dT,nT,zT] = detect_trans(dF,nF,zF)

[numframes, numcells] = size(dF);
wname  = 'sym4';               % Wavelet for analysis.
level  = 5;     

%compute dT
dT = zeros(size(dF)); %transients
% expgrid = linspace(0,0.032979438*120,120);
for c = 1:numcells
%     sm = adsmo(dF(:,c),1:numframes,1,50);
%     sm = adsmo(sm,1:numframes,1,25); %smooth flourescence traces
    [sm,~,~] = cmddenoise(dF(:,c),wname,level); %wavelet denoising
%     [~, locs] = findpeaks(-sm,'Threshold',0.5); %remove negative spikes
%     for l=1:length(locs), sm(locs(l))=mean([sm(locs(l)-1) sm(locs(l)+1)]);end
%     [dconv, ~] = deconv(sm,exp(-expgrid)); %exponential deconvolution
%     s = abs(min(dconv));
    s = 1 - min(sm-1);
%     dT(dconv>s,c) = 1; %spike detection
    dT(sm>s,c) = 1; %spike detection
end

%comput nT
nT = zeros(size(dF));
hthresh = nF>1; %upper threshold
lthresh = nF>0.5; %lower threshold
dhthresh = [zeros(1,numcells);diff(hthresh)];
dlthresh = [zeros(1,numcells);diff(lthresh)];
for c = 1:numcells
    hidx = find(dhthresh(:,c)==1); %high threshold crossings
    lup = find(dlthresh(:,c)==1); %low threshold upward crossings
    ldown = find(dlthresh(:,c)==-1); %low threshold downward crossings
    for h = 1:length(hidx)
        sidx = lup(find(lup<hidx(h),1,'last'));
        eidx = ldown(find(ldown>hidx(h),1,'first'));
        if ~isempty(sidx)>0 && ~isempty(eidx)>0
            w = eidx-sidx;
            if w >= 4 %delete transients that aren't wide enough
                nT(sidx,c) = 1;
            end
        end
    end
end
    
%compute zT
zT = zeros(size(dF));
hthresh = zF>3; %upper threshold
lthresh = zF>1; %lower threshold
dhthresh = [zeros(1,numcells);diff(hthresh)];
dlthresh = [zeros(1,numcells);diff(lthresh)];
for c = 1:numcells
    hidx = find(dhthresh(:,c)==1); %high threshold crossings
    lup = find(dlthresh(:,c)==1); %low threshold upward crossings
    ldown = find(dlthresh(:,c)==-1); %low threshold downward crossings
    for h = 1:length(hidx)
        sidx = lup(find(lup<hidx(h),1,'last'));
        eidx = ldown(find(ldown>hidx(h),1,'first'));
        if ~isempty(sidx)>0 && ~isempty(eidx)>0
            w = eidx-sidx;
            if w >= 4 %delete transients that aren't wide enough
                zT(sidx,c) = 1;
            end
        end
    end
end