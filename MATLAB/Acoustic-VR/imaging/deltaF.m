function [dF, nF, zF] = deltaF(raw, method)
%UNTITLED2 Summary of this function goes here
%   raw- rows= frames, columns = rois

IMG_FS=10.088;

remove_nan_rois=0;
remove_nan_points=0;

% clean NaNs
% set NaN to mean for that ROI
if (remove_nan_points==1)
[r,c]=find(isnan(raw)); % r = frame, c= roi
for i=1:length(r); raw(r(i),c(i))=mean(raw(:,c(i)), 'omitnan');end 
end

% clean blank ROIs, 01.04.2017, happens when ROIbuddy skips roi labels
if (remove_nan_rois==1)
    b=raw';
    nan_ind=~any(isnan(b),2)';
    raw=raw(:,nan_ind);
    labels=labels(:,nan_ind);
    tags=tags(nan_ind',:);
    
end

[num_frames, num_rois]=size(raw);

switch method 
    case 1
        % Bl = median of whole trace
        % if nan not omitted, the whole df/f becomes nan. not sure if I need to
        % clean nan values in subsequent step, now that I omit
        baseline=median(raw,1,'omitnan');

        %clean NaNs
        [r,c]=find(isnan(raw));
        for i=1:length(r); raw(r(i),c(i))=baseline(c(i));end 

        %finish d/F
        bl=repmat(baseline,size(raw,1),1);
        dF=((raw-bl)./bl)*100;


        % for trials, can put if statement hear to make it work with eb trials
        %longF=reshape(dF,num_frames,num_rois); % dont need for cont
    case 2
        
        baseline=[]; %not need for this case, but need blank for struct below
        %wSize= 3; % s, window size
        %win_size= round(100*wSize*1/IMG_FS); % # frames to grab
        win_size=50;
        
        dF=zeros(num_frames, num_rois);
        nF=zeros(num_frames, num_rois);
        for i=1:length(raw)
            if (i<=win_size)
                s=i+1; e=i+50;
                v=raw(s:e,:);
                sv=sort(v);
                f0=mean(sv(1:25,:),1);

                dF(i,:)=(raw(i,:)-f0)./f0;
            else
                s=i-50; e=i-1;
                v=raw(s:e,:);
                sv=sort(v);
                f0=mean(sv(1:25,:),1);

                dF(i,:)=(raw(i,:)-f0)./f0;
            end
        end
        
        
        % x = (dF-min(dF))/(max(dF)-min(dF));
        
        
    case 3
        % standard for EB like trials
        
        baseline_frames=6:25;

        baseline=mean(raw(baseline_frames,:,:));
        bl=repmat(baseline,num_frames,1);
        dF=((raw-bl)./bl)*100;
        
    case 4
        %robust lowess (robust to outliers i.e. calcium transients)
        %rlowess is a little faster than rloess and more stable at boundaries (a little less wiggly overall)
        %both much slower than methods 1-3, but more statistically standard (~20 min for 200 cells, 6000 frames)        
        baseline = zeros(size(raw));
        for c = 1:num_rois           
            baseline(:,c) = smooth(raw(:,c),1/15,'rlowess');            
        end
        dF = (raw - baseline)./baseline;
end

%normalize by largest negative transient, so everything >1 is probably real
%also z-score using std of negative transients
nF = dF;
zF = dF;
for c = 1:size(dF,2)    
    nF(:,c) = nF(:,c)/max(abs(nF(nF(:,c)<0,c)));
    zF(:,c) = zF(:,c)/abs(mean(zF(zF(:,c)<0,c)));
end



