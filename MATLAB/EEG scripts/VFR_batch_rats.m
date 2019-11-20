function [VFR] = VFR_batch(inFile,freqVec,width)

img_text = 'on';

fid = fopen(inFile,'r');
if fid == -1
    msgbox('Could not open the input file! Make sure the filname and path are correct.','ERROR');
end

% Get sessions and csc-file list from input file
fid = fopen(inFile,'r');
ii = -1;     
numsessions = 0;
while ~feof(fid)
    str = fgetl(fid);
    if ii == 0
%         str(3:7) = [];
%         str(1) = 'G';
        cscList = str;
    elseif ii == 1
        refList = str;
    elseif ii == 2
        ch4avg = str;
    elseif ii > 2
        numsessions  = numsessions+1;
        if ~strcmp(str(end),'\')
            str = strcat(str,'\');
        end
%         str(3:7) = [];
%         str(1) = 'G';
        sessions(numsessions) = {str};
    end
    ii = ii+1;
end

% read the file names from the csc-file list
cscid = fopen(cscList,'r');
jj = 1;
while ~feof(cscid)
       str = fgetl(cscid);
       channels(jj) = {str};
       jj = jj+1;
end
numchannels = jj-1;
cscid = fclose('all');

% Set the field selection for reading the video files. 1 = Add parameter, 0 = skip
% parameter
fieldSelection(1) = 1; % Timestamps
fieldSelection(2) = 1; % Extracted X
fieldSelection(3) = 1; % Extracted Y
fieldSelection(4) = 0; % Extracted Angel
fieldSelection(5) = 1; % Targets
fieldSelection(6) = 0; % Points
% Do we return header 1 = Yes, 0 = No.
extractHeader = 0;
% 5 different extraction modes, see help file for Nlx2MatVt
extractMode = 1; % Extract all data

V = [];
D = 90; %90 for circular track, 200 for linear

for ii = 1:numsessions
    disp(sprintf('%s%s','Reading data for session: ',sessions{ii}));
    % Check if subdir for storing images are present. If not, it is
    % created
    dirInfo = dir(sessions{ii});
    found = 0;
    for kk=1:size(dirInfo,1)
        if dirInfo(kk).isdir
            if strcmp(dirInfo(kk).name,strcat('VFR','\'))
                found = 1;
            end
        end
    end
    if found==0
        mkdir(cd,strcat('VFR','\'));
    end
    
    %compute running speed
    % Get position data
    file = strcat(sessions{ii},'vt1.nvt');
    [t{ii},x,y, targets] = Nlx2MatVT(file,fieldSelection,extractHeader,extractMode);
    % Convert timestamps to seconds
    t{ii} = t{ii}'/1000000;
    
    %interpolate (consider polar coordinates for circle track)
    x = pchip(t{ii}(x~=0),x(x~=0),t{ii});
    y = pchip(t{ii}(y~=0),y(y~=0),t{ii});
    
%     %(for linear track)
%     %linearize and scale
%     X = [ones(length(x),1) x];
%     b = X\y;
%     y = y - b(1);
%     alpha = atan2(b(2),1);    
%     Q = [cos(-alpha), -sin(-alpha); sin(-alpha) cos(-alpha)];
%     postemp = Q*[x';y'];
%     x = postemp(1,:);
%     y = postemp(2,:);
%     x = x - min(x);
%     q1 = quantile(x,0.025);
% %     q2 = quantile(x-q1,0.975);
%     pos = (x-q1)/q2*D;
    
    %(for circle track)
    %center track based on interior of trajectory 
    %start by plotting trajectory and grabbing image
    figure(ii)
    plot(x,y,'-k');axis equal;axis off
    im = frame2im(getframe(gcf));
    %binarize image
    bw = imbinarize(im(:,:,1));
    %get interior blob
    blobs = bwconncomp(bw,4);
    numPixels = cellfun(@numel,blobs.PixelIdxList);
    [~,idx] = sort(numPixels);    
    intblob = zeros(size(bw));
    intblob(blobs.PixelIdxList{idx(end-1)}) = 1;
    [inty,intx] = find(intblob);
    [bwy,bwx] = find(bw==0);
    %get blob center and convert back to x, y
    intxc = (mean([max(intx) min(intx)])-min(bwx))/(max(bwx)-min(bwx));
    intyc = (mean([max(inty) min(inty)])-min(bwy))/(max(bwy)-min(bwy));
    intRx = (max(intx)-min(intx))/(max(bwx)-min(bwx))/2;
    intRy = (max(inty)-min(inty))/(max(bwy)-min(bwy))/2;
    xc = min(x) + intxc*(max(x)-min(x));
    yc = max(y) - intyc*(max(y)-min(y));
    Rx = intRx*(max(x)-min(x));
    Ry = intRy*(max(y)-min(y));
    %center and scale trajectory to perfect circle
    x = (x - xc)/Rx;
    y = (y - yc)/Ry;
    figure(ii)
    plot(x,y,'-k');axis equal;hold on
    viscircles([0 0], 1,'EdgeColor','r');hold off
    im = getframe(gcf); %forces matlab to display the plot onscreen
    %position in unwrwapped cm
    pos = unwrap(atan2(y,x))*D/2;
    
    %smooth and get velocity   
    span = 1/(max(t{ii})-min(t{ii}));
    pos = smooth(pos,span,'loess');
    vel = abs(diff(pos))./diff(t{ii});
    vel = [vel(1);vel];
    vel = smooth(vel,span,'loess');
%     vel = exp(smooth(real(log(vel)),span,'lowess'));
    t{ii}(isnan(vel)) = [];
    vel(isnan(vel)) = [];
    V = [V;vel];    
end
    
% Load data from the .ncs files, make plots, and store them
for jj = 1:numchannels
    TFR = [];
    for ii = 1:numsessions
        xfile = [sessions{ii},channels{jj}];
        [x,ts,tt, Fs, bv, ir] = loadEEG2(xfile);  
        if strcmp(sessions{ii},'G:\Rats\rat113\2017-01-30_combine\begin2\')
            rfile = 'G:\Rats\rat113\2017-01-30_combine\begin2\Teth2_ref.ncs'
            [r,ts,tt, Fs, bv, ir] = loadEEG2(rfile); 
            x = x - r;
        end
        %get TFR indices for one sample per t{ii}
        if jj == 1
            ttidx{ii} = zeros(size(t{ii}));
            for it = 1:length(t{ii})
                [~,idx] = min((tt-t{ii}(it)).^2);
                ttidx{ii}(it) = idx;
            end
        end        
        [tfr] = traces2TFR_rev(x,freqVec,Fs,width);
        tfr = tfr(:,ttidx{ii});  
        tfr = zscore(tfr,0,2);
        TFR = [TFR tfr];
    end
    
    minv = 0;maxv = quantile(V,0.98);
    vbins = minv:1:maxv;
    sigma = 5;
    delta = repmat(V,1,length(vbins))-repmat(vbins,length(V),1);
    W = exp(-0.5*delta.*delta/sigma^2);
    D = ones(length(freqVec),length(V))*exp(-0.5*delta.*delta/sigma^2);
    VFR = TFR*W./D;
    figure(1);
    col = max(abs([min(min(VFR(freqVec>25,1:round(0.5*length(vbins))))) max(max(VFR(freqVec>25,1:round(0.5*length(vbins)))))]));
    surf(vbins,freqVec,zeros(size(VFR)),VFR);view([0 90]);axis tight;axis square;shading interp;
        caxis([-col col]);colorbar
        set(gca,'Ytick',0:10:250,'yticklabel',0:10:250,'YScale','log','XScale','linear');
        xlabel('running speed (cm/sec)');ylabel('Frequency (Hz)');
        title(strcat('Velocity Frequency Representation: ',channels{jj}));
    bmpImage = sprintf('%s%s%s%s',strcat('VFR','\'),channels{jj}(1:end-4),'.bmp');
    saveas(gcf,bmpImage,'bmp');
    
end