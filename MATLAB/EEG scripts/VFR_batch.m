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
fieldSelection(5) = 0; % Targets
fieldSelection(6) = 0; % Points
% Do we return header 1 = Yes, 0 = No.
extractHeader = 0;
% 5 different extraction modes, see help file for Nlx2MatVt
extractMode = 1; % Extract all data

V = [];

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
    [t{ii},x,y] = Nlx2MatVT(file,fieldSelection,extractHeader,extractMode);
    % Convert timestamps to seconds
    t{ii} = t{ii}'/1000000;
%     fidx = x==0&y==0; 
    
    %interpolate positions
    x = pchip(t{ii}(x~=0),x(x~=0),t{ii});
    y = pchip(t{ii}(y~=0),y(y~=0),t{ii});
    
    %center track
    nb = 250;
    xymap = zeros(nb,nb);
    v = 0.1;
    mx = linspace(min(x),max(x),nb);
    my = linspace(min(y),max(y),nb);    
    for xx = 1:length(mx)
        dx = mx(xx)-x;
        for yy = 1:length(my)
            dy = my(yy)-y;
            xymap(xx,yy) = nansum(v^2/sqrt(2*pi)*exp(-0.5*(dx.*dx*v^2+dy.*dy*v^2)));
        end
    end
    thmap = xymap+0.1;thmap(xymap>quantile(xymap(:),0.65))=0;
    
    regmax = imregionalmax(xymap,4);
    [xx,yy] = find(regmax);
    [z, a, b, alpha] = fitellipse([mx(xx);my(yy)]); %'linear' to speed up
    %Translate
    x = x-z(1);
    y = y-z(2);
    %Rotate
    Q = [cos(-alpha), -sin(-alpha); sin(-alpha) cos(-alpha)];
    temp = Q*[x';y'];
    %Scale
    D = 100;
    temp(1,:) = 0.5*D/a*temp(1,:);
    temp(2,:) = 0.5*D/b*temp(2,:);
    %Rotate back to orginal orientation
    Q = [cos(alpha), -sin(alpha); sin(alpha) cos(alpha)];
    temp = Q*temp;
    x = temp(1,:);y = temp(2,:);    

    [xc yc R] = circfit(x,y);
    x = x - xc;
    y = y - yc;
%     badidx = find(diff(t{ii})<0)+1;
%     x(badidx)=[];y(badidx)=[];
%     t{ii}(badidx)=[];
%     D=200;
%     X = [ones(length(x),1) x];
%     b = X\y;
%     y = y - b(1);
%     alpha = atan2(b(2),1);    
%     Q = [cos(-alpha), -sin(-alpha); sin(-alpha) cos(-alpha)];
%     postemp = Q*[x';y'];
%     x = postemp(1,:);
%     pos = (x-min(x))/quantile((x-min(x)),0.95)*D;
    
    %convert position to phase, smooth, and get velocity
    D = 100;
    pos = -unwrap(atan2(y,x));    
    span = 3/(max(t{ii})-min(t{ii}));
    pos = smooth(pos,span,'loess');
    vel = abs(diff(pos))./diff(t{ii});
    vel = [vel(1);vel];
    span = 3/(max(t{ii})-min(t{ii}));
    vel = exp(smooth(log(vel),span,'loess'));
%     vel = exp(smooth(real(log(vel)),span,'lowess'))*D/2;
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
        x = x*bv;
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
       tfr = abs(traces2Wx(x,freqVec,Fs,width));
        tfr = tfr(:,ttidx{ii}); 
%         tfr = zscore(tfr,0,2);
        TFR = [TFR tfr];
    end
    
    V = V/(2*pi)*pi*100;
    minv = 0;maxv = quantile(V,0.98);
    vbins = linspace(minv,maxv,50);
    sigma = 0.5;
    delta = repmat(V,1,length(vbins))-repmat(vbins,length(V),1);
    W = exp(-0.5*delta.*delta/sigma^2);
    D = ones(length(freqVec),length(V))*exp(-0.5*delta.*delta/sigma^2);
    VFR = TFR*W./D;keyboard
    figure(1);
    col = 0.7*max(max(VFR));
    surf(vbins,freqVec,zeros(size(VFR)),VFR);view([0 90]);axis tight;axis square;shading interp;
        caxis([0 col]);colorbar
        set(gca,'Ytick',0:10:250,'yticklabel',0:10:250,'YScale','linear','XScale','linear');
        xlabel('running speed (cm/sec)');ylabel('Frequency (Hz)');
        title(strcat('Velocity Frequency Representation: ',channels{jj}));colormap hot
        keyboard
    bmpImage = sprintf('%s%s%s%s',strcat('VFR','\'),channels{jj}(1:end-4),'.bmp');
    saveas(gcf,bmpImage,'bmp');
    
end