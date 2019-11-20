% CFC('inFile.txt')
%
% CFC generates cross frequency coherences, plots them and stores them to images and .mat files.
% Multiple sessions will be read from the CSC's specififed in
% 'CSCList.txt'. 'TTList.txt' is necessary for spike data scripts that use
% the same 'infile.txt'.
%
% The input file must be on the following format.
%
% C:\Data\TTList.txt
% C:\Data\CSCList.txt
% C:\Data\Begin 1
% C:\Data\Begin 2
% C:\Data\Begin 3
% C:\Data\Begin 4
% and so on ...
%
% 'CSCList.txt' contains a list of the Neuralynx .csc files to be analyzed.
% All plots will be stored to both bmp and eps imagefiles to a subdirectory in
% the data folder called CFC_plots.

function [CS] = traces2CCS_biv_rats(inFile,freqVec,width)

img_text = 'on';

fid = fopen(inFile,'r');
if fid == -1
    msgbox('Could not open the input file! Make sure the filname and path are correct.','ERROR');
end

parent = cd;

% Get sessions and csc-file list from input file
fid = fopen(inFile,'r');
ii = -1;
numsessions = 0;
while ~feof(fid)
    str = fgetl(fid);
    if ii == 0
        str(3:7) = [];
        str(1) = 'G';
        str(end-11) = '1';
        ca1List = str;
        str(end-11) = '3';
        ca3List = str;
    elseif ii == 1
        refList = str;
    elseif ii == 2
        ch4avg = str;
    elseif ii > 2
        numsessions  = numsessions+1;
        if ~strcmp(str(end),'\')
            str = strcat(str,'\');
        end
        str(3:7) = [];
        str(1) = 'G';
        sessions(numsessions) = {str};
    end
    ii = ii+1;
end

% read the file names from the csc-file list
cscid = fopen(ca1List,'r');
jj = 1;
while ~feof(cscid)
    str = fgetl(cscid);
    if str == -1
        continue
    end
    ca1channels(jj) = {str};
    jj = jj+1;
end
numCA1 = jj-1;
cscid = fclose('all');

if numCA1 == 0
    return
end

cscid = fopen(ca3List,'r');
jj = 1;
while ~feof(cscid)
    str = fgetl(cscid);
    if str == -1
        continue
    end
    ca3channels(jj) = {str};
    jj = jj+1;
end
numCA3 = jj-1;
cscid = fclose('all');

if numCA3 == 0
    return
end

% Set the field selection for reading the video files. 1 = Add parameter, 0 = skip
% parameter
fieldSelection(1) = 1; % Timestamps
fieldSelection(2) = 1; % Extracted X
fieldSelection(3) = 1; % Extracted Y
fieldSelection(4) = 0; % Extracted Angle
fieldSelection(5) = 1; % Targets
fieldSelection(6) = 0; % Points
% Do we return header 1 = Yes, 0 = No.
extractHeader = 0;
% 5 different extraction modes, see help file for Nlx2MatVt
extractMode = 1; % Extract all data

D = 200;
t = cell(numsessions,1);
tt = cell(1,numsessions);
phase = cell(numsessions,1);
Wx = cell(1,numsessions);
Wy = cell(1,numsessions);
sidx = cell(1,numsessions); %session indices

for ii = 1:numsessions
    disp(sprintf('%s%s','Reading data for session: ',sessions{ii}));
    
    % Check if subdir for storing images are present. If not, it is
    % created
    dirInfo = dir(sessions{ii});
    found = 0;
    for kk=1:size(dirInfo,1)
        if dirInfo(kk).isdir
            if strcmp(dirInfo(kk).name,strcat('CS','\'))
                found = 1;
            end
        end
    end
    if found==0
        mkdir(sessions{ii},strcat('CS','\'));
    end
    
    %compute running speed
    % Get position data
    file = strcat(sessions{ii},'vt1.nvt');
    [t{ii},x,y, targets] = Nlx2MatVT(file,fieldSelection,extractHeader,extractMode);
    % Convert timestamps to seconds
    t{ii} = t{ii}'/1000000;
    
    %interpolate don't extrapolate
    start = find(x~=0&y~=0,1,'first');
    stop = find(x~=0&y~=0,1,'last');
    xx = pchip(t{ii}(x~=0&y~=0),x(x~=0&y~=0),t{ii}(start:stop));
    yy = pchip(t{ii}(x~=0&y~=0),y(x~=0&y~=0),t{ii}(start:stop));
    x = xx;y = yy; t{ii} = t{ii}(start:stop);
    
    %(for linear track)
    %linearize and scale
    X = [ones(length(x),1) x];
    b = X\y;
    y = y - b(1);
    alpha = atan2(b(2),1);
    Q = [cos(-alpha), -sin(-alpha); sin(-alpha) cos(-alpha)];
    postemp = Q*[x';y'];
    x = postemp(1,:);
    y = postemp(2,:);
    x = x - min(x);
    q1 = quantile(x,0.025);
    q2 = quantile(x-q1,0.975);
    pos = (x-q1)/q2*D;
    phase{ii} = pos';
    
    %     %(for circle track)
    %     figure(ii)
    %     plot(x,y,'-k');axis equal;axis off
    %     im = frame2im(getframe(gcf));
    %     %binarize image
    %     bw = imbinarize(im(:,:,1));
    %     %get interior blob
    %     blobs = bwconncomp(bw,4);
    %     numPixels = cellfun(@numel,blobs.PixelIdxList);
    %     [~,idx] = sort(numPixels);
    %     intblob = zeros(size(bw));
    %     intblob(blobs.PixelIdxList{idx(end-1)}) = 1;
    %     [inty,intx] = find(intblob);
    %     [bwy,bwx] = find(bw==0);
    %     %get blob center and convert back to x, y
    %     intxc = (mean([max(intx) min(intx)])-min(bwx))/(max(bwx)-min(bwx));
    %     intyc = (mean([max(inty) min(inty)])-min(bwy))/(max(bwy)-min(bwy));
    %     intRx = (max(intx)-min(intx))/(max(bwx)-min(bwx))/2;
    %     intRy = (max(inty)-min(inty))/(max(bwy)-min(bwy))/2;
    %     xc = min(x) + intxc*(max(x)-min(x));
    %     yc = max(y) - intyc*(max(y)-min(y));
    %     Rx = intRx*(max(x)-min(x));
    %     Ry = intRy*(max(y)-min(y));
    %     %center and scale trajectory to perfect circle
    %     x = (x - xc)/Rx;
    %     y = (y - yc)/Ry;
    %     figure(ii)
    %     plot(x,y,'-k');axis equal;hold on
    %     viscircles([0 0], 1,'EdgeColor','r');hold off
    %     im = getframe(gcf); %forces matlab to display the plot onscreen
    %     saveas(gcf,strcat(sessions{ii},'trajectory.bmp'),'bmp');
    %     %position in rad
    %     phase{ii} = atan2(y,x);
    
    %load wavelet transforms for this session
    Wx{ii} = zeros(length(freqVec),length(t{ii}),2);
    Wy{ii} = zeros(length(freqVec),length(t{ii}),2);

    xfile = [sessions{ii},ca1channels{1}];
    [x,~,tt{ii}, Fs, bv, ~] = loadEEG2(xfile);
    x = bv*x; clear samples

    yfile = [sessions{ii},ca3channels{1}];
    [y,~,~, Fs, bv, ~] = loadEEG2(yfile);
    y = bv*y; clear samples
    
    %remove reference
    rfile = [sessions{ii},'Teth1_ref.ncs'];
    [r,~,~, Fs, bv, ~] = loadEEG2(rfile);
    r = bv*r; clear samples
    x = x + r; 
    y = y + r;
    
    wx = traces2Wx(x,freqVec,Fs,width);
    wx(:,tt{ii} < min(t{ii}) | tt{ii} > max(t{ii})) = [];
    keyboard
    wy = traces2Wx(y,freqVec,Fs,width);
    wy(:,tt{ii} < min(t{ii}) | tt{ii} > max(t{ii})) = [];

    ttidx = zeros(length(t{ii}),1);
    tt{ii}(:,tt{ii} < min(t{ii}) | tt{ii} > max(t{ii})) = [];
    for it = 1:length(t{ii})
        [~,ttidx(it)] = min((tt{ii}-t{ii}(it)).^2);
    end
    
    Wx{ii} = wx(:,ttidx);
    Wy{ii} = wy(:,ttidx);
    tt{ii} = tt{ii}(ttidx);
    clear wx wy
    
    %get session indices
    mt = min(t{ii});
    t{ii} = t{ii} - mt;
    tt{ii} = tt{ii} - mt;
    sidx{ii} = 1:length(t{ii});
    if ii>1
        sidx{ii} = sidx{ii} + sidx{ii-1}(end);
    end
    
    %%%%%%%%%%%%%%%% use this section to compute for each session separately %%%%%%%%%%%%%%%%%%
    %compute coherence and plv
    Wxy = Wx{ii}.*conj(Wy{ii});
    Wxx = Wx{ii}.*conj(Wx{ii});
    Wyy = Wy{ii}.*conj(Wy{ii});
    keep = phase{ii}>15&phase{ii}<185;
    CS = Wxy./sqrt(mean(Wxx(:,keep),2).*mean(Wyy(:,keep),2));
    PLV = Wxy./abs(Wxy);
    clear Wxy Wxx Wyy
    
    %rotate CCS and PLV so that Re(CCS) and Re(PLV) is in direction of coherence
    for f = 1:length(freqVec)
        cdir = angle(mean(CS(f,keep)));
        Q = [cos(-cdir),-sin(-cdir);sin(-cdir),cos(-cdir)];
        rot = Q*[real(CS(f,:));imag(CS(f,:))];
        CS(f,:) = rot(1,:);
        
        cdir = angle(mean(PLV(f,keep)));
        Q = [cos(-cdir),-sin(-cdir);sin(-cdir),cos(-cdir)];
        rot = Q*[real(PLV(f,:));imag(PLV(f,:))];
        PLV(f,:) = rot(1,:);
    end %frequency
    
    %makes plots
    plot(freqVec,mean(CS(:,keep),2),'-r','LineWidth',2);hold on
    plot(freqVec,mean(PLV(:,keep),2),'-b','LineWidth',2);hold off
    legend('coherence','plv');
    ylim([0 1]);xlim([min(freqVec) max(freqVec)])
    xlabel('Frequency (Hz)');ylabel('coherence'); axis square
    bmpImage = sprintf('%s%s%s%s',sessions{ii},strcat('CS','\coherence'),'.bmp');
    saveas(gcf,bmpImage,'bmp');
    
    %save csc data as .csv
    %CS
    names = cell(1,length(freqVec)+1);
    for f = 1:length(freqVec)
        names{f} = strcat('CS_f',num2str(f));
    end % frequencies
    names{length(freqVec)+1} = 'Time';
    
    rtable = array2table([CS',tt{ii}'],'VariableNames',names);
    rfilename = sprintf('%s%s%s%s',strcat(sessions{ii},'CS\'),'rtable_CSC_CS.csv');
    writetable(rtable,rfilename);clear rtable    
    %PLV
    names = cell(1,length(freqVec)+1);
    for f = 1:length(freqVec)
        names{f} = strcat('PLV_f',num2str(f));
    end % frequencies
    names{length(freqVec)+1} = 'Time';
    
    rtable = array2table([PLV',tt{ii}'],'VariableNames',names);
    rfilename = sprintf('%s%s%s%s',strcat(sessions{ii},'CS\'),'rtable_CSC_PLV.csv');
    writetable(rtable,rfilename);clear rtable
    
    %save tracking data as .csv
    trnames = cell(1,2);
    trnames{1} = 'Time';
    trnames{2} = 'Position';
    rtable = array2table([t{ii},phase{ii}],'VariableNames',trnames);
    rfilename = sprintf('%s%s%s%s',strcat(sessions{ii},'CS\'),'rtable_tracking.csv');
    writetable(rtable,rfilename);clear rtable
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end %sessions

%%%%%%%%%%%%%%%% use this section to compute for all sessions simultaneously %%%%%%%%%%%%%%%%%%
% %concatenate variables across sessions
% tt = cell2mat(tt);
% t = cell2mat(t);
% Wx = cell2mat(Wx);
% Wy = cell2mat(Wy);
% phase = cell2mat(phase);
% 
% %compute coherence and plv
% Wxy = Wx.*conj(Wy);
% Wxx = Wx.*conj(Wx);
% Wyy = Wy.*conj(Wy);
% CS = Wxy./sqrt(mean(Wxx(:,phase>15&phase<185),2).*mean(Wyy(:,phase>15&phase<185),2));
% PLV = Wxy./abs(Wxy);
% clear Wx Wy Wxy Wxx Wyy
% 
% %rotate CCS and PLV so that Re(CCS) and Re(PLV) is in direction of coherence
% for f = 1:length(freqVec)
%     cdir = angle(mean(CS(f,phase>15&phase<185)));
%     Q = [cos(-cdir),-sin(-cdir);sin(-cdir),cos(-cdir)];
%     rot = Q*[real(CS(f,:));imag(CS(f,:))];
%     CS(f,:) = rot(1,:);
%     
%     cdir = angle(mean(PLV(f,phase>15&phase<185)));
%     Q = [cos(-cdir),-sin(-cdir);sin(-cdir),cos(-cdir)];
%     rot = Q*[real(PLV(f,:));imag(PLV(f,:))];
%     PLV(f,:) = rot(1,:);
% end %frequency
% 
% %now split data by session and save in each session's CS folder
% for ii = 1:numsessions
%     
%     %makes plots
%     bool = zeros(1,size(CS,2)); bool(sidx{ii}) = 1;
%     bool(phase<=15|phase>=185) = 0;bool = logical(bool);
%     plot(freqVec,mean(CS(:,bool),2),'-r','LineWidth',2);hold on
%     plot(freqVec,mean(PLV(:,bool),2),'-b','LineWidth',2);hold off
%     legend('coherence','plv');
%     ylim([0 1]);
%     xlabel('Frequency (Hz)');ylabel('coherence'); axis square
%     bmpImage = sprintf('%s%s%s%s',sessions{ii},strcat('CS','\coherence'),'.bmp');
%     saveas(gcf,bmpImage,'bmp');
%     
%     %save csc data as .csv
%     %CS
%     names = cell(1,length(freqVec)+1);
%     for f = 1:length(freqVec)
%         names{f} = strcat('CS_f',num2str(f));
%     end % frequencies
%     names{length(freqVec)+1} = 'Time';
%     
%     rtable = array2table([CS(:,sidx{ii})',tt(sidx{ii})'],'VariableNames',names);
%     rfilename = sprintf('%s%s%s%s',strcat(sessions{ii},'CS\'),'rtable_CSC_CS.csv');
%     writetable(rtable,rfilename);clear rtable    
%     %PLV
%     names = cell(1,length(freqVec)+1);
%     for f = 1:length(freqVec)
%         names{f} = strcat('PLV_f',num2str(f));
%     end % frequencies
%     names{length(freqVec)+1} = 'Time';
%     
%     rtable = array2table([PLV(:,sidx{ii})',tt(sidx{ii})'],'VariableNames',names);
%     rfilename = sprintf('%s%s%s%s',strcat(sessions{ii},'CS\'),'rtable_CSC_PLV.csv');
%     writetable(rtable,rfilename);clear rtable
%     
%     %save tracking data as .csv
%     trnames = cell(1,2);
%     trnames{1} = 'Time';
%     trnames{2} = 'Position';
%     rtable = array2table([t(sidx{ii}),phase(sidx{ii})],'VariableNames',trnames);
%     rfilename = sprintf('%s%s%s%s',strcat(sessions{ii},'CS\'),'rtable_tracking.csv');
%     writetable(rtable,rfilename);clear rtable
%         
% end
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%% Other Functions %%%%%%%%%%%%%%

function [dTargets,trackingColour] = decodeTargets(targets)

% Number of samples
numSamp = size(targets,2);

% Allocate memory to the array. 9 fields per sample: X-coord, Y-coord and
% 7 colour flag.
% Colour flag: 3=luminance, 4=rawRed, 5=rawGreen, 6=rawBlue, 7=pureRed,
% 8=pureGreen, 9=pureBlue.
dTargets = int16(zeros(numSamp,50,9));

for ii = 1:numSamp
    for jj = 1:50
        bitField = bitget(targets(jj,ii),1:32);
        if bitField(13)% Raw blue
            % Set the x-coord to the target
            ind = find(bitField(1:12));
            dTargets(ii,jj,1) = sum(2.^(ind-1));
            % Set the y-coord to the target
            ind = find(bitField(17:28));
            dTargets(ii,jj,2) = sum(2.^(ind-1));
            dTargets(ii,jj,6) = 1;
        end
        if bitField(14) % Raw green
            ind = find(bitField(1:12));
            dTargets(ii,jj,1) = sum(2.^(ind-1));
            ind = find(bitField(17:28));
            dTargets(ii,jj,2) = sum(2.^(ind-1));
            dTargets(ii,jj,5) = 1;
        end
        if bitField(15) % Raw red
            ind = find(bitField(1:12));
            dTargets(ii,jj,1) = sum(2.^(ind-1));
            ind = find(bitField(17:28));
            dTargets(ii,jj,2) = sum(2.^(ind-1));
            dTargets(ii,jj,4) = 1;
        end
        if bitField(16) % Luminance
            ind = find(bitField(1:12));
            dTargets(ii,jj,1) = sum(2.^(ind-1));
            ind = find(bitField(17:28));
            dTargets(ii,jj,2) = sum(2.^(ind-1));
            dTargets(ii,jj,3) = 1;
        end
        if bitField(29) % Pure blue
            ind = find(bitField(1:12));
            dTargets(ii,jj,1) = sum(2.^(ind-1));
            ind = find(bitField(17:28));
            dTargets(ii,jj,2) = sum(2.^(ind-1));
            dTargets(ii,jj,9) = 1;
        end
        if bitField(30) % Puregreen
            ind = find(bitField(1:12));
            dTargets(ii,jj,1) = sum(2.^(ind-1));
            ind = find(bitField(17:28));
            dTargets(ii,jj,2) = sum(2.^(ind-1));
            dTargets(ii,jj,8) = 1;
        end
        if bitField(31) % Pure red
            ind = find(bitField(1:12));
            dTargets(ii,jj,1) = sum(2.^(ind-1));
            ind = find(bitField(17:28));
            dTargets(ii,jj,2) = sum(2.^(ind-1));
            dTargets(ii,jj,7) = 1;
        end
    end
end

% Find out what colours were used in the tracking
trackingColour = zeros(1,7);
if ~isempty(find(dTargets(:,:,3),1)) % Luminance
    trackingColour(1) = 1;
end
if ~isempty(find(dTargets(:,:,7),1)) % Pure Red
    trackingColour(2) = 1;
end
if ~isempty(find(dTargets(:,:,8),1)) % Pure Green
    trackingColour(3) = 1;
end
if ~isempty(find(dTargets(:,:,9),1)) % Pure Blue
    trackingColour(4) = 1;
end
if ~isempty(find(dTargets(:,:,4),1)) % Raw Red
    trackingColour(5) = 1;
end
if ~isempty(find(dTargets(:,:,5),1)) % Raw Green
    trackingColour(6) = 1;
end
if ~isempty(find(dTargets(:,:,6),1)) % Raw Blue
    trackingColour(7) = 1;
end

% Exctracts the individual coordinates for the centre of mass of each
% tracking diode. The red LEDs are assumed to be at the front and the green
% diodes are assumed to be at the back.
function [frontX,frontY,backX,backY] = extractPosition(targets,tracking)

ind = find(tracking(2:end));
if length(ind) <= 1
    % Need at least two colours to get head direction
    disp('ERROR: To few LED colours have been tracked. Not possible to find head direction')
    frontX = NaN;
    frontY = NaN;
    backX = NaN;
    backY = NaN;
    return
else
    if ~tracking(2) && ~tracking(5)
        disp('ERROR: Red LED has not been tracked')
        frontX = NaN;
        frontY = NaN;
        backX = NaN;
        backY = NaN;
        return
    end
    if ~tracking(3) && ~tracking(6)
        disp('ERROR: Green LED has not been tracked')
        frontX = NaN;
        frontY = NaN;
        backX = NaN;
        backY = NaN;
        return
    end
end

% Number of samples in the data
numSamp = size(targets,1);

% Allocate memory for the arrays
frontX = zeros(1,numSamp);
frontY = zeros(1,numSamp);
backX = zeros(1,numSamp);
backY = zeros(1,numSamp);

% Exctract the front coordinates (red LED)
if tracking(2) && ~tracking(5)
    % Pure red but not raw red
    for ii = 1:numSamp
        ind = find(targets(ii,:,7));
        if ~isempty(ind)
            frontX(ii) = mean(targets(ii,ind,1));
            frontY(ii) = mean(targets(ii,ind,2));
        end
    end
end
if ~tracking(2) && tracking(5)
    % Not pure red but raw red
    for ii = 1:numSamp
        ind = find(targets(ii,:,4));
        if ~isempty(ind)
            frontX(ii) = mean(targets(ii,ind,1));
            frontY(ii) = mean(targets(ii,ind,2));
        end
    end
end
if tracking(2) && tracking(5)
    % Both pure red and raw red
    for ii = 1:numSamp
        ind = find(targets(ii,:,7) | targets(ii,:,4));
        if ~isempty(ind)
            frontX(ii) = mean(targets(ii,ind,1));
            frontY(ii) = mean(targets(ii,ind,2));
        end
    end
end

% Exctract the back coordinates (green LED)
if tracking(3) && ~tracking(6)
    % Pure green but not raw green
    for ii = 1:numSamp
        ind = find(targets(ii,:,8));
        if ~isempty(ind)
            backX(ii) = mean(targets(ii,ind,1));
            backY(ii) = mean(targets(ii,ind,2));
        end
    end
end
if ~tracking(3) && tracking(6)
    % Not pure green but raw green
    for ii = 1:numSamp
        ind = find(targets(ii,:,5));
        if ~isempty(ind)
            backX(ii) = mean(targets(ii,ind,1));
            backY(ii) = mean(targets(ii,ind,2));
        end
    end
end
if tracking(3) && tracking(6)
    % Both pure green and raw green
    for ii = 1:numSamp
        ind = find(targets(ii,:,8) | targets(ii,:,5));
        if ~isempty(ind)
            backX(ii) = mean(targets(ii,ind,1));
            backY(ii) = mean(targets(ii,ind,2));
        end
    end
end

% Estimates lacking position samples using linear interpolation. When more
% than timeTreshold sec of data is missing in a row the data is left as
% missing.
%
% Raymond Skjerpeng 2006.
function [theta,r] = interporPos(x,y,timeTreshold,sampRate)

theta = unwrap(atan2(y,x));
r = sqrt(x.^2+y.^2);

% Turn off warning
% warning('off','MATLAB:divideByZero');

% Number of samples that corresponds to the time threshold.
sampTreshold = floor(timeTreshold * sampRate);

% number of samples
numSamp = length(x);
% Find the indexes to the missing samples
% temp1 = 1./x;
% indt1 = isinf(temp1);
ind = isnan(x);
ind2 = find(ind==1);
% Number of missing samples
N = length(ind2);

if N == 0
    % No samples missing, and we return
    return
end

change = 0;

% Remove NaN in the start of the path
if ind2(1) == 1
    change = 1;
    count = 0;
    while 1
        count = count + 1;
        if ind(count)==0
            break
        end
    end
    theta(1:count) = theta(count);
    r(1:count) = r(count);
end

% Remove NaN in the end of the path
if ind2(end) == numSamp
    change = 1;
    count = length(x);
    while 1
        count = count - 1;
        if ind(count)==0
            break
        end
    end
    theta(count:numSamp) = theta(count);
    r(count:numSamp) = r(count);
end

if change
    % Recalculate the missing samples
    %     temp1 = 1./r;
    %     indt1 = isinf(temp1);
    ind = isnan(r);
    % Missing samples are where both x and y are equal to zero
    ind2 = find(ind==1);
    % Number of samples missing
    N = length(ind2);
end

for ii = 1:N
    % Start of missing segment (may consist of only one sample)
    start = ind2(ii);
    % Find the number of samples missing in a row
    count = 0;
    while 1
        count = count + 1;
        if ind(start+count)==0
            break
        end
    end
    % Index to the next good sample
    stop = start+count;
    if start == stop
        % Only one sample missing. Setting it to the last known good
        % sample
        theta(start) = theta(start-1);
        r(start) = r(start-1);
    else
        if count < sampTreshold
            % Last good position before lack of tracking
            theta1 = theta(start-1);
            r1 = r(start-1);
            % Next good position after lack of tracking
            theta2 = theta(stop);
            r2 = r(stop);
            % Calculate the interpolated positions
            theta(start:stop) = interp1([1,2],[theta1,theta2],1:1/count:2);
            r(start:stop) = interp1([1,2],[r1,r2],1:1/count:2);
            % Increment the counter (avoid estimating allready estimated
            % samples)
            ii = ii+count;
        else
            % To many samples missing in a row and they are left as missing
            ii = ii+count;
        end
    end
end
theta = wraptopi(theta);

% Calculates the direction of the head stage from the two set of
% coordinates. If one or both coordinate sets are missing for one samle the
% direction is set to NaN for that sample. Direction is also set to NaN for
% samples where the two coordinate set are identical. Returns the
% direction in degrees
function direct = headDirection(frontX,frontY,backX,backY)

% Number of position samples in data set
N = length(frontX);
direct = zeros(N,1);

for ii = 1:N
    
    if frontX(ii)==0 || backX(ii)==0 || isnan(frontX(ii)) || isnan(backX(ii))
        % One or both coordinates are missing. No angle.
        direct(ii) = NaN;
        continue
    end
    
    % Calculate the difference between the coordinates
    xd = frontX(ii) - backX(ii);
    yd = frontY(ii) - backY(ii);
    
    if xd==0
        if yd==0
            % The two coordinates are at the same place and it is not
            % possible to calculate the angle
            direct(ii) = NaN;
            continue
        elseif yd>0
            direct(ii) = 90;
            continue
        else
            direct(ii) = 270;
            continue
        end
    end
    if yd==0
        if xd>0
            % Angle is zero
            continue
        else
            direct(ii) = 180;
            continue
        end
    end
    
    if frontX(ii)>backX(ii) && frontY(ii)>backY(ii)
        % Angle between 0 and 90 degrees
        direct(ii) = atan(yd/xd) * 360/(2*pi);
        
    elseif frontX(ii)<backX(ii) && frontY(ii)>backY(ii)
        % Angle between 90 and 180 degrees
        direct(ii) = 180 - atan(yd/abs(xd)) * 360/(2*pi);
        
    elseif frontX(ii)<backX(ii) && frontY(ii)<backY(ii)
        % Angle between 180 and 270 degrees
        direct(ii) = 180 + atan(abs(yd)/abs(xd)) * 360/(2*pi);
        
    else
        % Angle between 270 and 360 degrees
        direct(ii) = 360 - atan(abs(yd)/xd) * 360/(2*pi);
    end
end

function [sptimes] = spikeTimes(ts,loct)

N = length(ts);
sptimes = zeros(N,4);

for ii = 1:N
    tdiff = loct-ts(ii);
    sptimes(ii,1) = find(tdiff<=0,1,'last');
    sptimes(ii,2) = tdiff(sptimes(ii,1));
    sptimes(ii,3) = find(tdiff>0,1,'first');
    sptimes(ii,4) = tdiff(sptimes(ii,3));
end

