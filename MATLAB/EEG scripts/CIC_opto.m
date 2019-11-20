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

function CIC_opto(inFile,numperms)  % width - number of cycles in wavelet (> 5 advisable); order - whitening

img_text = 'on';

fid = fopen(inFile,'r');
if fid == -1
    msgbox('Could not open the input file! Make sure the filname and path are correct.','ERROR');
end

% Get sessions and csc-file list from input file
% Get sessions and csc-file list from input file
fid = fopen(inFile,'r');
ii = -1;
numsessions = 0;
while ~feof(fid)
    str = fgetl(fid);
    if ii == 0
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
        sessions(numsessions) = {str};
    end
    ii = ii+1;
end

%make directory to save figures
savedir = cd;
dirInfo = dir(cd);
found = 0;
for kk=1:size(dirInfo,1)
    if dirInfo(kk).isdir
        if strcmp(dirInfo(kk).name,strcat('CIC_opto','\'))
            found = 1;
        end
    end
end
if found==0
    mkdir(savedir,strcat('CIC_opto','\'));
end

veloff = [];
velon = [];
ImXoff = [];
ImXon = [];

for jj=1:numsessions   
        
    if jj == 1 || jj == 2        
        
        %resize ImX to one sample per vel
        data = load([sessions{jj},strcat('CS','\CS'),'.mat']);
        t = data.time-min(data.time);
        ImX = imag(data.PCS).*repmat(sign(mean(imag(data.PCS),2)),1,size(data.PCS,2));
        vel = data.vel;
        freqVec = data.freqVec;
        
        efile = strcat(sessions{jj},'Events.nev');
        FS = [1 0 1 0 1]; EH = 0; EM = 1;
        [TimeStamps, TTLs, EventStrings] = Nlx2MatEV(efile,FS,EH,EM);
        TTLs(1) = []; TimeStamps(1) = []; TimeStamps = TimeStamps/1000000;
        if isempty(TTLs)
            msgbox('No TTLs in input file!','ERROR');
        end
        TS = TimeStamps - min(data.time);
        
        tidxoff = [];
        tidxon = [];
        for ts=2:length(TS)
            tidx = find(t>TS(ts-1)&t<=TS(ts));
            if TTLs(ts) == 0
                tidxoff = [tidxoff tidx'];
            else
                tidxon = [tidxon tidx'];
            end
        end

        velon = [velon vel(tidxon)'];
        veloff = [veloff vel(tidxoff)'];
        
        ImXon = [ImXon ImX(:,tidxon)];
        ImXoff = [ImXoff ImX(:,tidxoff)];
        
        clear ImX
    end
end

cic(:,1) = mean(ImXoff,2);
cic(:,2) = mean(ImXon,2);

for jj=3:numsessions
    
    %resize ImX to one sample per vel
    data = load([sessions{jj},strcat('CS','\CS'),'.mat']);
    t = data.time-min(data.time);
    ImX = imag(data.PCS).*repmat(sign(mean(imag(data.PCS),2)),1,size(data.PCS,2));
    vel = data.vel;
    freqVec = data.freqVec;
    clear data
    
end

firstlight = linspace(0,30,numperms);
p_oncic = zeros(length(freqVec),numperms);
p_offcic = zeros(length(freqVec),numperms);

for b = 1:numperms
    ImXon = []; ImXoff = [];
    velon = []; veloff = [];
    
    for jj = 3:numsessions
        tidxon = []; tidxoff = [];
        TS = firstlight(b):30:max(t);
        pTTLs = ones(size(TS));
        pTTLs(1:2:end) = 0;
        
        for ts=2:length(TS)
            tidx = find(t>TS(ts-1)&t<=TS(ts));
            if pTTLs(ts) == 0
                tidxoff = [tidxoff tidx'];
            else
                tidxon = [tidxon tidx'];
            end
        end
        
        velon = [velon vel(tidxon)'];
        veloff = [veloff vel(tidxoff)'];
        ImXon = [ImXon ImX(:,tidxon)];
        ImXoff = [ImXoff ImX(:,tidxoff)];
        
    end
    
    p_offcic(:,b) = mean(ImXoff,2);
    p_oncic(:,b) = mean(ImXon,2);
end

pstd = std([p_oncic-p_offcic p_offcic-p_oncic],0,2);

close all
figure(1);
subplot(1,2,1)
[ci] = prctile([(p_oncic-p_offcic)./p_offcic*100 (p_offcic-p_oncic)./p_oncic*100]',[2.5 97.5]);
shadedplot(freqVec, ci(2,:), ci(1,:),[1 0 0],[1 0.25 0.25]); alpha(0.5);hold on
plot(freqVec,mean([p_oncic-p_offcic p_offcic-p_oncic],2),'-r')
plot(freqVec,(cic(:,2)-cic(:,1))./cic(:,1)*100,'-k')
xlabel('frequency (Hz)')
ylabel('percent change from no light')

subplot(1,2,2)
[ci] = prctile([p_oncic-p_offcic p_offcic-p_oncic]'./repmat(pstd,1,2*numperms)',[2.5 97.5]);
shadedplot(freqVec, ci(2,:), ci(1,:),[1 0 0],[1 0.25 0.25]); alpha(0.5);hold on
plot(freqVec,mean([p_oncic-p_offcic p_offcic-p_oncic]'./repmat(pstd,1,2*numperms)'),'-r')
plot(freqVec,(cic(:,2)-cic(:,1))./pstd,'-k')
xlabel('frequency (Hz)')
ylabel('z-score from control')

figImage = sprintf('%s%s%s%s',savedir,strcat('\CIC_opto','\'),'CIC_opto.fig');
bmpImage = sprintf('%s%s%s%s',savedir,strcat('\CIC_opto','\'),'CIC_opto.bmp');
saveas(gcf,figImage,'fig');
saveas(gcf,bmpImage,'bmp');



%_______________________________________________________________________
%other functions
%_______________________________________________________________________

% Centre the path/box for the two sessions. Both boxes are set according to
% the path/coordinates to the box for the first session.
function [posx,posy] = centreBox(posx,posy)

% Find border values for box for session 1
maxX = max(posx);
minX = min(posx);
maxY = max(posy);
minY = min(posy);

% Set the corners of the reference box
NE = [maxX, maxY];
NW = [minX, maxY];
SW = [minX, minY];
SE = [maxX, minY];

% Get the centre coordinates of the box
centre = findCentre(NE,NW,SW,SE);
% Centre both boxes according to the coordinates to the first box
posx = posx - centre(1);
posy = posy - centre(2);


% Calculates the centre of the box from the corner coordinates
function centre = findCentre(NE,NW,SW,SE);

% The centre will be at the point of interception by the corner diagonals
a = (NE(2)-SW(2))/(NE(1)-SW(1)); % Slope for the NE-SW diagonal
b = (SE(2)-NW(2))/(SE(1)-NW(1)); % Slope for the SE-NW diagonal
c = SW(2);
d = NW(2);
x = (d-c+a*SW(1)-b*NW(1))/(a-b); % X-coord of centre
y = a*(x-SW(1))+c; % Y-coord of centre
centre = [x,y];




