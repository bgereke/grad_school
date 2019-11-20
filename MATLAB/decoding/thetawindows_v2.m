function [windows,bp,eegTS,eeg,bpd,thetadelta] = thetawindows_v3
%version to specifically list which csc to use below
%also to find minima and then build window around that

% tet = gettetnumbers(TTfile);
% disp('Using this Tetrode for theta windows:');
% disp(strcat('CSC',int2str(tet),'.ncs'));
% eegfile = strcat('CSC',int2str(tet),'.ncs');

eegfile = 'CSC5.ncs';

%ratio cutoff
ratiocut = 3;
%smooth ratio by
smo = 1000; %in index points 
%theta
f1 = 6;
f2 = 12;
%delta
f1d = 2;
f2d = 4;

minwindow = 1/f2;%in seconds
maxwindow = 1/f1;

[eeg,eegTS] = loadEeg8(eegfile);
TFR = (TFR_frequency_band(eeg,2000,5,f1,f2)); %theta
TFRd = (TFR_frequency_band(eeg,2000,5,f1d,f2d));%delta
thetadelta = smooth(TFR./TFRd,smo);

%zTFRr = (TFR_frequency_band(eeg,2000,5,150,250));%ripple
bp = fftbandpass(eeg,2000,f1-2,f1,f2,f2+2);%theta
bpd = fftbandpass(eeg,2000,f1d-1,f1d,f2d,f2d+1);%delta
[xmx,imx,~,~] =extrema(bp);

imx = sort(imx);
times = nan(length(xmx),2);

%get time windows between each max
for i = 1:length(xmx)-1
times(i,1) = eegTS(imx(i));
times(i,2) = eegTS(imx(i+1));    
end

%time window lengths
times(:,3) = times(:,2) - times(:,1); 

%copy over the windows that have power above the ratio cutoff and are
%correct length
windows = [];
for i = 1:length(times)
    if mean(thetadelta(eegTS>=times(i,1) & eegTS<=times(i,2))) >= ratiocut
        if times(i,3) >= minwindow && times(i,3) <= maxwindow
            windows = [windows;  times(i,1:2)];
        end
    end
end
[~,~,t] = readVideoData('VT1.nvt',.275);
ind = find(windows(:,2)>max(t) | windows(:,1)<min(t));
windows (ind,:) = [];

function tets = gettetnumbers(file)

tets = nan(numel(textread(file,'%1c%*[^\n]')),1); %#ok<REMFF1>

fid = fopen(file,'r');
if fid == -1
    msgbox('Could not open the input file','ERROR');
end

for i = 1:size(tets,1)
  tline = fgetl(fid);
  if ~ischar(tline) 
      break 
  else
      if tline(4)=='_'
        tets(i) = str2double(tline(3));
      elseif tline(5) == '_'
        tets(i) = nan;
      end
  end          
end
       
fclose(fid);

tets(isnan(tets)) = [];
tets(tets>6) = [];
tets = mode(tets);
%tets = unique(tets);
%tets = num2cell(tets);

