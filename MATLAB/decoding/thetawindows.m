function [data] = thetawindows(data)

%for choosing eeg w/most cells
% tet = gettetnumbers(TTfile);%returns the tetrode number with msot cells
% disp('Using this Tetrode for theta windows:');
% disp(strcat('CSC',int2str(tet),'.ncs'));
% eegfile = strcat('CSC',int2str(tet),'.ncs');
% [data.eeg,eegTS] = loadEeg8(eegfile);

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

% [eeg,eegTS] = loadEeg8(eegfile);
% TFR = mean(data.TFR(data.freqVec>f1 & data.freqVec<f2,:),1); %theta
% TFRd = mean(data.TFR(data.freqVec>f1d & data.freqVec<f2d,:),1);%delta
% data.thetadelta = smooth(TFR./TFRd,smo);

%zTFRr = (TFR_frequency_band(eeg,2000,5,150,250));%ripple
%bp = fftbandpass(data.eeg,2000,f1-2,f1,f2,f2+2);%theta
DTAS = hilbert(data.bp);
data.theta_phase = angle(DTAS);
%bpd = fftbandpass(data.eeg,2000,f1d-1,f1d,f2d,f2d+1);%delta
[~,~,xmx,imx] =extrema(data.bp); %left two for maxima

imx = sort(imx);
times = nan(length(xmx),2);
wimx = nan(length(xmx),2);

%get time windows between each max
for i = 1:length(xmx)-1
    times(i,1) = data.tsi(imx(i));
    times(i,2) = data.tsi(imx(i+1));
    wimx(i,1) = imx(i);
    wimx(i,2) = imx(i+1);
end

%time window lengths
times(:,3) = times(:,2) - times(:,1); 

%copy over the windows that have power above the ratio cutoff and are
%correct length
data.w = [];
data.widx = [];

for i = 1:length(times)
    if mean(data.thetadelta(data.tsi>=times(i,1) & data.tsi<=times(i,2))) >= ratiocut
        if times(i,3) >= minwindow && times(i,3) <= maxwindow
            data.w = [data.w;  times(i,1:2)];
            data.widx = [data.widx; wimx(i,1:2)];
        end
    end
end
[~,~,t] = readVideoData('VT1.nvt',.275);
data.w(data.w(:,2)>max(t) | data.w(:,1)<min(t),:) = [];
data.widx(data.w(:,2)>max(t) | data.w(:,1)<min(t),:) = [];


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

