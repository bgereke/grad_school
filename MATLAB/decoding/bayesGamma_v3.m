function [scores] =  bayesGamma_v3(scores,listfile)

%version gets trace, then prewhitens, then TFR
%also changed prewhitening calue to three and wideded the trace window by
%300/2000 seconds on either side

%modified to save some extra stuff: tetrode numbers used, pwd, tfr matrix
%for each window and traces

%also, a prewhitening filter was incorporated into the TFR

%f1 = 60;
f2 = 100;
s1 = 25;
%s2 = 55;

tetnums = gettetnumbers(listfile);

%load the tetrodes of interest and get their TFRs
disp('Calculating power for the following Tetrodes:');
for i = 1:size(tetnums,1)
    disp(strcat('CSC',int2str(tetnums{i,1}),'.ncs'));
    temptet = strcat('CSC',int2str(tetnums{i,1}),'.ncs');
    [tetnums{i,2},~,tetnums{i,3}] = loadEeg8(temptet);       
end

%go through each window listed and get the power for fast and slow
for i = 1:size(scores,1)
    if ~isnan(scores{i,22})
        startEegInd = spikeTStoEegIndex(scores{i,2}(1,1),tetnums{1,3},2000)-100;
        stopEegInd =  spikeTStoEegIndex(scores{i,2}(1,2),tetnums{1,3},2000)+100;
        if startEegInd < 1
            startEegInd = 1;
        end
        if stopEegInd > length(tetnums{1,2})
            stopEegInd = length(tetnums{1,2});
        end

        rawtraces = [];
        tfrWhole = [];

        for j = 1:size(tetnums,1)            
            rawtraces(:,j) = tetnums{j,2}(startEegInd:stopEegInd);
            tfrWhole(:,:,j) = TFR_frequency_band_bayes(rawtraces(:,j),2000,5,s1,f2);
        end        
        scores{i,27} = mean(tfrWhole,3);
        scores{i,28} = rawtraces;
        scores{i,29} = tetnums(:,1);
    end
end

 %get mean slow/fast power for the window from each tet
%     meanF = 0;meanFz = 0;
%     meanS = 0;meanSz = 0;
%     tfrWhole = [];

%tetnums{i,8} = TFR_frequency_band_bayes(tetnums{i,2},2000,5,s1,f2);
    %tetnums{i,4} = mean(tetnums{i,8}(19:end,:),1);  
    %tetnums{i,5} = mean(tetnums{i,8}(1:16,:),1); %gives the mean across frequencies, not time
    %tetnums{i,6} = zscore(tetnums{i,4});
    %tetnums{i,7} = zscore(tetnums{i,5}); 
    
    % meanF = meanF + mean(tetnums{j,4}(startEegInd:stopEegInd));
%             meanS = meanS + mean(tetnums{j,5}(startEegInd:stopEegInd));
%             meanFz = meanFz + mean(tetnums{j,6}(startEegInd:stopEegInd));
%             meanSz = meanSz + mean(tetnums{j,7}(startEegInd:stopEegInd));
%             tfrWhole(:,:,j) = tetnums{j,8}(:,startEegInd:stopEegInd); 

% 
% scores{i,23} = meanF/j;
%         scores{i,24} = meanS/j;
%         scores{i,25} = meanFz/j;
%         scores{i,26} = meanSz/j;
%         scores{i,27} = tfrWhole;




%figure out which tetrodes to use
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
tets = unique(tets);
tets = num2cell(tets);

function [TFR] = TFR_frequency_band_bayes(W,Fs,width, f1, f2)

W = prewhitening(W,3);

freq(1,:) = [f1 f2];

n = 0;
for l = freq(1,1):2:freq(1,2)
    n = n + 1;
    temp(n,:) = zeros(1,size(W,1));
    temp(n,:) = energyvec(l,detrend(W),Fs,width);
end
%TFR(1,1:size(W,1)) = mean(temp, 1);
TFR = (temp); %return whole thing


function y = energyvec(f,s,Fs,width)
% function y = energyvec(f,s,Fs,width)
%
% Return a vector containing the energy as a
% function of time for frequency f. The energy
% is calculated using Morlet's wavelets. 
% s : signal
% Fs: sampling frequency
% width : width of Morlet wavelet (>= 5 suggested).

dt = 1/Fs;
sf = f/width;
st = 1/(2*pi*sf);

t=-3.5*st:dt:3.5*st;
m = morlet(f,t,width);
y = convfft(s,m);
%y = conv(s,m);
y = (2*abs(y)/Fs).^2;
y = y(ceil(length(m)/2):length(y)-floor(length(m)/2));

function y = morlet(f,t,width)
% function y = morlet(f,t,width)
% 
% Morlet's wavelet for frequency f and time t. 
% The wavelet will be normalized so the total energy is 1.
% width defines the ``width'' of the wavelet. 
% A value >= 5 is suggested.
%
% Ref: Tallon-Baudry et al., J. Neurosci. 15, 722-734 (1997)
%
% See also: PHASEGRAM, PHASEVEC, WAVEGRAM, ENERGY 
%
% Ole Jensen, August 1998 

sf = f/width;
st = 1/(2*pi*sf);
A = (st*sqrt(pi))^(-0.5);

y = A*exp(-t.^2/(2*st^2)).*exp(i*2*pi*f.*t);