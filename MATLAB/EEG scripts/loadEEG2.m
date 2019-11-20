function [samples,ts,tt, Fs, bv, ir] = loadEEG2(file)

% Set the field selection for reading CSC files. 1 = Add parameter, 0 = skip
% parameter
fieldSelection(1) = 1; % Timestamps
fieldSelection(2) = 0; % Channel number
fieldSelection(3) = 1; % Sample Frequency
fieldSelection(4) = 0; % Number of valid samples
fieldSelection(5) = 1; % Samples (EEG data)
% Do we return header, 1 = Yes, 0 = No.
extractHeader = 1;
% 5 different extraction modes, see help file for Nlx2MatCSC_v3

extractMode = 1; % Extract all data

[ts,Fs,samp, header] =...
   Nlx2MatCSC(file,fieldSelection,extractHeader,extractMode);

% Transform the 2-D samples array to an 1-D array
M = size(samp,2);
samples = zeros(512*M,1);
%convert ts to seconds
%ts = ts/1000000;

%Convert the sampling rate to a scalar if there is only one rate
if(min(Fs)==max(Fs))
	Fs=Fs(1);
else
	printf('\nNon-constant sampling rate detected in %s\n',file);
end

for jj = 1:M
    samples(((jj-1)*512)+1:512*jj) = samp(1:512,jj);
end

bv = str2num(header{12}(13:end))*10^6;  %bit-microvolt conversion
ir = str2num(header{15}(13:end));  %Input range in microvolts

for aa=1:length(ts)-1
	fact = 512*(aa-1);
	temp = linspace(ts(aa),ts(aa+1),513);
	tt(1+fact:512+fact) = temp(1:end-1);
end

endbuf = max(tt)+1000000/Fs:1000000/Fs:max(tt)+1000000*512/Fs;
tt = [tt endbuf];
tt = tt/1000000;
ts = ts/1000000;
