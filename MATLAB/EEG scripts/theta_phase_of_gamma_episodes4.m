function [phaseBin_theta_slow_gamma_norm, phaseBin_theta_fast_gamma_norm] = theta_phase_of_gamma_episodes4(EEG, non_overlap_peak_ind_unique_slow, non_overlap_peak_ind_unique_fast, Fs)

EEG_theta_filter = bandpass(EEG, 6, 10, Fs);
DTAS = hilbert(EEG_theta_filter);
theta_phase = angle(DTAS);


theta_phase_non_overlap_peak_ind_unique_slow = theta_phase(non_overlap_peak_ind_unique_slow);
theta_phase_non_overlap_peak_ind_unique_fast = theta_phase(non_overlap_peak_ind_unique_fast);




[phaseBin_theta_slow_gamma] = make_spike_time_histogram4(theta_phase_non_overlap_peak_ind_unique_slow);
phaseBin_theta_slow_gamma_norm = phaseBin_theta_slow_gamma/length(non_overlap_peak_ind_unique_slow);

[phaseBin_theta_fast_gamma] = make_spike_time_histogram4(theta_phase_non_overlap_peak_ind_unique_fast);
phaseBin_theta_fast_gamma_norm = phaseBin_theta_fast_gamma/length(non_overlap_peak_ind_unique_fast);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[dataFiltered] = bandpass( data, lowFreq, highFreq, Fs )
%
% Filter data with FFT based bandpass. Bandpass frequencies are defined by
% 'lowFreq' and 'highFreq' (in Hz), and assuming the sampling frequency defined
% by Fs (in kHz).  Note from LLC on Jan 9 2006- I changed this to read the
% sampling frequency in Hz.
%
% data: each column is assumed to be a continuous trace.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [dataFiltered] = bandpass( data, lowFreq, highFreq, Fs )


%--- reshape data to be vertical ---

if( size(data,1) < size(data,2) )
	data = data';
end 


%--- Filter parmeters ---

%Fs = Fs * 1000;		% convert to Hz
lengthTrace = size( data, 1 );
nCutLow  = (lowFreq * lengthTrace) / Fs + 1;
nCutHigh = (highFreq * lengthTrace) / Fs + 1;


%===== Filter the data =====

F = fft( data ); 

%fig(1)
%plot(abs(F))

cx0 = complex(0);

halfLength = floor( lengthTrace / 2 );


%--- low cut ---

if (1 < nCutLow)					% the zero component ...
   F( 1, : ) = cx0;          	    
end

for k = 2:nCutLow-1		
    F( k, : ) = cx0;    	      	    % Set this component to zero...
    F( lengthTrace - k + 2, : ) = cx0;	% ...for all channels
end


%--- high cut ---

for k = nCutHigh+1:halfLength+1			% the + 1 takes care of the even/odd cases
   F( round(k), : ) = cx0;                	    % Set this component to zero...
   F( lengthTrace - round(k) + 2, : ) = cx0;		% ...for all channels
end


%hold on, plot(abs(F),'r'), hold off


dataFiltered = real( ifft(F) );




function [phaseBin] = make_spike_time_histogram4(sPhase);

% Convert sPhase to degrees
sPhase = (sPhase+pi) * 360/(2*pi);

start = 0;
stop = 30;
phaseBin = zeros(12,1);
for ii=1:12
ind = find(sPhase >= start & sPhase < stop);
phaseBin(ii) = length(ind);
start = start + 30;
stop = stop + 30;
end
