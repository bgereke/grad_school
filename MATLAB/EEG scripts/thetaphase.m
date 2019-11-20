function [theta_phase,wave_phase,sym,slowBP] = thetaphase(S,Fs)

Ftheta = 9;    % Hz
slowBP = fftbandpass(S,Fs,Ftheta-5,Ftheta-3,Ftheta+3,Ftheta+5);
fastBP = fftbandpass(S,Fs,Ftheta-5,Ftheta-4,80,82);
DTAS = hilbert(slowBP);
theta_phase = angle(DTAS);
[henv,lenv] = myenvelope(fastBP,1,'peak');
waveBP = (henv+lenv)/2;
% h = hilbert(waveBP);
% wave_phase = angle(h);
sym = [];

% %Belluscio waveform-based phase estimate
% tridx = find(diff(theta_phase)<-1.5*pi)+1;
% pkidx = find(diff(sign(theta_phase))==2)+1;
% % [~,tridx] = findpeaks(-slowBP,'MinPeakProminence',25,'MinPeakDistance',0.08*2000);
% % [~,pkidx] = findpeaks(slowBP,'MinPeakProminence',25,'MinPeakDistance',0.08*2000);
% for i = 1:length(pkidx)-1
%    if isempty(tridx(tridx>pkidx(i)&tridx<pkidx(i+1))) 
%        [~,newmin] = min(slowBP(pkidx(i):pkidx(i+1)));
%        tridx = [tridx newmin+pkidx(i)-1];
%    end
% end
% tridx = sort(tridx);
% % pkidx(abs(theta_phase(pkidx))>0.05) = [];
% maxima = zeros(size(tridx)-1);
% minima = zeros(size(tridx));
% for i = 1:length(tridx)-1
%     wigglewin = tridx(i):tridx(i+1);
%     lw = round(length(wigglewin)/4);
%     wigglewin = wigglewin(lw:end-lw);
%     [~,maxima(i)] = max(waveBP(tridx(i):tridx(i+1)));
%     maxima(i) = maxima(i) + tridx(i) - 1;
% end
% maxima = [1 maxima];
% for i = 1:length(minima)-1
%     [~,minima(i)] = min(waveBP(maxima(i):maxima(i+1)));
%     minima(i) = minima(i) + maxima(i) - 1;
% end
% [~,minima(end)] = min(waveBP(maxima(end):end)); 
% minima(end) = minima(end) + maxima(end) - 1;
% if minima(1) ~= 1
%     if ~isempty(findpeaks(slowBP(1:minima(1))))
%         [~,maxima(1)] = max(waveBP(1:minima(1)));
%     else
%         maxima(1) = [];
%     end
% end
% if minima(end) ~= length(theta_phase)
%     if ~isempty(findpeaks(slowBP(minima(end):end)))
%         [~,midx] = max(waveBP(minima(end):end));
%         maxima = [maxima midx+minima(end)-1];
%     end
% end

%interpolate between peaks and troughs
wave_phase = zeros(size(theta_phase));
sym = [];
% if min(maxima)<min(minima)
%     for i = 1:length(minima)-1
%         wave_phase(minima(i):maxima(i+1)) = linspace(-pi,0,length(minima(i):maxima(i+1))); %trough-to-peak
%         wave_phase(maxima(i+1):minima(i+1)) = linspace(0,pi,length(maxima(i+1):minima(i+1))); %peak-to-trough
%         sym = [sym;round(mean(minima(i:i+1))), log(length(minima(i):maxima(i+1))/length(maxima(i+1):minima(i+1)))];
%     end
%     wave_phase(1:maxima(1)) = linspace(theta_phase(1),0,length(1:maxima(1)));
%     wave_phase(maxima(1):minima(1)) = linspace(0,pi,length(maxima(1):minima(1)));
% else
%     for i = 1:length(minima)-1
%         wave_phase(minima(i):maxima(i)) = linspace(-pi,0,length(minima(i):maxima(i))); %trough-to-peak
%         wave_phase(maxima(i):minima(i+1)) = linspace(0,pi,length(maxima(i):minima(i+1))); %peak-to-trough
%         sym = [sym;round(mean(minima(i:i+1))), log(length(minima(i):maxima(i))/length(maxima(i):minima(i+1)))];
%     end
%     wave_phase(1:minima(1)) = linspace(theta_phase(1),pi,length(1:minima(1)));
% end
% if max(maxima)>max(minima)
%     wave_phase(minima(end):maxima(end)) = linspace(-pi,0,length(minima(end):maxima(end)));
%     wave_phase(maxima(end):end) = linspace(0,theta_phase(end),length(maxima(end):length(wave_phase)));
% else
%     wave_phase(minima(end):end) = linspace(-pi,theta_phase(end),length(minima(end):length(wave_phase)));
% end