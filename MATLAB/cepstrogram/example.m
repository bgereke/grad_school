clear, clc, close all

% load a sound file
[x, fs] = audioread('sample1.wav');   % read the file
x = x(:, 1);                        % get the first channel
x = x/max(abs(x));                  % normalize

% analysis parameters
wlen = 100;                        % window length
h = 1;                            % window hop size

% calculate the cepstrogram
[C, q, t] = cepstrogram(x, wlen, h, fs);

% some conditioning
q = q*1000;                         % convert s to ms
q = q(q > 0.25);                    % ignore all quefrencies bellow 0.25 ms                    
C = C(end-length(q)+1:end, :);      % ignore all cepstrum coefficient for 
                                    % quefrencies bellow 0.25 ms  

% plot the signal spectrogram
figure(1)
subplot(4, 1, 1) 
plot((0:length(x)-1)/fs, x, 'r')
grid on
axis([0 max(t) -1.1 1.1])
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
xlabel('Time, s')
ylabel('Amplitude')
title('The signal in the Time domain')

subplot(4, 1, 2:4) 
[T, Q] = meshgrid(t, q);
surf(T, Q, C)
shading interp
box on
axis([0 max(t) 0 max(q)])
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
xlabel('Time, s')
ylabel('Quefrency, ms')
title('Cepstrogram of the signal')
view(0, 90)

% represent all cepstrum coefficients bellow 0 with 
% one and the same color for better visualization
colormap jet(8)
caxis([0 max(C(:))])