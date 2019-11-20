%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Cepstrogram with MATLAB Implementation     %
%                                                %
% Author: Ph.D. Eng. Hristo Zhivomirov  08/25/16 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ceps, q, t] = cepstrogram(x, wlen, h, fs)

% function: [ceps, q, t] = cepstrogram(x, wlen, h, fs)
% x - signal in the time domain
% wlen - length of the hamming window
% h - hop size
% fs - sampling frequency, Hz
% ceps - cepstrogram matrix (only unique points, time across columns, quefrency across rows)
% q - quefrency vector, s
% t - time vector, s

% represent x as column-vector
x = x(:);

% length of the signal
xlen = length(x);

% form a periodic hamming window
win = hamming(wlen, 'periodic');

% form the stft matrix
rown = ceil((1+wlen)/2);            % calculate the total number of rows                      
coln = 1+fix((xlen-wlen)/h);        % calculate the total number of columns
ceps = zeros(rown, coln);           % form the cepstrogram matrix

% initialize the indexes
indx = 0;
col = 1;

% cepstrogram computing
while indx + wlen <= xlen
    % windowing
    xw = x(indx+1:indx+wlen).*win;
    
    % compute the cepstrum
    c = real(ifft(log(abs(fft(xw)))));
    
    % update the cepstrogram matrix
    ceps(:, col) = c(1:rown);
    
    % update the indexes
    indx = indx + h;
    col = col + 1;
end

% calculate the time and quefrency vectors
t = (wlen/2:h:wlen/2+(coln-1)*h)/fs;
q = (0:rown-1)/fs;

end