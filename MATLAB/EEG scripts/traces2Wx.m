function [Wx] = traces2Wx(x,freqVec,Fs,width,type,normalization)

%Inputs:
%x - n x 1 column vector (i.e. your signal)
%freqVec - vector of frequencies at which to evaluate (in Hz)
%Fs - a scalar specifying the sampling rate (in Hz)
%width - a scalar specifying the width of the Morlet wavelet (usually between 5-6, not <5!)
%        for STFT this fixed and taken to be Gaussian variance in seconds
%type - type of transform... either 'Morlet' or 'STFT'
%normalization - can be 'energy' or 'area' 
%
%Outputs:
%Wx - nfreq x nsamp matrix of complex-valued wavelet coefficients
%
%Notes:
%Amplitude = abs(Wx)
%Phase = angle(Wx)
%Power ~ abs(Wx).^2 

x = x';
Wx = zeros(length(freqVec),size(x,2)); 

%detrend based on lowest frequency
g = gaussian(freqVec(1),width,Fs);
n = length(x)+length(g)-1;
y = real(ifft(fft(x,n).*fft(g,n)));
y = y(ceil(length(g)/2):length(y)-floor(length(g)/2));
x = x-y;
    
for j=1:length(freqVec)
 
    %reflect boundaries to reduce edge effects
    npnts = round(Fs/freqVec(j)*width/2);
    xr = [x(npnts+1:-1:2) x x(end-1:-1:end-npnts)];
    
    if strcmp(type,'Morlet')
        m = morlet(freqVec(j),width,Fs,normalization);
    elseif strcmp(type,'STFT')
        m = stft(freqVec(j),width,Fs,normalization);
    else
        error('type must be Morlet or STFT');
    end
    n = length(xr)+length(m)-1;
    y = ifft(fft(xr,n).*fft(m,n));
    y = y(ceil(length(m)/2):length(y)-floor(length(m)/2));
    Wx(j,:) = y(npnts+1:end-npnts);
end

function m = morlet(f,width,Fs,normalization)

dt = 1/Fs;
sf = f/width;
st = 1/(2*pi*sf);
t = -3.5*st:dt:3.5*st;
if strcmp(normalization,'energy')
    A = sqrt(dt)*(st*sqrt(pi))^(-0.5); %classic energy-normalized scaling
elseif strcmp(normalization,'area')
    A = 2*dt/(st*sqrt(2*pi)); %better (for visualization) area-normalized scaling
else
    error('normalization must be energy or area');
end
m = A*exp(-t.^2/(2*st^2)).*exp(1i*2*pi*f.*t);

function m = stft(f,width,Fs,normalization)

dt = 1/Fs;
st = width;
t = -3.5*st:dt:3.5*st;
if strcmp(normalization,'energy')
    A = sqrt(dt)*(st*sqrt(pi))^(-0.5); 
elseif strcmp(normalization,'area')
    A = 2*dt/(st*sqrt(2*pi));
else
    error('normalization must be energy or area');
end
m = A*exp(-2*pi*t.^2/(2*st)).*exp(1i*2*pi*f.*t);

function g = gaussian(f,width,Fs)

dt = 1/Fs;
sf = f/width;
st = 1/(2*pi*sf);
t = -3.5*st:dt:3.5*st;
A = 1/(ones(size(t))*exp(-t.^2/(2*st^2))');
g = A*exp(-t.^2/(2*st^2));