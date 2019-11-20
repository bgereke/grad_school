% Forward continuous wavelet transform, discretized, as described
% in Sec. 4.3.3 of [1] and Sec. IIIA of [2].  This algorithm uses
% the FFT and samples the wavelet atoms in the fourier domain.
% Options such as padding of the original signal are allowed.
% Returns the vector of scales and, if requested, the analytic
% time-derivative of the wavelet transform (as described in
% Sec. IIIB of [2].
%
% 1. Mallat, S., Wavelet Tour of Signal Processing 3rd ed.
%
% 2. G. Thakur, E. Brevdo, N.-S. Fučkar, and H.-T. Wu,
% "The Synchrosqueezing algorithm for time-varying spectral analysis: robustness
%  properties and new paleoclimate applications," Signal Processing, 93:1079-1094, 2013.
% 
% Inputs:
%     x: input signal vector, length n (need not be dyadic length)
%  type: wavelet type, string (see help wfiltfn)
%    nv: number of voices (suggest nv >= 32)
%    dt: sampling period (default, dt = 1)
%   opt: options structure
%    opt.padtype: type of padding, options: 'symmetric',
%                 'replicate', 'circular' (default = 'symmetric')
%    opt.rpadded: return padded Wx and dWx?  (default = 0)
%    opt.type, opt.s, opt,mu, etc: wavelet options (see help wfiltfn)
%
% Outputs:
%    Wx: [na x n] size matrix (rows = scales, cols = times)
%        containing samples of the CWT of x.
%    as: na length vector containing the associated scales
%   dWx (calculated if requested): [na x n] size matrix containing
%       samples of the time-derivatives of the CWT of x.
%	xMean: mean of padded x
%
%---------------------------------------------------------------------------------
%    Synchrosqueezing Toolbox
%    Authors: Eugene Brevdo, Gaurav Thakur
%---------------------------------------------------------------------------------
function [Wx,as,dWx,xMean] = cwt_fw(x, type, nv, dt, opt)
    if nargin<5, opt = struct(); end
    if nargin<4, dt = 1; end

    % Output derivative of Wx too?
    if nargout>2,
        dout = 1;
    else
        dout = 0;
    end

    % options: symmeric, replicate, circular
    if ~isfield(opt, 'padtype'), opt.padtype = 'symmetric'; end
    if ~isfield(opt, 'rpadded'), opt.rpadded = 0; end

    x = x(:); % Turn into column vector
    n = length(x);

    % Pad x first
    [x,N,n1,n2] = padsignal(x, opt.padtype);
%	x([1:n1,n1+n+1:N])=0;
	
	%output mean of padded x
	xMean = mean(x);
	x = x-xMean;

    % Choosing more than this means the wavelet window becomes too short
    noct = log2(N)-1;
    assert(noct > 0 && mod(noct,1) == 0);
    assert(nv>0 && mod(nv,1)==0); % integer
    assert(dt>0);
    assert(~any(isnan(x)));
    
    na = noct*nv;
    as = 2^(1/nv) .^ (1:1:na);
    
    Wx = zeros(na, N);
    dWx = Wx;
	opt.dt = dt;
    
    x = x(:).';
    xh = fft(x);
    
    % for each octave
	% reworked this part to not use wfilth, which slows things down a lot due to
	% branching and temp objects; see that function for more comments
	k = 0:(N-1);
    xi = zeros(1, N);
    xi(1:N/2+1) = 2*pi/N*[0:N/2];
    xi(N/2+2:end) = 2*pi/N*[-N/2+1:-1];
    psihfn = wfiltfn(type, opt);

    for ai = 1:na
        a = as(ai);
		psih = psihfn(a*xi) * sqrt(a) / sqrt(2*pi) .* (-1).^k;
		dpsih = (i*xi/opt.dt) .* psih;

        xcpsi = ifftshift(ifft(psih .* xh));
        Wx(ai, :) = xcpsi;
		dxcpsi = ifftshift(ifft(dpsih .* xh));
		dWx(ai, :) = dxcpsi;
%        if dout,
%            [psih, dpsih] = wfilth(type, N, a, opt);
%            dxcpsi = ifftshift(ifft(dpsih .* xh));
%            dWx(ai, :) = dxcpsi;
%        else
%            psih = wfilth(type, N, a, opt);
%        end
    end
	
%	Wx = Wx + xMean*repmat((as.^1/2).',[1,N]);
	
    % Shorten W to proper size (remove padding)
    if (~opt.rpadded)
        Wx = Wx(:, n1+1:n1+n);
 %       if dout,
            dWx = dWx(:, n1+1:n1+n);
 %       end
    end

    % Output a for graphing purposes, scale by dt
    as = as * dt;
end % cwt_fw
