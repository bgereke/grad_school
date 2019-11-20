% RPSEMD: 
% the outer loop of RPSEMD to decompose signal Y into a series of pure IMFs
%
% By W.Chenxing, 2015.11.7, at Singapore. w.chenxing@gmail.com. Copyright reserved.

function imf_rps = Rpsemd(Y)
 
% initialization                    
ext_f_max_init = 0.5;            % the threshold to filter the invalid high-frequency extrema
f_sin_pre = 0.75*ones(1,2);      % the threshold to filter the invalid high-frequency cluster
imf_rps = [];                    % the final pure IMFs
nimf = 1;

imf = emd(Y);                                    % initial IMFs
% k=size(imf,1);  figure;  for i=1:k   subplot(k,1,i);  plot(imf(i,:));  ylabel(['IMF' num2str(i)]); end;

% get f_sin and A_sin to design the sinusoid; get the fixed number of phase shifting 
[ f_sin, A_sin, nshift, ext_f_max] = Design_Asin_Fsin(imf, ext_f_max_init, f_sin_pre, nimf);     
% if nshift>10
%     nshift = 10;
% end

while 1
    % if no need of iteration
    if f_sin==1   imfr=emd(Y);  break;  end

    % produce a pure final IMF assisted by phase shifted sinusoid
    [imf_sin, Y] = Phaseshift(Y, A_sin(1), f_sin(1), nshift);  % figure,plot(imf_sin)
    
    % range the final IMFs
    imf_rps = [imf_rps; imf_sin];   
    nimf = nimf + 1;

    % update the threshold used in "Design_Asin_Fsin"
    if length(f_sin) == 1   f_sin(2) = 0;  end
    f_sin_pre = f_sin;      
    
    % identify whether the iteration stops
    [imax, imin, iminmax]=Detect_extrema(Y);        
    if length(iminmax) <= 4
        imfr = emd(Y);
        if imfr == 0   imfr = [];  end
        break;
    end
     
    imf = emd(Y);
   
    % get f_sin and A_sin to design a new sinusoid sk(t)
    [f_sin, A_sin, ns, ext_f_max] = Design_Asin_Fsin(imf, ext_f_max, f_sin_pre, nimf);       
end

imf_rps = [imf_rps; imfr];    % range the final IMFs
end






% PHASESHIFT: produce a pure final IMF through phase shift the sinusoid
%
% OUTPUT:      
% imf_sin,imf_res: a final pure IMF and the corresponding residue
%
% INPUT:           
%           Y: original signal
% A_sin,f_sin: the amplitude and frequency to design a auxiliary sinusoid
%        Npha: the number of phase shifting 

function [imf_sin, imf_res] = Phaseshift( Y, A_sin, f_sin, Npha )

L=length(Y);  imf_sin=0;  imf_res=0;
 
for i = 1 : Npha, 
    S = 0;
    S = A_sin*sin( 2*pi*f_sin*[0:L-1] + (i-1)*2*pi/Npha );   % the designed sinusoid    
   
    Y1 = Y + S;                          % the new signal
     
    imf = emd( Y1 );                     % IMFs of the new sigal
    
    imf_sin = imf_sin + imf(1,:);        % the summation of all IMF1
      
    imf_res = imf_res + (Y1-imf(1,:));   % the summation of all residues
end

% the final IMFk after average of the results in all phase shifting
imf_sin = imf_sin/Npha;  

% the residue
imf_res = imf_res/Npha;  
end
