function [P, I] = myMatrixProfile(T, m)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%a vectorized version of STOMP%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = length(T);
if n == size(T, 2)
    T = T';
end

%set vars
l = n - m + 1;
P = inf(l, 1);
I = zeros(l, 1);

%compute first row of QT and D
T(n + 1:(m + n)) = 0;
Tfft = fft(T);
TCumsum = cumsum(T);
T2Cumsum =  cumsum(T .^ 2);
T2Sum = T2Cumsum(m:n) - [0; T2Cumsum(1:n - m)];
TSum = TCumsum(m:n) - [0; TCumsum(1:n - m)];
TMu = TSum ./ m;
T2Sig = (T2Sum ./ m) - (TMu .^ 2);
TSig = sqrt(T2Sig);
Q = T(1:m);
Q = Q(end:-1:1);
Q(m+1:(m+n)) = 0;
Qfft = fft(Q);
QT = ifft(Tfft .* Qfft); 
QT = QT(m:n);
D = sqrt(abs(2*(m-(QT-m*TMu*TMu(1))./(TSig * TSig(1)))));

%update matrix profile/index
[P(1), I(1)] = min(D(1+m:end));
I(1) = I(1) + m;
I(1+m:end) = 1;
P(1+m:end) = D(1+m:end);

%main STOMP
QT = QT(1+m:end);
for i = 2:l-m
    
    %get new QT
    QT = QT(1:end-1) - T(i-1)*T(i-1+m:l-1) + T(i-1+m)*T(i-1+2*m:l-1+m);
    
    %get new distance profile
    D = sqrt(abs(2*(m-(QT-m*TMu(i+m:l)*TMu(i))./(TSig(i+m:l) * TSig(i)))));
    
    %update matrix profile/index
    [mD, midx] = min(D);
    if mD < P(i)
        P(i) = mD;
        I(i) = midx+m+i-1;
    end    
    updatePos = D < P(i+m:l);
    I(logical([zeros(m+i-1,1);updatePos])) = i;
    P(logical([zeros(m+i-1,1);updatePos])) = D(updatePos);
    
end




    
  
