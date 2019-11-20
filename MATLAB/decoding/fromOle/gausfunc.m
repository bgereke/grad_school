function y = gausfunc(x,m,d)


y = exp(-(1/(2*d^2)).*((x-m).^2))./(d*sqrt(2*pi));
