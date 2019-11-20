function r = vm_kernel(x,kappa)

C = 1/(2*pi*besseli(0,kappa));
r = C * exp(kappa*cos(x));