function A = calcium(A,dt)

tau = 0.020;     
A = A * exp(-dt/tau);
A = A + 1; 

