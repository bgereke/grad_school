function z = makeratio(T,G);
% function z = makeratio(T,G)  ;  
%
% z = [time Ttheta Ttheta/Tgamma];
%
% Find the # of gamma cycles per theta cycles
% extracting from the files ttheta and tgamma


% Acceptable gamma periods
Gmin = 0.0125;
Gmax = 0.0333; 

% Acceptable theta periods
Tmax = 0.167;    
Tmin = 0.100;  

z = zeros(length(T),3); 
j = 2; 
l = 0;  
sum = 0; 
rold = 0; 
for i=1:size(T)-1
    T1 = T(i);
    T2 = T(i+1);
    Tper = T2 - T1; 
    k = 0;

    sum = 0; 
    while G(j) <= T2
        Gper = G(j) - G(j-1); 
        if Gper < Gmax & Gper > Gmin
            sum = sum + Gper;  
            k = k + 1;
        end 
        j = j + 1; 
    end 
    if k ~= 0
         r = sum/k; 
         rold = r; 
    else
         r = rold;  
    end
    if Tper < Tmax % Tper >  Tmin
        l = l + 1;
        z(l,1:3) = [T1 Tper Tper/r]; 
    end 

    if (mod(i,1000) == 0)
       fprintf(1,'%d\n',i/length(T));
    end
end

z = z(1:l,:); 

