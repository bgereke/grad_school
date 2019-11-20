function w = ext_gam
% function w = ext_gam
%
% w = [time Ttheta gamma_in_theta];
%
% Find the # of gamma cycles per theta cycles
% extracting from the files ttheta and tgamma



fid = fopen('ttheta','r');
T = fscanf(fid,'%e'); 
fclose(fid);

fid = fopen('tgamma','r');
G = fscanf(fid,'%e'); 
fclose(fid);

size(T) ;
size(G) ; 
c = zeros(1,(size(T)-1)); 
size(c)
j = 1; 
for i=1:size(T)-1
    T1 = T(i);
    T2 = T(i+1);
    k = 0;
    while G(j) <= T2
        k = k + 1;
        j = j + 1; 
    end 
    if (mod(i,1000) == 0)
       fprintf(1,'%d\n',i/length(T));
    end
    c(i) = k; 
end

dT = diff(T); 
z = [T(1:length(c)) dT(1:length(c))  c'];

Tmax = 0.167;
Tmin = 0.100;

w = zeros(length(z),3); 
j = 0;
for i=1:length(z)
    if  z(i,2) < Tmax & z(i,2) > Tmin
         j = j + 1; 
         w(j,:) = z(i,:);
         if (mod(i,1000) == 0)
             fprintf(1,'%d\n',i/length(z));
         end
    end
end
w = w(1:j,:); 
