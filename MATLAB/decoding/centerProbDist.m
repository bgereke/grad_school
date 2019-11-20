Pf = zeros(600,1); 

for i = 1:size(fPxn,2)
     center = f(i,2);
     temp = fPxn(:,i)/max(fPxn(:,i));
     temp(isnan(temp)) = 0;
     Pf(300-center:300-center+266,1) = Pf(300-center:300-center+266,1) + temp;
end
 
Ps = zeros(600,1); 

for i = 1:size(sPxn,2)
     center = s(i,2);
     temp = sPxn(:,i)/max(sPxn(:,i));
     temp(isnan(temp)) = 0;
     Ps(300-center:300-center+266,1) = Ps(300-center:300-center+266,1) + temp;
end

plot(Ps/max(Ps))
hold on;
plot(Pf/max(Pf),'r')