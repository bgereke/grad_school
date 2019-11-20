function returnmap(x)



dx = diff(x);

for i=1:length(dx)-1
    plot(dx(i),dx(i+1)); 
    hold on
end
