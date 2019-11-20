function t = flippos(p); 

c1 = [38 192];
c2 = [132 66];
c3 = [195 209];
plot(c1(1),c1(2),'ro');
hold on
plot(c2(1),c2(2),'ro');
plot(c3(1),c3(2),'ro');
t = [];
for j=1:length(p)/1000; 
    clf       
    plot(c1(1),c1(2),'ro');
    hold on
    plot(c2(1),c2(2),'ro');
    plot(c3(1),c3(2),'ro');
    i1 = (j-1)*1000+1;
    i2 = j*1000;
    fprintf(1,'%d %d\n',i1,i2); 
    plot(p(i1:i2,2),p(i1:i2,3),'.'); 
    answer = input(' ','s');
    if answer == '-'
        fprintf(1,'Reject!\n');
    else
        fprintf(1,'Accept!\n');
        t = [t [p(i1,1) p(i2,1)]']; 
    end
end


