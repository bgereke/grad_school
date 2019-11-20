function tf5plot(c1,c2,c3,c4,c5,xmin,xmax) 

w = 100;       
c1 =c1(:,find(c1(1,:) <= xmax & c1(1,:) >= xmin));
c2 =c2(:,find(c2(1,:) <= xmax & c2(1,:) >= xmin));
c3 =c3(:,find(c3(1,:) <= xmax & c3(1,:) >= xmin));
c4 =c4(:,find(c4(1,:) <= xmax & c4(1,:) >= xmin));
c5 =c5(:,find(c5(1,:) <= xmax & c5(1,:) >= xmin));


[t1,cs1] = interpsmooth(c1,w);
[t2,cs2] = interpsmooth(c2,w);
[t3,cs3] = interpsmooth(c3,w);
[t4,cs4] = interpsmooth(c4,w);
[t5,cs5] = interpsmooth(c5,w);


plot(t1,cs1,'b-'),
hold on
plot(t2,cs2,'g-'),
plot(t3,cs3,'r-'),
plot(t4,cs4,'c-'),
plot(t5,cs5,'m-'),
legend('1','2','3','4','5'); 
