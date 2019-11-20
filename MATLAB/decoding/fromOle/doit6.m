load velo2.mat
load pos.mat
load ratio.mat
load ttheta

1
load c6_1     
z = extract(c6_1,ttheta,p,v,ratio); 
save c6_1ex.mat z

2
load c6_2     
z = extract(c6_2,ttheta,p,v,ratio); 
save c6_2ex.mat z

3
load c6_3     
z = extract(c6_3,ttheta,p,v,ratio); 
save c6_3ex.mat z

4
load c6_4     
z = extract(c6_4,ttheta,p,v,ratio); 
save c6_4ex.mat z

5
load c6_5     
z = extract(c6_5,ttheta,p,v,ratio); 
save c6_5ex.mat z

6
load c6_6     
z = extract(c6_6,ttheta,p,v,ratio); 
save c6_6ex.mat z

7
load c6_7     
z = extract(c6_7,ttheta,p,v,ratio); 
save c6_7ex.mat z



