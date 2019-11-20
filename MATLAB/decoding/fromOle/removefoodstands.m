function [ev,xv,tv,pv]=removefoodstand(ev,xv,tv,pv);
s1 = 24;
s2 = 75;
s3 = 125;
ds = 5;

i = find((tv < s1) | (tv > s1+ds & tv < s2) | (tv > s2+ds & tv < s3) | (tv > s3+ds))  ;

tv = tv(i);
pv = pv(i);
xv = xv(i);
ev = ev(i);
