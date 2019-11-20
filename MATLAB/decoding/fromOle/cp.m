function cp

load pos
load velosmooth
load nall
load fall
load location


sig = 10;
dt = 0.1;
z = rec_bay2(p,vs,nall,fall,loc,sig,dt);

sig = 12;
dt = 0.12;
z = rec_bay2(p,vs,nall,fall,loc,sig,dt);

sig = 14;
dt = 0.14;
z = rec_bay2(p,vs,nall,fall,loc,sig,dt);

sig = 16;
dt = 0.16;
z = rec_bay2(p,vs,nall,fall,loc,sig,dt);

sig = 18;
dt = 0.18;
z = rec_bay2(p,vs,nall,fall,loc,sig,dt);

sig = 20;
dt = 0.2;
z = rec_bay2(p,vs,nall,fall,loc,sig,dt);

sig = 22;
dt = 0.22;
z = rec_bay2(p,vs,nall,fall,loc,sig,dt);

sig = 24;
dt = 0.24;
z = rec_bay2(p,vs,nall,fall,loc,sig,dt);

sig = 26;
dt = 0.26;
z = rec_bay2(p,vs,nall,fall,loc,sig,dt);

sig = 28;
dt = 0.28;
z = rec_bay2(p,vs,nall,fall,loc,sig,dt);

sig = 30;
dt = 0.30;
z = rec_bay2(p,vs,nall,fall,loc,sig,dt);


