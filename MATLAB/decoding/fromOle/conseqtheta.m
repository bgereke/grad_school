function [xvo,tvo,pvo,svo]=conseqtheta(xv,tv,pv,sv)

xvo = [];
pvo = [];
tvo = [];
svo = [];
xvr = [];
pvr = [];
tvr = [];
theta_thres =4 
theta_per = 0.167; 
for k=1+theta_thres:length(xv)
    if xv(k) - xv(k-theta_thres) < theta_thres*theta_per
        xvo = [xvo xv(k)];
        tvo = [tvo tv(k)];
        pvo = [pvo pv(k)];
        svo = [svo sv(k)];
    else
        xvr = [xvr xv(k)];
        tvr = [tvr tv(k)];
        pvr = [pvr pv(k)];

    end
end
