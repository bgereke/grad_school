function velocity = movePointer(vr)

velocity = [0 0 0 0];
scr = get(0,'screensize');
ptr = get(0,'pointerlocation')-scr(3:4)/2;
velocity(1) = ptr(2)/2*sin(-vr.position(4));
velocity(2) = ptr(2)/2*cos(-vr.position(4));
velocity(4) = ptr(1)/100;

velocity = velocity / 5;