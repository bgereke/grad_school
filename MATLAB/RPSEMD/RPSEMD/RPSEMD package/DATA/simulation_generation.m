% the simulation

r=1200;     %%%%length
t=meshgrid(1:r,1:r); %meshgrid

f=4.5/1200;
phi_f=2*pi*f*t;
Ih1=sin(phi_f);    %%% the IM1

Ih2=1/10*cos(phi_f*16+5*2*pi/16);   %%% the IM2, small oscillations
noi=Ih2(1,:);   
noi(1:292)=0;
noi(374:559)=0;
noi(641:826)=0;
noi(908:end)=0;

Ih3=1/5*cos(phi_f*4+5*2*pi/8);   %%%% the IM3
noi1=Ih3(1,:);
noi1(1:429)=0;
noi1(530:687)=0;
noi1(688:953)=0;
noi1(1038:end)=0;

Y=Ih1(1,:)+noi+noi1;   %%%% the composed signal

figure;plot(Y)
