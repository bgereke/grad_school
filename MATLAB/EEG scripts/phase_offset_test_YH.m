clear;close all
cd('G:\R69');
list=textread('analyzelist.txt','%s','headerlines',0);
clock;
disp(ans(1,1:5))
f1=25; %low pass slow gamma
f2=55; %high pass slow gamma
f3=60; %low pass fast gamma
f4=100; %high pass fast gamma
f5=6;f6=12;%theta passes
Fs=2000;
j=0;
for i=1:size(list,1);
    path=char(list(i))
    cd (path);cd('begin1');
    for DGCSClist=[2]; %your DG CSC numbers
        
        file=strcat(path,'\begin1\CSC',num2str(DGCSClist),'.ncs')
        [samples1,ts,tt, Fs, bv, ir] = loadEEG2(file);
        samples1=bv.*samples1;
        for CA3CSClist=[6]; %your CA3 CSC numbers
            file=strcat(path,'\begin1\CSC',num2str(CA3CSClist),'.ncs')
            [samples2,ts,tt, Fs, bv, ir] = loadEEG2(file);
            samples2=bv.*samples2;
            params.Fs=2000; % sampling frequency
            params.fpass=[f1 f2]; % band of frequencies to be kept
            params. tapers=[3 5]; % taper parameters
            params.pad=2; % pad factor for fft
            params.err=[2 0.05];
            params.trialave=1;
            movingwin=[0.5 0.05];
            [C,phisl,S12,S1,S2,f,confC,phierr,Cerr]=coherencyc(samples1,samples2,params);
            params.fpass=[f3 f4]; % band of frequencies to be kept
            [C,phifs,S12,S1,S2,f,confC,phierr,Cerr]=coherencyc(samples1,samples2,params);
            params.fpass=[f5 f6]; % band of frequencies to be kept
            [C,phith,S12,S1,S2,f,confC,phierr,Cerr]=coherencyc(samples1,samples2,params);
            j=j+1;
            [b1x1(j,:),b1y1(j,:)] = rose((phisl),36);
            [b1x2(j,:),b1y2(j,:)] = rose((phifs),36);
            [b1x3(j,:),b1y3(j,:)] = rose((phith),36);
        end
    end
    j=0;cd('../');cd('begin2');
    for DGCSClist=[2]; %your DG CSC numbers
        
        file=strcat(path,'\begin2\CSC',num2str(DGCSClist),'.ncs')
        [samples1,ts,tt, Fs, bv, ir] = loadEEG2(file);
        samples1=bv.*samples1;
        for CA3CSClist=[6]; %your CA3 CSC numbers
            file=strcat(path,'\begin2\CSC',num2str(CA3CSClist),'.ncs')
            [samples2,ts,tt, Fs, bv, ir] = loadEEG2(file);
            samples2=bv.*samples2;
            params.Fs=2000; % sampling frequency
            params.fpass=[f1 f2]; % band of frequencies to be kept
            params. tapers=[3 5]; % taper parameters
            params.pad=2; % pad factor for fft
            params.err=[2 0.05];
            params.trialave=1;
            movingwin=[0.5 0.05];
            [C,phisl,S12,S1,S2,f,confC,phierr,Cerr]=coherencyc(samples1,samples2,params);
            params.fpass=[f3 f4]; % band of frequencies to be kept
            [C,phifs,S12,S1,S2,f,confC,phierr,Cerr]=coherencyc(samples1,samples2,params);
            params.fpass=[f5 f6]; % band of frequencies to be kept
            [C,phith,S12,S1,S2,f,confC,phierr,Cerr]=coherencyc(samples1,samples2,params);
            j=j+1;
            [b2x1(j,:),b2y1(j,:)] = rose((phisl),36);
            [b2x2(j,:),b2y2(j,:)] = rose((phifs),36);
            [b2x3(j,:),b2y3(j,:)] = rose((phith),36);
            
        end
    end
    j=0;cd('../');cd('begin3');
    for DGCSClist=[2]; %your DG CSC numbers
        
        file=strcat(path,'\begin3\CSC',num2str(DGCSClist),'.ncs')
        [samples1,ts,tt, Fs, bv, ir] = loadEEG2(file);
        samples1=bv.*samples1;
        for CA3CSClist=[6]; %your CA3 CSC numbers
            file=strcat(path,'\begin3\CSC',num2str(CA3CSClist),'.ncs')
            [samples2,ts,tt, Fs, bv, ir] = loadEEG2(file);
            samples2=bv.*samples2;
            params.Fs=2000; % sampling frequency
            params.fpass=[f1 f2]; % band of frequencies to be kept
            params. tapers=[3 5]; % taper parameters
            params.pad=2; % pad factor for fft
            params.err=[2 0.05];
            params.trialave=1;
            movingwin=[0.5 0.05];
            [C,phisl,S12,S1,S2,f,confC,phierr,Cerr]=coherencyc(samples1,samples2,params);
            params.fpass=[f3 f4]; % band of frequencies to be kept
            [C,phifs,S12,S1,S2,f,confC,phierr,Cerr]=coherencyc(samples1,samples2,params);
            params.fpass=[f5 f6]; % band of frequencies to be kept
            [C,phith,S12,S1,S2,f,confC,phierr,Cerr]=coherencyc(samples1,samples2,params);
            j=j+1;
            [b3x1(j,:),b3y1(j,:)] = rose((phisl),36);
            [b3x2(j,:),b3y2(j,:)] = rose((phifs),36);
            [b3x3(j,:),b3y3(j,:)] = rose((phith),36);
            
        end
    end
        cd('../');
        x1=[b1x1;b2x1;b3x1];x2=[b1x2;b2x2;b3x2];x3=[b1x3;b2x3;b3x3];
        y1=[b1y1;b2y1;b3y1];y2=[b1y2;b2y2;b3y2];y3=[b1y3;b2y3;b3y3];
        Mx1=mean(x1);My1=mean((y1/sum(b1y1(1,:))));Mx2=mean(x2);My2=mean((y2/sum(b1y2(1,:)))); 
        Mx3=mean(x3);My3=mean((y3/sum(b1y3(1,:))));%normalize y2 counts to tt points
        close all;
        figure;
        rose3=polar(Mx1,My1);
        set(rose3,'LineWidth',3,'Color','b');
        title('Normalized DG-CA3 slow gamma phase offset')
        
        figure;
        rose3=polar(Mx2,My2);
        set(rose3,'LineWidth',3,'Color','r');
        title('Normalized DG-CA3 fast gamma phase offset')
        
        figure;
        rose3=polar(Mx3,My3);
        set(rose3,'LineWidth',3,'Color',[0 .4 0]);
        title('Normalized DG-CA3 theta phase offset')
        
        saveas(1, 'slow_phase_offset', 'fig');saveas(1, 'slow_phase_offset', 'jpg')
        saveas(2, 'fast_phase_offset', 'fig');saveas(2, 'fast_phase_offset', 'jpg')
        saveas(3, 'theta_phase_offset', 'fig');saveas(3, 'theta_phase_offset', 'jpg')
    end
    
    


clock;
disp(ans(1,1:5))


%