
efile = 'C:\Data\mouse23\2014-09-29_10-21-24\begin1\Events.nev';
FS = [1 0 1 0 1]; EH = 0; EM = 1;
[TimeStamps, TTLs, EventStrings] = Nlx2MatEV(efile,FS,EH,EM);
TTLs(1) = []; TimeStamps(1) = []; TimeStamps = TimeStamps/1000000*1000; %milliseconds

file1 = 'C:\Data\mouse23\2014-09-29_10-21-24\begin1\CSC5.ncs';
file2 = 'C:\Data\mouse23\2014-09-29_10-21-24\begin1\CSC9.ncs';
gfile = 'C:\Data\mouse23\2014-09-29_10-21-24\begin1\CSC13.ncs';

[samples,ts,tt, Fs, bv, ir] = loadEEG2(file1);
ch1 = bv*samples;
tt = tt*1000; %milliseconds
tt = tt-min(tt);
TimeStamps = TimeStamps - min(tt);

[samples,ts,tt, Fs, bv, ir] = loadEEG2(file2);
ch2 = bv*samples;
tt = tt*1000; %milliseconds
tt = tt-min(tt);
TimeStamps = TimeStamps - min(tt);

[samples,ts,tt, Fs, bv, ir] = loadEEG2(gfile);
gch = bv*samples;
tt = tt*1000; %milliseconds
tt = tt-min(tt);
TimeStamps = TimeStamps - min(tt);

ch1 = ch1+gch; ch2 = ch2+gch;
[ImX] = traces2ImX(ch1,ch2,[9 30 70],Fs,7); 

w = 2000;
idx = [8];

for i = idx % or whatever, depending on the length of your recording
    subplot(4,1,1);hold off
    plot(tt,ch1,'-r','LineWidth',1);hold on
    plot(tt,ch2,'-b','LineWidth',1);
    xlim([w*(i-1) w*i]); % scroll through your recording to find a good example, one second at a time
%     ylim([-max(abs([ch1(tt>w*(i-1)&tt<w*i);ch2(tt>w*(i-1)&tt<w*i)]))-300 max(abs([ch1(tt>w*(i-1)&tt<w*i);ch2(tt>w*(i-1)&tt<w*i)]))])
    ylim([-max(abs([ch1(tt>w*(i-1)&tt<w*i);ch2(tt>w*(i-1)&tt<w*i)]))-300 max(abs([ch1(tt>w*(i-1)&tt<w*i);ch2(tt>w*(i-1)&tt<w*i)]))])
    hold on
    y = -max(abs([ch1(tt>w*(i-1)&tt<w*i);ch2(tt>w*(i-1)&tt<w*i)]))-200; x = (w*(i-1)+w*i)/2;
    line([x-200 x], [y y], 'color','k');  % x calibration bar- get the values from your figure depending on the length and location you want- you can move it later in Illustrator
    line([x x], [y y+250], 'color','k');  % y calibration bar
    
    text([x-125], [y-75], '200 ms');  % label for your x calibration bar
    text([x-20], [y+350], '250 uV');  % label for your y calibration bar
    title(i);axis off;  % which second-long trace you have plotted
    
    subplot(4,1,2);hold off
    if ~isempty(tt(tt>w*(i-1)&tt<w*i&ImX(1,:)>0))
        shadedplot(tt(tt>w*(i-1)&tt<w*i&ImX(1,:)>0), ImX(1,tt>w*(i-1)&tt<w*i&ImX(1,:)>0),zeros(size(ImX(1,tt>w*(i-1)&tt<w*i&ImX(1,:)>0))),[1 0 0],[0 0 0]); alpha(1);hold on
    end
    if ~isempty(tt(tt>w*(i-1)&tt<w*i&ImX(1,:)<=0))
        shadedplot(tt(tt>w*(i-1)&tt<w*i&ImX(1,:)<=0),zeros(size(ImX(1,tt>w*(i-1)&tt<w*i&ImX(1,:)<=0))), ImX(1,tt>w*(i-1)&tt<w*i&ImX(1,:)<=0),[0 0 1],[0 0 0]); alpha(1);hold on
    end
    xlim([w*(i-1) w*i]); % scroll through your recording to find a good example, one second at a time
%     ylim([-max(abs(ImX(1,tt>w*(i-1)&tt<w*i))) max(abs(ImX(1,tt>w*(i-1)&tt<w*i)))])
    ylim([-max(abs(ImX(1,:))) max(abs(ImX(1,:)))])
    axis off;
    
    subplot(4,1,3);hold off
    if ~isempty(tt(tt>w*(i-1)&tt<w*i&ImX(2,:)>0))
        shadedplot(tt(tt>w*(i-1)&tt<w*i&ImX(2,:)>0), ImX(2,tt>w*(i-1)&tt<w*i&ImX(2,:)>0),zeros(size(ImX(2,tt>w*(i-1)&tt<w*i&ImX(2,:)>0))),[1 0 0],[0 0 0]); alpha(1);hold on
    end
    if ~isempty(tt(tt>w*(i-1)&tt<w*i&ImX(2,:)<=0))
        shadedplot(tt(tt>w*(i-1)&tt<w*i&ImX(2,:)<=0),zeros(size(ImX(2,tt>w*(i-1)&tt<w*i&ImX(2,:)<=0))), ImX(2,tt>w*(i-1)&tt<w*i&ImX(2,:)<=0),[0 0 1],[0 0 0]); alpha(1);hold on
    end
    xlim([w*(i-1) w*i]); % scroll through your recording to find a good example, one second at a time
%     ylim([-max(abs(ImX(2,tt>w*(i-1)&tt<w*i))) max(abs(ImX(2,tt>w*(i-1)&tt<w*i)))])
    ylim([-max(abs(ImX(2,:))) max(abs(ImX(2,:)))])
    axis off;
    
    subplot(4,1,4);hold off
    if ~isempty(tt(tt>w*(i-1)&tt<w*i&ImX(3,:)>0))
        shadedplot(tt(tt>w*(i-1)&tt<w*i&ImX(3,:)>0), ImX(3,tt>w*(i-1)&tt<w*i&ImX(3,:)>0),zeros(size(ImX(3,tt>w*(i-1)&tt<w*i&ImX(3,:)>0))),[1 0 0],[0 0 0]); alpha(1);hold on
    end
    if ~isempty(tt(tt>w*(i-1)&tt<w*i&ImX(3,:)<=0))
        shadedplot(tt(tt>w*(i-1)&tt<w*i&ImX(3,:)<=0),zeros(size(ImX(3,tt>w*(i-1)&tt<w*i&ImX(3,:)<=0))), ImX(3,tt>w*(i-1)&tt<w*i&ImX(3,:)<=0),[0 0 1],[0 0 0]); alpha(1);hold on
    end
    xlim([w*(i-1) w*i]); % scroll through your recording to find a good example, one second at a time
%     ylim([-max(abs(ImX(3,tt>w*(i-1)&tt<w*i))) max(abs(ImX(3,tt>w*(i-1)&tt<w*i)))])
    ylim([-max(abs(ImX(3,:))) max(abs(ImX(3,:)))])
    axis off;
    keyboard  % hit any key to proceed to the next trace
    
end

% after you have picked the trace you want, plot it
i=1;
plot(tt,ch1,'k');
xlim([w*(i-1) w*i]);
%xlabel('time (ms)')
hold on
line([w*(i-1) w*(i-1)+100], [-300 -300], 'color','k');  % x calibration bar- get the values from your figure depending on the length and location you want- you can move it later in Illustrator
line([w*(i-1)+100 w*(i-1)+100], [-300 -100], 'color','k');  % y calibration bar

text([w*(i-1)], [-300], '200 ms');  % label for your x calibration bar
text([w*(i-1)+100], [-100], '200 uV');  % label for your y calibration bar


axis off;  


% then just save it as an .eps file, and you can work on it in
% illustrator