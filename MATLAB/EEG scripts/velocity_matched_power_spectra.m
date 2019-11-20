%%%%% Velocity Matched Power Spectrum %%%%%

velbins = 5:25;
pre = 15;
post_start_one = 33;
post_end_one = 43;
% post_start_two = 50;
% post_end_two = 80;
TFRpre = []; TFRpost_one = [];% TFRpost_two = [];

for vb = 1:size(velbins,2)-1
    
    %cumpute number of velocity samples for each period and find min period
    vidx = find(vel >= velbins(vb) & vel < velbins(vb+1));
    numvels_pre = sum(w(vidx,2)/60<=pre);
    numvels_post_one = sum(w(vidx,2)/60>=post_start_one & w(vidx,2)/60<=post_end_one);
%     numvels_post_two = sum(w(vidx,2)/60>=post_start_two & w(vidx,2)/60<=post_end_two);
    
    minvels = min([numvels_pre numvels_post_one]); %numvels_post_two]);
    
    %pick random subsample from each period
    pre_ss = randsample(numvels_pre,minvels);
    post_one_ss = randsample(numvels_post_one,minvels);
%     post_two_ss = randsample(numvels_post_two,minvels);
    
    %take the subsample
    w_pre = widx(w(vidx,2)/60<=pre,:);
    w_pre = w_pre(pre_ss,:);
    w_post_one = widx(w(vidx,2)/60>=post_start_one & w(vidx,2)/60<=post_end_one,:);
    w_post_one = w_post_one(post_one_ss,:);
%     w_post_two = widx(w(vidx,2)/60>=post_start_two & w(vidx,2)/60<=post_end_two,:);
%     w_post_two = w_post_two(post_two_ss,:);    
    
    for ss = 1:minvels
       TFRpre = [TFRpre TFR(:,w_pre(ss,1):w_pre(ss,2))];
       TFRpost_one = [TFRpost_one TFR(:,w_post_one(ss,1):w_post_one(ss,2))];
%        TFRpost_two = [TFRpost_two TFR(:,w_post_two(ss,1):w_post_two(ss,2))];
    end
    
end

meanTFRpre = mean(log(TFRpre),2);
seTFRpre = std(log(TFRpre),0,2)/sqrt(size(TFRpre,2)-1);
meanTFRpost_one = mean(log(TFRpost_one),2);
seTFRpost_one = std(log(TFRpost_one),0,2)/sqrt(size(TFRpost_one,2)-1);
% meanTFRpost_two = mean(log(TFRpost_two),2);
% seTFRpost_two = std(log(TFRpost_two),0,2)/sqrt(size(TFRpost_two,2)-1);

plot(freqvec,meanTFRpre,'-k'); hold on
plot(freqvec,meanTFRpost_one,'-b');
% plot(freqvec,meanTFRpost_two,'-r');

xlabel('frequency (Hz)');ylabel('log(power)');title('power spectra equalized for running speed')
legend('first 15 min saline','33-48 min saline')

plot(freqvec,meanTFRpre+seTFRpre,'--k')
plot(freqvec,meanTFRpre-seTFRpre,'--k')

plot(freqvec,meanTFRpost_one+seTFRpost_one,'--b')
plot(freqvec,meanTFRpost_one-seTFRpost_one,'--b')

% plot(freqvec,meanTFRpost_two+seTFRpost_two,'--r')
% plot(freqvec,meanTFRpost_two-seTFRpost_two,'--r')