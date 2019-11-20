%function [result] = doall

%be in a session

%get fields that will be used to score the sequences
%[f,~,~,~] = getfields('TTDecodeList.txt',2);

%get all the theta windows
[w] = thetawindows('TTDecodeList.txt');

%convert w into cell scores
scores = cell(size(w,1),2);
for p = 1:size(w,1)
       scores{p,2} = w(p,:);
end

%get sequence scores for the windows using the fields
%[scores] = sequencewindows(f,w,3); %must have 3 different cells to score it

%decode all the windows
[scores] = decodewindows('TTDecodeList.txt',scores,w);

%characterize the decoded windows
 for i = 1:size(scores,1)
% [scores{i,16}] = getbestline(scores{i,12});
% if ~isempty(scores{i,16})
% [scores{i,17}] = scores{i,14}(scores{i,16}(1,2)) - scores{i,14}(scores{i,16}(end,2)); %pathlength
% best = scores{i,14}(scores{i,16}(:,2));
% actualind = round(scores{i,13}(scores{i,16}(:,1)));
% actualind(actualind > size(scores{1,14}(:,1),1)) = size(scores{1,14}(:,1),1);
% actualind(actualind<1) = 1;
% actual = scores{i,14}(actualind);
% scores{i,18} = best - actual;
% fit = polyfit(1:length(scores{i,18}),scores{i,18}',1);
% scores{i,19} = fit(1,1);
% R = (corrcoef(1:length(scores{i,18}),scores{i,18}')).^2;
% scores{i,20} = R(1,2);
%%% % column 21 is the normalized decoded so each time window is maxed at 1
mx = nansum(scores{i,12}); %adjust so that each time bin, the sum of P(x|n) (for all x) = 1;
scores{i,21} = scores{i,12};
for o = 1:size(scores{i,12},2)
scores{i,21}(:,o) = scores{i,21}(:,o)/mx(o);
% if mx(o) == 0
%     scores{i,21}(:,o) = 0; % so not dividing by zero and getting Nana
% end
end
%%%%% for getting cumalitvie score - see scriptscore.m
% for g = 1:size(scores{i,21},2) %going through time windows
%     all = 0; 
%     for j = 1:size(scores{i,21},1) %going through each xbin
%          all = all + j*scores{i,21}(j,g);
%     end
%     com = all/sum(scores{i,21}(:,g));
%     allcom (1,g) = round(com);
% end
% allcom'
% floor(scores{i,13})
% for z = 1:length(scores{i,21})
%     if isnan(allcom(z))
%         scores{i,22}(1,z) = 0; %neutral score if no spikes had occurred
%     else
%         scores{i,22}(1,z) = scores{i,14}(allcom(z))-scores{i,14}(floor(scores{i,13}(z)));
%     end
% end
% clear actual best allcom;
end
%end

scriptscore;  %fix to minimize |x* - x| (i.e. take median of posterior), and maybe weight cycle estimate by posterior variances
 
[scores] = bayesGamma_v3(scores,'TTDecodeList.txt'); %fix to calculate whole TFR first instead of individual windows

%scriptscore_plotter;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save('bayes_3_04_01_2_2_15_Pxxu1000_1_uV_tfr_w3_onwindows','scores');
%numbers in filename:
%binsize (cm)
%dt
%step time
%mincell
%minspk
%velocity smoothing
%Pxv or Pvx, unidirectional, window velocity range (set inside probXV.m)
%theta filtered template = 1, not = 0.
%unit eeg is in if nothing, it is in default bits
%if tfr was saved
%whitening of signal before tfr calculated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get gamma info during all the windwos
%[scores] = bayesGamma(scores);

% 
% for i = 1:size(scores,1)
%    if scores{i,11} == 1 && ~isempty(scores{i,17})
%    hold on;
%    plot(scores{i,14}(round(scores{i,13}(1,1)),1),  scores{i,17},'ko'); 
%    end
% end