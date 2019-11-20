%unpaires specs
numboots = 10000;
onboots = zeros(numboots, 50);
offboots = zeros(numboots, 50);

for b = 1:numboots   
    ridx = randsample(size(spec,2),size(spec,2),'true');
    onboots(b,:) = mean(spec(1:50,ridx),2);
    offboots(b,:) = mean(spec(51:100,ridx),2);   
end

[onci] = prctile(onboots,[2.5 97.5]);
[offci] = prctile(offboots,[2.5 97.5]);

p1=plot(2:2:100,log10(mean(spec(51:100,:),2)),'-k','LineWidth',2);hold on
plot(2:2:100,log10(offci)','-k','LineWidth',0.5)
p2=plot(2:2:100,log10(mean(spec(1:50,:),2)),'-','LineWidth',2,'Color',[0 0 0.75]);
plot(2:2:100,log10(onci)','-','Color',[0 0 0.75],'LineWidth',0.5)
xlim([2,100])
xlabel('Frequency (Hz)')
ylabel('log(power)')
legend([p1 p2],'light off','light on')
title('Unpaired - Mouse 77 GACR2')
axis square
txt = '1 mouse, 6 days';
text(64, 3.25, txt,'FontSize',12)

%paired specs
numboots = 10000;
boots = zeros(numboots, 50);

for b = 1:numboots   
    ridx = randsample(size(spec,2),size(spec,2),'true');
    boots(b,:) = (mean(spec(1:50,ridx),2)-mean(spec(51:100,ridx),2))./mean(spec(51:100,ridx),2)*100;  
end

[ci] = prctile(boots,[2.5 97.5]);

figure
plot([2 100],[0 0],'--k');hold on
shadedplot(2:2:100,ci(2,:),ci(1,:),[0 0 0.5],[1 1 1]);alpha(0.5);hold on
plot(2:2:100,(mean(spec(1:50,:),2)-mean(spec(51:100,:),2))./mean(spec(51:100,:),2)*100,'-','LineWidth',2,'Color',[0 0 0.75]);
% plot(2:2:100,ci','-','Color',[0 0.75 0],'LineWidth',0.5)
xlim([2,100])
xlabel('Frequency (Hz)')
ylabel('% change from light off')
title('Paired Differences - Mouse 77 GACR2')
axis square
txt = '1 mouse, 6 days';
text(64, 10, txt,'FontSize',12)


%unpaired ONS/OFFS
numboots = 10000;
onboots = zeros(numboots, 50);
offboots = zeros(numboots, 50);

for b = 1:numboots   
    onidx = randsample(size(ONS,2),size(ONS,2),'true');
    offidx = randsample(size(OFFS,2),size(OFFS,2),'true');
    onboots(b,:) = mean(ONS(:,onidx),2);
    offboots(b,:) = mean(OFFS(:,offidx),2);   
end

[onci] = prctile(onboots,[2.5 97.5]);
[offci] = prctile(offboots,[2.5 97.5]);

p1=plot(2:2:100,log10(mean(OFFS,2)),'-k','LineWidth',2);hold on
plot(2:2:100,log10(offci)','-k','LineWidth',0.5)
p2=plot(2:2:100,log10(mean(ONS,2)),'-','LineWidth',2,'Color',[0.3 0.3 1]);
plot(2:2:100,log10(onci)','-','Color',[0.3 0.3 1],'LineWidth',0.5)
xlim([2,100])
xlabel('Frequency (Hz)')
ylabel('log(power)')
legend([p1 p2],'light off','light on')
title('Unpaired - Mouse 88 GACR2 MEC')
axis square
txt = '1 mouse, 2 days';
text(64, 3.25, txt,'FontSize',12)

%unpaired fields on/off
numboots = 10000;
onboots = zeros(numboots, 60);
offboots = zeros(numboots, 60);
nFields = Fields./max(Fields(:,61:120),[],2);

for b = 1:numboots   
    bidx = randsample(size(Fields,1),size(Fields,1),'true');
    onboots(b,:) = mean(nFields(bidx,61:120));
    offboots(b,:) = mean(nFields(bidx,1:60));   
end

[onci] = prctile(onboots,[2.5 97.5]);
[offci] = prctile(offboots,[2.5 97.5]);

p1=plot(linspace(0,1,60),mean(nFields(:,1:60)),'-k','LineWidth',2);hold on
plot(linspace(0,1,60),offci','-k','LineWidth',0.5)
p2=plot(linspace(0,1,60),mean(nFields(:,61:120)),'-','LineWidth',2,'Color',[0.3 0.3 1]);
plot(linspace(0,1,60),onci','-','Color',[0.3 0.3 1],'LineWidth',0.5)
xlim([0,1])
xlabel('position (normalized)')
ylabel('firing rate (normalized))')
legend([p1 p2],'light off','light on')
