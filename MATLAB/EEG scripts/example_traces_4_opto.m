efile = 'C:\Data\mouse27\2015-02-12_09-47-36\begin2\Events.nev';
FS = [1 0 1 0 1]; EH = 0; EM = 1;
[TimeStamps, TTLs, EventStrings] = Nlx2MatEV(efile,FS,EH,EM);
TTLs(1) = []; TimeStamps(1) = []; TimeStamps = TimeStamps/1000000;

% read the file names from the csc-file list
cscid = fopen('C:\Data\mouse27\2015-02-12_09-47-36\cscList.txt','r');
jj = 1;
while ~feof(cscid)
    str = fgetl(cscid);
    channels(jj) = {str};
    jj = jj+1;
end
numchannels = jj-1;
cscid = fclose('all');

%plot full session for each shank
trspace = 3000; %trace spacing in uV
shnum = 1;
for jj = 1:numchannels
    % Load data from the .ncs files, make plots, and store them
    newshnum = str2num(channels{jj}(2));
    if newshnum ~= shnum
        figure
    end
    elnum = str2num(channels{jj}(6:end-4));
    xfile = [cd,'\',channels{jj}];
    [samples,~,tt, Fs, bv, ~] = loadEEG2(xfile);
    x = bv*samples; clear samples
    x = x(1:16:end)-elnum*trspace; %downsample traces to 2000 Hz
    tt = tt(1:16:end);
    mt = min(tt); tt = tt-mt;
    if jj == 1
        TimeStamps = TimeStamps - mt;
    end
    
    %make light vector
    L  = zeros(size(tt)); tidx = 1;
    hold on
    for i = 1:length(TTLs)
        L(tidx:length(tt(tt<=TimeStamps(i)))) = TTLs(i);
        if TTLs(i) ~= 0
            plot(tt(tidx:length(tt(tt<=TimeStamps(i)))),x(tidx:length(tt(tt<=TimeStamps(i)))),'b')
            %scatter(tt(tidx:length(tt(tt<=TimeStamps(i)))),ch_X(tidx:length(tt(tt<=TimeStamps(i)))),'g')
        else
            plot(tt(tidx:length(tt(tt<=TimeStamps(i)))),x(tidx:length(tt(tt<=TimeStamps(i)))),'k')
            %scatter(tt(tidx:length(tt(tt<=TimeStamps(i)))),ch_X(tidx:length(tt(tt<=TimeStamps(i)))),'k')
        end
        tidx = length(tt(tt<=TimeStamps(i)))+1;
    end
    L(tidx:end) = isequal(1,TTLs(end))*ones(size(L(tidx:end)));
    if TTLs(end) ~= 0
        plot(tt(tidx:end),x(tidx:end),'k')
        %scatter(tt(tidx:end),ch_X(tidx:end),'g')
    else
        plot(tt(tidx:end),x(tidx:end),'b')
        %scatter(tt(tidx:end),ch_X(tidx:end),'k')
    end
    xlabel('time (sec)');ylabel('potential (uV)');xlim([min(tt) max(tt)])
    
    shnum = newshnum;
end

%plot off/on transitions for each shank
FS = [1 0 1 0 1]; EH = 0; EM = 1;
[TimeStamps, TTLs, EventStrings] = Nlx2MatEV(efile,FS,EH,EM);
TTLs(1) = []; TimeStamps(1) = []; TimeStamps = TimeStamps/1000000;
figure
shnum = 1; trspace = 4000; %trace spacing in uV
amin = 0.15;cmax = 0.8;a=0.5;
for jj = 1:numchannels
    % Load data from the .ncs files, make plots, and store them
    newshnum = str2num(channels{jj}(2));
    subplot(1,4,shnum)
    title(strcat('shank ',num2str(shnum)))
    elnum = str2num(channels{jj}(6:end-4));
    xfile = [cd,'\',channels{jj}];
    [samples,~,tt, Fs, bv, ~] = loadEEG2(xfile);
    x = bv*samples; clear samples
    x = x(1:16:end)-elnum*trspace; %downsample traces to 2000 Hz
    tt = tt(1:16:end);
    mt = min(tt); tt = tt-mt;
    if jj == 1
        TimeStamps = TimeStamps - mt;
    end
    
    %make light vector
    L  = zeros(size(tt));
    hold on
    trans = find(TTLs==0);
    for i = 1:length(trans)
        if trans(i)==1
            offidx = find(tt<TimeStamps(trans(i)),1,'last');
            startidx = find(tt>=TimeStamps(trans(i)),1,'first');
            endidx = find(tt<TimeStamps(trans(i)+1),1,'last');
            plot(tt(1:offidx)-TimeStamps(trans(i)),x(1:offidx),'color',[0 0 0 a])
            plot(tt(startidx:endidx)-TimeStamps(trans(i)),x(startidx:endidx),'color',[0 0 1 a])
        elseif trans(i)==length(TTLs)
            offidx = find(tt<TimeStamps(trans(i)) & tt>=TimeStamps(trans(i)-1));
            onidx = find(tt>=TimeStamps(trans(i)));
            plot(tt(offidx)-TimeStamps(trans(i)),x(offidx),'color',[0 0 0 a])
            plot(tt(onidx)-TimeStamps(trans(i)),x(onidx),'color',[0 0 1 a])
        else
            offidx = find(tt<TimeStamps(trans(i)) & tt>=TimeStamps(trans(i)-1));
            onidx = find(tt>=TimeStamps(trans(i)) & tt<TimeStamps(trans(i)+1));
            plot(tt(offidx)-TimeStamps(trans(i)),x(offidx),'color',[0 0 0 a])
            plot(tt(onidx)-TimeStamps(trans(i)),x(onidx),'color',[0 0 1 a])
        end                
    end
    xlabel('time (sec)');ylabel('potential (uV)');xlim([-1 3]);ylim([-17*trspace 0])
%     set(gca,'Color',[1 0.5 0.5]);
    shnum = newshnum;
end

%plot on/off transitions for each shank
FS = [1 0 1 0 1]; EH = 0; EM = 1;
[TimeStamps, TTLs, EventStrings] = Nlx2MatEV(efile,FS,EH,EM);
TTLs(1) = []; TimeStamps(1) = []; TimeStamps = TimeStamps/1000000;
figure
shnum = 1; trspace = 4000; %trace spacing in uV
amin = 0.15;cmax = 0.8;a=0.5;
for jj = 1:numchannels
    % Load data from the .ncs files, make plots, and store them
    newshnum = str2num(channels{jj}(2));
    subplot(1,4,shnum)
    title(strcat('shank ',num2str(shnum)))
    elnum = str2num(channels{jj}(6:end-4));
    xfile = [cd,'\',channels{jj}];
    [samples,~,tt, Fs, bv, ~] = loadEEG2(xfile);
    x = bv*samples; clear samples
    x = x(1:16:end)-elnum*trspace; %downsample traces to 2000 Hz
    tt = tt(1:16:end);
    mt = min(tt); tt = tt-mt;
    if jj == 1
        TimeStamps = TimeStamps - mt;
    end
    
    %make light vector
    L  = zeros(size(tt));
    hold on
    trans = find(TTLs~=0); %~0 switches it to the on/off transitions
    for i = 1:length(trans)
        if trans(i)==1
            offidx = find(tt<TimeStamps(trans(i)),1,'last');
            startidx = find(tt>=TimeStamps(trans(i)),1,'first');
            endidx = find(tt<TimeStamps(trans(i)+1),1,'last');
            plot(tt(1:offidx)-TimeStamps(trans(i)),x(1:offidx),'color',[0 0 1 a])
            plot(tt(startidx:endidx)-TimeStamps(trans(i)),x(startidx:endidx),'color',[0 0 0 a])
        elseif trans(i)==length(TTLs)
            offidx = find(tt<TimeStamps(trans(i)) & tt>=TimeStamps(trans(i)-1));
            onidx = find(tt>=TimeStamps(trans(i)));
            plot(tt(offidx)-TimeStamps(trans(i)),x(offidx),'color',[0 0 1 a])
            plot(tt(onidx)-TimeStamps(trans(i)),x(onidx),'color',[0 0 0 a])
        else
            offidx = find(tt<TimeStamps(trans(i)) & tt>=TimeStamps(trans(i)-1));
            onidx = find(tt>=TimeStamps(trans(i)) & tt<TimeStamps(trans(i)+1));
            plot(tt(offidx)-TimeStamps(trans(i)),x(offidx),'color',[0 0 1 a])
            plot(tt(onidx)-TimeStamps(trans(i)),x(onidx),'color',[0 0 0 a])
        end                
    end
    xlabel('time (sec)');ylabel('potential (uV)');xlim([-1 3]);ylim([-17*trspace 0])
%     set(gca,'Color',[1 0.5 0.5]);
    shnum = newshnum;
end

%plot full session for each tetrode
trspace = 1500; %trace spacing in uV
shnum = 1;
FS = [1 0 1 0 1]; EH = 0; EM = 1;
[TimeStamps, TTLs, EventStrings] = Nlx2MatEV(efile,FS,EH,EM);
TTLs(1) = []; TimeStamps(1) = []; TimeStamps = TimeStamps/1000000;
for jj = 1:numchannels
    % Load data from the .ncs files, make plots, and store them
    xfile = [cd,'\',channels{jj}];
    [samples,~,tt, Fs, bv, ~] = loadEEG2(xfile);
    x = bv*samples; clear samples
    x = x-jj*trspace; %downsample traces to 2000 Hz
    mt = min(tt); tt = tt-mt;
    if jj == 1
        TimeStamps = TimeStamps - mt;
    end
    
    %make light vector
    L  = zeros(size(tt)); tidx = 1;
    hold on
    for i = 1:length(TTLs)
        L(tidx:length(tt(tt<=TimeStamps(i)))) = TTLs(i);
        if TTLs(i) ~= 0
            plot(tt(tidx:length(tt(tt<=TimeStamps(i)))),x(tidx:length(tt(tt<=TimeStamps(i)))),'g')
            %scatter(tt(tidx:length(tt(tt<=TimeStamps(i)))),ch_X(tidx:length(tt(tt<=TimeStamps(i)))),'g')
        else
            plot(tt(tidx:length(tt(tt<=TimeStamps(i)))),x(tidx:length(tt(tt<=TimeStamps(i)))),'k')
            %scatter(tt(tidx:length(tt(tt<=TimeStamps(i)))),ch_X(tidx:length(tt(tt<=TimeStamps(i)))),'k')
        end
        tidx = length(tt(tt<=TimeStamps(i)))+1;
    end
    L(tidx:end) = isequal(1,TTLs(end))*ones(size(L(tidx:end)));
    if TTLs(end) ~= 0
        plot(tt(tidx:end),x(tidx:end),'k')
        %scatter(tt(tidx:end),ch_X(tidx:end),'g')
    else
        plot(tt(tidx:end),x(tidx:end),'g')
        %scatter(tt(tidx:end),ch_X(tidx:end),'k')
    end
    xlabel('time (sec)');ylabel('potential (uV)');xlim([min(tt) max(tt)])
end

%plot off/on transitions for tetrodes
FS = [1 0 1 0 1]; EH = 0; EM = 1;
[TimeStamps, TTLs, EventStrings] = Nlx2MatEV(efile,FS,EH,EM);
TTLs(1) = []; TimeStamps(1) = []; TimeStamps = TimeStamps/1000000;
figure
trspace = 1500; %trace spacing in uV
amin = 0.15;cmax = 0.8;a=0.5;
for jj = 1:numchannels
    % Load data from the .ncs files, make plots, and store them
    elnum = str2num(channels{jj}(4:end-4));
    xfile = [cd,'\',channels{jj}];
    [samples,~,tt, Fs, bv, ~] = loadEEG2(xfile);
    x = bv*samples; clear samples
    x = x-jj*trspace; %downsample traces to 2000 Hz
    mt = min(tt); tt = tt-mt;
    if jj == 1
        TimeStamps = TimeStamps - mt;
    end
    
    %make light vector
    L  = zeros(size(tt));
    hold on
    trans = find(TTLs==0);
    for i = 1:length(trans)
        if trans(i)==1
            offidx = find(tt<TimeStamps(trans(i)),1,'last');
            startidx = find(tt>=TimeStamps(trans(i)),1,'first');
            endidx = find(tt<TimeStamps(trans(i)+1),1,'last');
            plot(tt(1:offidx)-TimeStamps(trans(i)),x(1:offidx),'color',[0 0 0 a])
            plot(tt(startidx:endidx)-TimeStamps(trans(i)),x(startidx:endidx),'color',[0 0 1 a])
        elseif trans(i)==length(TTLs)
            offidx = find(tt<TimeStamps(trans(i)) & tt>=TimeStamps(trans(i)-1));
            onidx = find(tt>=TimeStamps(trans(i)));
            plot(tt(offidx)-TimeStamps(trans(i)),x(offidx),'color',[0 0 0 a])
            plot(tt(onidx)-TimeStamps(trans(i)),x(onidx),'color',[0 0 1 a])
        else
            offidx = find(tt<TimeStamps(trans(i)) & tt>=TimeStamps(trans(i)-1));
            onidx = find(tt>=TimeStamps(trans(i)) & tt<TimeStamps(trans(i)+1));
            plot(tt(offidx)-TimeStamps(trans(i)),x(offidx),'color',[0 0 0 a])
            plot(tt(onidx)-TimeStamps(trans(i)),x(onidx),'color',[0 0 1 a])
        end                
    end
    xlabel('time (sec)');ylabel('potential (uV)');xlim([-1 2]);ylim([-(jj+1)*trspace 0])
end
axis square

%plot on/off transitions for tetrodes
FS = [1 0 1 0 1]; EH = 0; EM = 1;
[TimeStamps, TTLs, EventStrings] = Nlx2MatEV(efile,FS,EH,EM);
TTLs(1) = []; TimeStamps(1) = []; TimeStamps = TimeStamps/1000000;
figure
trspace = 1500; %trace spacing in uV
amin = 0.15;cmax = 0.8;a=0.5;
for jj = 1:numchannels
    % Load data from the .ncs files, make plots, and store them
    elnum = str2num(channels{jj}(4:end-4));
    xfile = [cd,'\',channels{jj}];
    [samples,~,tt, Fs, bv, ~] = loadEEG2(xfile);
    x = bv*samples; clear samples
    x = x-jj*trspace; %downsample traces to 2000 Hz
    mt = min(tt); tt = tt-mt;
    if jj == 1
        TimeStamps = TimeStamps - mt;
    end
    
    %make light vector
    L  = zeros(size(tt));
    hold on
    trans = find(TTLs~=0);
    for i = 1:length(trans)
        if trans(i)==1
            offidx = find(tt<TimeStamps(trans(i)),1,'last');
            startidx = find(tt>=TimeStamps(trans(i)),1,'first');
            endidx = find(tt<TimeStamps(trans(i)+1),1,'last');
            plot(tt(1:offidx)-TimeStamps(trans(i)),x(1:offidx),'color',[0 0 1 a])
            plot(tt(startidx:endidx)-TimeStamps(trans(i)),x(startidx:endidx),'color',[0 0 0 a])
        elseif trans(i)==length(TTLs)
            offidx = find(tt<TimeStamps(trans(i)) & tt>=TimeStamps(trans(i)-1));
            onidx = find(tt>=TimeStamps(trans(i)));
            plot(tt(offidx)-TimeStamps(trans(i)),x(offidx),'color',[0 0 1 a])
            plot(tt(onidx)-TimeStamps(trans(i)),x(onidx),'color',[0 0 0 a])
        else
            offidx = find(tt<TimeStamps(trans(i)) & tt>=TimeStamps(trans(i)-1));
            onidx = find(tt>=TimeStamps(trans(i)) & tt<TimeStamps(trans(i)+1));
            plot(tt(offidx)-TimeStamps(trans(i)),x(offidx),'color',[0 0 1 a])
            plot(tt(onidx)-TimeStamps(trans(i)),x(onidx),'color',[0 0 0 a])
        end                
    end
    xlabel('time (sec)');ylabel('potential (uV)');xlim([-1 2]);ylim([-(jj+1)*trspace 0])
end
axis square