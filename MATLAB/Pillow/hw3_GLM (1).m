%% #5 - Generalized Linear Models - 50 Presentations for Stim 1
%     close all
%     clear all

%     load ('Stim1.mat')
%     load ('Stim2.mat')
%     load ('GLMparams.mat')

    % Instantaneous Firing Rate at time t = exp(k * x(t) + h * y(t) +b)
    % given: k, h, b
    % x = stimulus
    % y is the spike history at time t, [y(t-1),y(t-2),...,y(t-nh)]'

    y = [];
    yhist=[];
    figure
    hold on
    
    for j=1:50 
        for i=1:10  % We can't use a history filter if there is no history yet

            frate(i) = exp(k' * Stim1(i:i+39) + b);
            y(i) = rand(1)<=frate(i);
        end
        for i=11:(length(Stim1)-40) %Now we can incorporate a history for these trials

            yhist = y(i-1:-1:i-10);
            frate(i) = exp(k' * Stim1(i:i+39) + (yhist*h) + b);

            y(i) = rand(1)<=frate(i);    
        end
        spikes(j,:)=y;
        stimes=find(y); % When in y is there a 1 (spike)
         
        plot(stimes,j*ones(size(stimes)),'o')
    end    
    
    title('Raster for GLM Model of Stim1')
    xlabel('Spike Times')
    ylabel('Trial Number')
    
    % 
%% #6 - Generalized Linear Models- 200 presentations for Stim 2

% Stim2 Amplitude = 1
%     close all
%     clear all

%     load ('Stim1.mat')
%     load ('Stim2.mat')
%     load ('GLMparams.mat')
    
    y = [];
    yhist=[];
%     figure
%     hold on
    
    for j=1:200 
        for i=1:10

            frate(i) = exp(k' * Stim2(i:i+39) + b);
            y(i) = rand(1)<=frate(i);
        end
        for i=11:(length(Stim2)-40)

            yhist = y(i-1:-1:i-10);
            frate(i) = exp(k' * Stim2(i:i+39) + (yhist*h) + b);

            y(i) = rand(1)<=frate(i);    
        end
        spikes2(j,:)=y;
        stimes2=find(y);
         
%         plot(stimes2,j*ones(size(stimes2)),'o')
    end    
    
%     title('Raster for GLM Model of Stim2')
%     xlabel('Spike Times')
%     ylabel('Trial Number')

 	spkcnt1 = sum(spikes2,2);
    meanspkcnt1 = mean(spkcnt1);
    varspkcnt1 = var(spkcnt1);
    
 % Stim2 Amplitude = 0.125 
 
    Stim125 = .125 * Stim2;
    y = [];
    yhist=[];

    
    for j=1:200 
        for i=1:10

            frate(i) = exp(k' * Stim125(i:i+39) + b);
            y(i) = rand(1)<=frate(i);
        end
        for i=11:(length(Stim125)-40)

            yhist = y(i-1:-1:i-10);
            frate(i) = exp(k' * Stim125(i:i+39) + (yhist*h) + b);

            y(i) = rand(1)<=frate(i);    
        end
        spikes125(j,:)=y;
        stimes125=find(y);
         
    end 
    
	spkcnt125 = sum(spikes125,2);
    meanspkcnt125 = mean(spkcnt125);
    varspkcnt125 = var(spkcnt125);
 
% Stim2 Amplitude = 0.25 
 
    Stim25 = .25 * Stim2;
    y = [];
    yhist=[];

    
    for j=1:200 
        for i=1:10

            frate(i) = exp(k' * Stim25(i:i+39) + b);
            y(i) = rand(1)<=frate(i);
        end
        for i=11:(length(Stim25)-40)

            yhist = y(i-1:-1:i-10);
            frate(i) = exp(k' * Stim25(i:i+39) + (yhist*h) + b);

            y(i) = rand(1)<=frate(i);    
        end
        spikes25(j,:)=y;
        stimes25=find(y);
         
    end 
    
	spkcnt25 = sum(spikes25,2);
    meanspkcnt25 = mean(spkcnt25);
    varspkcnt25 = var(spkcnt25);
    
 % Stim2 Amplitude = 0.5 
 
    Stim5 = .5 * Stim2;
    y = [];
    yhist=[];

    
    for j=1:200 
        for i=1:10

            frate(i) = exp(k' * Stim5(i:i+39) + b);
            y(i) = rand(1)<=frate(i);
        end
        for i=11:(length(Stim5)-40)

            yhist = y(i-1:-1:i-10);
            frate(i) = exp(k' * Stim5(i:i+39) + (yhist*h) + b);

            y(i) = rand(1)<=frate(i);    
        end
        spikes5(j,:)=y;
        stimes5=find(y);
         
    end 
    
	spkcnt5 = sum(spikes5,2);
    meanspkcnt5 = mean(spkcnt5);
    varspkcnt5 = var(spkcnt5);
    
 % Stim2 Amplitude = 2
 
    Stim02 = 2 * Stim2;
    y = [];
    yhist=[];

    
    for j=1:200 
        for i=1:10

            frate(i) = exp(k' * Stim02(i:i+39) + b);
            y(i) = rand(1)<=frate(i);
        end
        for i=11:(length(Stim02)-40)

            yhist = y(i-1:-1:i-10);
            frate(i) = exp(k' * Stim02(i:i+39) + (yhist*h) + b);

            y(i) = rand(1)<=frate(i);    
        end
        spikes2(j,:)=y;
        stimes2=find(y);
         
    end 
    
	spkcnt2 = sum(spikes2,2);
    meanspkcnt2 = mean(spkcnt2);
    varspkcnt2 = var(spkcnt2);
    
% Stim2 Amplitude = 4
 
    Stim4 = 4 * Stim2;
    y = [];
    yhist=[];

    
    for j=1:200 
        for i=1:10

            frate(i) = exp(k' * Stim4(i:i+39) + b);
            y(i) = rand(1)<=frate(i);
        end
        for i=11:(length(Stim4)-40)

            yhist = y(i-1:-1:i-10);
            frate(i) = exp(k' * Stim4(i:i+39) + (yhist*h) + b);

            y(i) = rand(1)<=frate(i);    
        end
        spikes4(j,:)=y;
        stimes4=find(y);
         
    end 
    
	spkcnt4 = sum(spikes4,2);
    meanspkcnt4 = mean(spkcnt4);
    varspkcnt4 = var(spkcnt4);
    
    meanspkcnts = [meanspkcnt125,meanspkcnt25,meanspkcnt5,meanspkcnt1,meanspkcnt2,meanspkcnt4];
    varspkcnts = [varspkcnt125,varspkcnt25,varspkcnt5,varspkcnt1,varspkcnt2,varspkcnt4];
    amplitudes = [.125,.25,.5,1,2,4];
    
    figure
    hold on
    subplot(2,1,1)
    plot(amplitudes,meanspkcnts)
    title('Mean of Spike Counts')
    xlabel('Amplitudes')
    ylabel('Mean of Spike Counts')
    subplot(2,1,2)
    plot(amplitudes,varspkcnts)
    ylabel('Variance of Spike Counts')
    xlabel('Amplitudes')
    title('Variance of Spike Counts')