%mean vector length distributions across days

mice = string(['md229';'md230';'md231']);
dates = string(['02262017';'02272017';'02282017';'03012017';...
                '03022017';'03032017';'03042017']);
imdir = 'Z:\imaging\wheel_run_5\signals\';
rundir = 'Z:\imaging\wheel_run_5\rpm\';

nummice = size(mice,1);
numdays = size(dates,1);
numbins = 25;
hgrid = linspace(-1,1,numbins);
pgrid = [0.0517,0.1509,0.2501,0.3492,0.4484,0.5476,0.6467,0.7459,0.8451,0.9443];
MVL_d = zeros(numdays,numbins);
MVLp_d = zeros(numdays,length(pgrid));
frac_pos = zeros(numdays,3); %fraction of positive dbmvl for each day with bootstrapped confidence intervals
frac_sig = zeros(numdays,3); %fraction of significant mvl's for each day with bootstrapped confidence intervals
numboots = 10000; %numbootstraps

for d = 1:numdays
    DBMVL = cell(1,1); %store all the dbmvl's for bootstrapping 
    MVLP = cell(1,1); %store all the mvlp's for bootstrapping 
    mcount = 0;
    for m = 1:nummice
        filestart = strcat(mice(m,:),'_',dates(d,:));
        %check for existence and load detected transients
        exists = 0;
        cd(imdir)
        imfiles = dir('*.mat'); 
        for f = 1:numel(imfiles)
            if ~isempty(regexp(imfiles(f).name,filestart))
               load(imfiles(f).name,'-mat','nT');%add variables to load
               exists = 1;
            end
        end
        if exists
            mcount = mcount + 1;
            %load run data
            cd(rundir)
            runfiles = dir('*.txt');
            for f = 1:numel(runfiles)
                if ~isempty(regexp(runfiles(f).name,filestart))
                    runtable = readtable(runfiles(f).name,'Delimiter',',','ReadVariableNames',true);
                end
            end
            %get mvl for day mouse combo
            cdt=1/10.088781275221955;
            ct = cdt*ones(size(nt,1),1);
            ct = cumsum(ct);
            ard_timestamp = (runtable.ard_timestamp - min(runtable.ard_timestamp))/1000 + cdt;
            wheel_position = runtable.lap_position/64;
            wheel_position = wheel_position - floor(wheel_position);keyboard
%             [dF,nF,zF] = deltaF(raw,2);
            [maps,pmap,grid,mvl,dbmvl,mvlp,bpsl,bptl] = ratemap(nT,runtable.lap_position,ct,ard_timestamp,100,'vonMises',30,false);
            MVL_d(d,:) = MVL_d(d,:) + hist(dbmvl,hgrid);
            MVLp_d(d,:) = MVLp_d(d,:) + hist(mvlp,pgrid);
            DBMVL{mcount} = dbmvl;
            MVLP{mcount} = mvlp;
        end          
    end
    %do multilevel bootstrap
        nummice_d = numel(DBMVL);
        numpos = 0;numsig=0;numtot = 0;
        for nm = 1:nummice_d
            numpos = numpos + sum(DBMVL{nm}>0);
            numsig = numsig + sum(MVLP{nm}<0.05);
            numtot = numtot + length(DBMVL{nm});
        end
        frac_pos(d,1) = numpos/numtot;
        frac_sig(d,1) = numsig/numtot;
        bootfracpos = zeros(numboots,1);
        bootfracsig = zeros(numboots,1);
        for b = 1:numboots
            dbsamp = [];psamp = [];
            mboot = randsample(nummice_d,nummice_d,true); %random draw from this day's mice
            %random draw from mouse's dbmvl's
            for bm = 1:nummice_d
               dbsamp = [dbsamp; randsample(DBMVL{mboot(bm)},length(DBMVL{mboot(bm)}),true)]; 
               psamp = [psamp; randsample(MVLP{mboot(bm)},length(MVLP{mboot(bm)}),true)];
            end
            bootfracpos(b) = sum(dbsamp>0)/length(dbsamp);
            bootfracsig(b) = sum(psamp<0.05)/length(psamp);
        end
        frac_pos(d,2:3) = quantile(bootfracpos,[0.025 0.975]);
        frac_sig(d,2:3) = quantile(bootfracsig,[0.025 0.975]);
    disp(sprintf('%s%s','Completed ',num2str(d),' of ',num2str(numdays)))
end
%covert to probabilities
MVL_d = MVL_d./repmat(sum(MVL_d,2),1,size(MVL_d,2));
%make plot
c = colormap(copper(numdays));
for d = 1:numdays
 hold on
 plot(hgrid,MVL_d(d,:),'LineWidth',2,'Color',c(d,:)');
end
legend(strread(num2str(1:numdays),'%s'))
plot([0 0],[0 max(MVL_d(:))+0.01],'--k')
xlabel('mean vector length')
ylabel('probability')
xlim([-1 1])

figure
for d = 1:numdays
   hold on
   plot([d d],frac_pos(d,2:3),'-k');
   plot(d,frac_pos(d,1),'ok');
end
plot([0 numdays+1],[0.5 0.5],'--k');
xlabel('day')
ylabel('fraction of positive debiased mean vector lengths')
xlim([0 numdays+1])
ylim([0 1])

figure
for d = 1:numdays
   hold on
   plot([d d],frac_sig(d,2:3),'-k');
   plot(d,frac_sig(d,1),'ok');
end
plot([0 numdays+1],[0.05 0.05],'--k');
xlabel('day')
ylabel('fraction of mean vector lengths with p < 0.05')
xlim([0 numdays+1])
% ylim([0 1])