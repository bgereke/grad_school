%ACF by day

nm = 5;
nd = 5;
licks = [];
pos = [];
nl = 400;
acfs = zeros(nm*nd,nl);
nb = 30;
grid = linspace(1,104,nb);
pmaps = zeros(nm*nd,nb);
lmaps = zeros(nm*nd,nb);
count = 1;

df = []; %time, position, licks, mouse, day

for m = 1:nm   
    for d = 1:nd

        %get licks & acfs
        temp_l = DATA((m-1)*nd+d).behavior.lick;
        flidx = find(diff(temp_l)==1) + 1;
        temp_fl = zeros(size(temp_l));
        temp_fl(flidx) = 1;
        acfs(count,:) = acf(temp_l,nl);
        
        %get position and maps
        temp_p = DATA((m-1)*nd+d).behavior.position_in_lap;
        pmaps(count,:) = hist(temp_p,grid);
        lmaps(count,:) = hist(temp_p(flidx),grid);
        
        %get time
        temp_t = DATA((m-1)*nd+d).behavior.run.timestamp;
        temp_t = (temp_t - min(temp_t)) / 1000;
        
        %get speed
        speed = DATA((m-1)*nd+d).behavior.run.rpm_filt;
        
        %get rewards
        rew = DATA((m-1)*nd+d).behavior.reward;
        ridx = find(diff(rew)==1)+1;
        rew = zeros(size(rew));
        rew(ridx) = 1;
        
        %get time since reward
        rt = temp_t(ridx);
        time_rew = temp_t;
        for t=1:length(temp_t)
            if temp_t(t) <= min(rt)
                 temp_t(t) = nan0;
            else
                time_rew(t) = temp_t(t) - max(rt(rt<temp_t(t)));
            end
        end
        
        %create dataframe
        df = [df;temp_t temp_p temp_l temp_fl speed rew time_rew repmat(m,length(temp_l),1) repmat(d,length(temp_l),1)];
        
        count = count + 1;
        
    end    
end

names = cell(1,9);
names{1} = 'Time';
names{2} = 'Position';
names{3} = 'Licks';
names{4} = 'FirstLicks';
names{5} = 'Speed';
names{6} = 'Reward';
names{7} = 'RTime';
names{8} = 'Mouse';
names{9} = 'Day';

rtable = array2table(df,'VariableNames',names);
rfilename = sprintf('%s%s%s%s','rtable_Licking.csv');
writetable(rtable,rfilename);
