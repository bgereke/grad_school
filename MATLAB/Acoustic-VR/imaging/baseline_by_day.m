%baseline correction and transient detection across days

mice = string(['md229';'md230';'md231']);
dates = string(['02262017';'02272017';'02282017';'03012017';...
                '03022017';'03032017';'03042017']);
imdir = 'Z:\imaging\wheel_run_5\signals\';
rundir = 'Z:\imaging\wheel_run_5\rpm\';

nummice = size(mice,1);
numdays = size(dates,1);

for d = 1:numdays
    mcount = 0;
    for m = 1:nummice
        data = [];
        filestart = strcat(mice(m,:),'_',dates(d,:));
        %check for existence and load calcium data
        exists = 0;
        cd(imdir)
        imfiles = dir('*signals.mat');
        for f = 1:numel(imfiles)
            if ~isempty(regexp(imfiles(f).name,filestart))
               %load raw data
               load(imfiles(f).name)
               raw = double(squeeze(raw)');
               %do baseline correction and transient detection
               [dF,nF,zF] = deltaF(raw,4);
               [dT,nT,zT] = detect_trans(dF,nF,zF);
               %save all data to .mat struct
               data.dF = dF;data.nF = nF;data.zF = zF;
               data.dT = dT;data.nT = nT;data.zT = zT;               
               filename = char(strcat(filestart,'.mat'));
               save(filename,'-struct','data');
            end
        end                
    end
    disp(sprintf('%s%s','Completed ',num2str(d),' of ',num2str(numdays)))
end
