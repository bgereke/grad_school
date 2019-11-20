function Analyze_All(mice)

img_text = 'on';

fid = fopen(mice,'r');
if fid == -1
    msgbox('Could not open the input file! Make sure the filname and path are correct.','ERROR');
end

% Get the mice to be analyzed from input file
fid = fopen(mice,'r');
ii = 1;     
while ~feof(fid)
    str = fgetl(fid);
    if ~strcmp(str(end),'\')
        str = strcat(str,'\');
    end
    Mice(ii) = {str};
    ii = ii+1;
end
numMice = ii-1;   
cycles = [];
lags = 8000;
sgauto = zeros(lags+1,5);
fgauto = zeros(lags+1,5);
cross = zeros(2*lags+1,5);
count = 1;
spec = [];ONS = []; OFFS = [];
Cells = [];Multiu = [];Fields = [];
for ii = 1:numMice
%     Mice{ii}(1) = 'E';
    disp(sprintf('%s%s','Reading data for: ',Mice{ii}));
    
    % find and open the dates.txt file
    dates_txt = strcat(Mice{ii},'dates.txt');   
    fid = fopen(dates_txt,'r');
    jj = 1; 
    while ~feof(fid)
    	str = fgetl(fid);
        if ~strcmp(str(end),'\')
            str = strcat(str,'\');
        end
        Dates(jj) = {str};
        jj = jj+1;
    end
    numDates = jj-1;
    
%     Tot_Rates = cell(numDates,1);

    % Perform all analyses for each date
    for jj=1:numDates
        disp(Dates{jj});
        newpath = strcat(Mice{ii},Dates{jj});
        cd(newpath);
%         inFile = 'infile.txt';
        nvoice = 32;
        freq = 2.^(1:1/nvoice:7);freq = freq(1:end-1);
%         freq = 2.^linspace(log2(2),log2(100),50);
        inFile = 'InFile.txt';
%         inFile = 'infile - Copy.txt';
        
%         % Check if subdir for storing data is present. If not, it is
%         % created
%         dirInfo = dir(newpath);
%         found = 0;
%         
%         for ss=1:size(dirInfo,1)
%             if dirInfo(ss).isdir
%                 if strcmp(dirInfo(ss).name,strcat('Cell_Rates','\'))
%                     found = 1;
%                 end
%             end
%         end
%         
%         if found==0
%             mkdir(newpath,strcat('Cell_Rates','\'));
%         end
%         
%         [data.Rates] = Rates(inFile);
%         filename = sprintf('%s%s%s%s',newpath,strcat('Cell_Rates','\'),'.mat');
%         save(filename,'-struct','data');
%         
%         Tot_Rates{jj} = data.Rates;
        
        % All the analyses to be performed
%         equalPlot_max_original('inFile.txt',0)
%         powerspectrum(inFile);
%         STA(inFile);
%         STA_filt(inFile);
%         CFC(inFile,2:2:100);
%         CFC2ch(inFile,2:2:100);
%         Coherence(inFile,2:2:100);
%         WPLI(inFile,2:2:100,6);
%         WPLI_one(inFile,2:2:100,6);
%         crossWPLI(inFile,2:2:100,6,1,2);
%         bestCSCs_batch(inFile,2:2:100,7,5);
%         MCWPLI_batch(inFile,2:2:100,7,5);  
%         traces2CS_rats(inFile,freq,6);
    traces2CS(inFile,freq,6);
%     traces2CICS_rats(inFile,freq,6);
%     traces2CCS_biv_rats(inFile,freq,6);
%         ratTraces2CS(inFile,2:2:100,7);
%         mat2Rcsv(inFile,'CS','CS.mat')
%         GIMt_batch(inFile,20,300);
%         VFR_GIM_batch(inFile,2.5,300);
%         AFR_GIM_batch_rats(inFile,8,300);
%         VFR_GIM_batch_rats(inFile,5,300);
%         PFR_GIM_batch_rats(inFile,10,300);
%         WPLI_bias_comp(inFile,2:2:100,7,5);
%         traces2MCWPLI(inFile,2:2:100,7);
%         SpecCorr(inFile,2:2:100,6);
%         XFSpecCorr(inFile,2:2:100,6);
%         RPD_batch(inFile,2:2:100,50);  
%         VFR_batch_rats(inFile,exp(linspace(log(2),log(250),65)),7);
%           VFR_batch(inFile,freq,6);
%         [cycle]=PFR_batch_rats(inFile,exp(linspace(log(10),log(250),65)),linspace(-pi,pi,60));
%         cycles = [cycles; cycle];
%           GMD(inFile,exp(linspace(log(10),log(250),65)),linspace(-pi,pi,60));
%           [sgauto(:,count), fgauto(:,count), cross(:,count)] = gamma_xcorr(inFile,lags,23,25,50,52,58,60,100,102);
%           count = count+1;
%           readHeader(inFile);
%         PFR_batch_2ch(inFile,2:2:100,-pi:10/360*2*pi:pi);
%         PFR_WPLI_batch(inFile,20:2:100,50);
%         PxWPLI_batch(inFile,2:2:100,50);
%         PFR_WRFR_batch(inFile,2:2:100,-pi:10/360*2*pi:pi);
%         VFR_WPLI_batch(inFile,18:2:100);
%         xVFR_WPLI_batch(inFile,2:2:100);
%         VFR_WPLI_joint(inFile,2:2:100);
%         VFR_ImX_batch(inFile,2:2:100);
%         PFR_ImX_batch(inFile,2:2:100,-pi:10/360*2*pi:pi);
%         PFR_noPi_batch(inFile,2:2:100,-pi:10/360*2*pi:pi);
%         LFPmap_batch(inFile,2:2:100,0,3,5);
%         LFPmap_batch(inFile,20:2:40,0,10,6);
%         LFPmap_batch(inFile,50:2:90,0,10,6);   
%         WPLImap_batch(inFile,20:2:40,0,10,6);
%         WPLImap_batch(inFile,50:2:90,0,10,6);
%         ImXmap_batch(inFile,20:2:40,0,10,6);
%         ImXmap_batch(inFile,50:2:90,0,10,6);
%         CFC_filt(inFile,1:5:100);
%         trace_filt(inFile);
%         gammawindows(inFile,25,45,55,95);
%         Spike_Phase_Hist(inFile,20,40,50,90);
%         Spike_Phase_Hist_MultiU(inFile,20,40,50,90);
%         TFR(inFile,2:2:100,6);
%         TFR_w_light(inFile,2:2:100,6);
%         TFR_opto(inFile,2:2:100,6,0:30);
%         WPLI_opto_ts(inFile,2:2:100,6,0:30);
%         powerspectrum_opto(inFile,2:2:100,6,250);
%         [ons offs] = load_specs_opto(inFile);
%         ONS = [ONS ons];
%         OFFS = [OFFS offs];
%         MCWPLI_opto(inFile,250);
%         CIC_opto(inFile,250);
%         rates_opto(inFile,10,0:0.2:30,linspace(-pi,pi,100),0.5,100);
%         [fields] = fields_opto(inFile,250,0:0.2:30,linspace(-pi,pi,100),0.5,100);
%         Fields = [Fields; fields];
%         [cells, multiu] = multiu_opto(inFile,0.5:1:29.5,0.25);
%         Cells = [Cells;cells];
%         Multiu = [Multiu;multiu];
%         spike_angles(inFile,0);
%         wpli_opto(inFile,2:2:100,6,250);
%         theta_triggered_TFR(inFile,1:2:100,7,3);     
%         bayes_decode(inFile,2:2:150,-pi:10/360*2*pi:pi,7);
%         bayes_decode(inFile,2:2:150,-pi:10/360*2*pi:pi,7);
%         cycle_vels(inFile,2:2:100,-pi:10/360*2*pi:pi,7);
%         PLVspec(inFile,2:2:100,6);=
    end
    
%     % Check if subdir for storing data is present. If not, it is
%     % created
%     dirInfo = dir(Mice{ii});
%     found = 0;
% 
%     for ss=1:size(dirInfo,1)
%         if dirInfo(ss).isdir
%             if strcmp(dirInfo(ss).name,strcat('Tot_Cell_Rates','\'))
%                 found = 1;
%             end
%         end
%     end
% 
%     if found==0
%         mkdir(Mice{ii},strcat('Tot_Cell_Rates','\'));
%     end
%     
%     Data.Rates = Tot_Rates;
%     filename = sprintf('%s%s%s%s',Mice{ii},strcat('Tot_Cell_Rates','\'),'.mat');
%     save(filename,'-struct','Data');
end

keyboard



