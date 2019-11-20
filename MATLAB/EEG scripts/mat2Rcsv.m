% CFC('inFile.txt')
%
% CFC generates cross frequency coherences, plots them and stores them to images and .mat files.
% Multiple sessions will be read from the CSC's specififed in
% 'CSCList.txt'. 'TTList.txt' is necessary for spike data scripts that use
% the same 'infile.txt'.
%
% The input file must be on the following format.
%
% C:\Data\TTList.txt
% C:\Data\CSCList.txt
% C:\Data\Begin 1
% C:\Data\Begin 2
% C:\Data\Begin 3
% C:\Data\Begin 4
% and so on ...
%
% 'CSCList.txt' contains a list of the Neuralynx .csc files to be analyzed.
% All plots will be stored to both bmp and eps imagefiles to a subdirectory in
% the data folder called CFC_plots.

function mat2Rcsv(inFile,folder,filename)

img_text = 'on';

fid = fopen(inFile,'r');
if fid == -1
    msgbox('Could not open the input file! Make sure the filname and path are correct.','ERROR');
end

% Get sessions and csc-file list from input file
fid = fopen(inFile,'r');
ii = -1;     
numsessions = 0;
while ~feof(fid)
    str = fgetl(fid);
    if ii == 0
        cscList = str;
    elseif ii == 1
        refList = str;
        elseif ii == 2
        ch4avg = str;
    elseif ii > 2
        numsessions  = numsessions+1;
        if ~strcmp(str(end),'\')
            str = strcat(str,'\');
        end
        sessions(numsessions) = {str};        
    end
    ii = ii+1;
end

for ii = 1:numsessions
    disp(sprintf('%s%s','Reading data for session: ',sessions{ii}));
    data = load(strcat(sessions{ii},folder,'\',filename));
    names = cell(1,4*length(data.freqVec)+4);
    for f = 1:length(data.freqVec)
        names{f} = strcat('PICS_f',num2str(f));
        names{f+length(data.freqVec)} = strcat('CICS_f',num2str(f));
        names{f+2*length(data.freqVec)} = strcat('RePCP_f',num2str(f));
        names{f+3*length(data.freqVec)} = strcat('ImPCP_f',num2str(f));
    end
    names{4*length(data.freqVec)+1} = 'Time';
    names{4*length(data.freqVec)+2} = 'Position';
    names{4*length(data.freqVec)+3} = 'Running_Speed';
    names{4*length(data.freqVec)+4} = 'Acceleration';
    mt = min(data.time);
    data.time = data.time - mt;
    rtable = array2table([data.PA',imag(data.PCS').*repmat(sign(mean(imag(data.PCS'))),size(data.PCS,2),1),...
        real(data.phaseCSC'),imag(data.phaseCSC'),data.time,data.phase,data.vel,data.acc],'VariableNames',names);
    rfilename = sprintf('%s%s%s%s',strcat(sessions{ii},folder,'\'),'rtable_CSC.csv');
    writetable(rtable,rfilename)
    
    tsnames = cell(1,length(data.freqVec)+4);
    tsnames{1} = 'CellID';
    tsnames{2} = 'SpTime';
    tsnames{3} = 'SpPosition';
    tsnames{4} = 'SpRunningSpeed';
    for f = 1:length(data.freqVec)
        tsnames{4+f} = strcat('PCP_f',num2str(f));
    end
    data.TS(:,2) = data.TS(:,2) - mt;
    rtable = array2table([data.TS],'VariableNames',tsnames);
    rfilename = sprintf('%s%s%s%s',strcat(sessions{ii},folder,'\'),'rtable_TS.csv');
    writetable(rtable,rfilename)
end



