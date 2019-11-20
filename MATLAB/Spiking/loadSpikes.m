function [TS] = loadSpikes(tFile)

tfp = fopen(tFile,'rb','b');
if tfp == -1
    keyboard
    disp('ERROR');
    disp('  Couldn''t load the t-file, please check the file name and path');
    return;
end

% Read header if present
ReadHeader(tfp);
S = cell(1,1);
% Load data into a cell array of length 1
S{1} = fread(tfp,inf,'uint32');
S{1} = ts(S{1});
fclose(tfp);
% Convert the timestamps to seconds
TS = Data(S{1}) / 10000;



function H = ReadHeader(fp)
% H = ReadHeader(fp)
%  Reads NSMA header, leaves file-read-location at end of header
%  INPUT: 

%      fid -- file-pointer (i.e. not filename)
%  OUTPUT: 
%      H -- cell array.  Each entry is one line from the NSMA header
% Now works for files with no header.
% ADR 1997
% version L4.1
% status: PROMOTED
% v4.1 17 nov 98 now works for files sans header
%---------------

% Get keys
beginheader = '%%BEGINHEADER';
endheader = '%%ENDHEADER';

iH = 1; H = {};
curfpos = ftell(fp);

% look for beginheader
headerLine = fgetl(fp);
if strcmp(headerLine, beginheader)
   H{1} = headerLine;
   while ~feof(fp) & ~strcmp(headerLine, endheader)     
      headerLine = fgetl(fp);
      iH = iH+1;
      H{iH} = headerLine;
   end
else % no header
   fseek(fp, curfpos, 'bof');
end