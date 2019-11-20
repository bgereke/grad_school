function [ rpm_fn ] = match_sig_to_rpm_fn( sig_path )
% give full path to a signal file, this will return the matching rpm file
%
% e.g.
% sig_path='X:\imaging\wheel_run_5\signals\md229_02262017_r2_003_signals.mat';
% match_sig_to_rpm_fn( sig_path );
% returns correct rpm file
% return blank string if two or more rpm files are found, or no rpm files





% make matching string ('animal_day')
[fp,sig_fn]=fileparts(sig_path);
ss=strsplit(sig_fn,'_');
animal=ss{1};
day=ss{2};
matching_ss=strcat(animal,'_',day);

% get list of rpm files
sf=fp; % sig folder
ind=strfind(sf,'signals');
sf(ind:end)='';
rf=strcat(sf,'rpm\');

files = dir(fullfile(rf, '*.txt')); 
fns={files.name}; % put fns in cell array
rpm_files=strcat(rf,fns); % add dir name 


% find match
index_cell = strfind(rpm_files,matching_ss);
index_match = find(not(cellfun('isempty', index_cell)));

% report
if (length(index_match) > 1)
    display('Multiple matches found, check signals and rpm files');
    rpm_fn='';
elseif ~isempty(index_match)
    rpm_fn=rpm_files{index_match};
else
    rpm_fn='';
    display('No Rpm file found for');
end
end

