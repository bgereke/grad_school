function [] = batchTWriter_general();
if ispc
	dir_delim='\';
else
	dir_delim = '/';
end

%Get the directories

dirlist = dir;
nn=0;
for aa=1:length(dirlist)
	if(dirlist(aa).isdir) 
		temp=dirlist(aa).name;

		if(length(temp) < 5)
			continue;
		end
		if( strcmp(temp(1:5), 'sleep') || strcmp(temp(1:5),'begin'))   %
			nn=nn+1;
			goodList{nn} = [pwd,dir_delim,temp];
		end
	end
end
%ttList =[]; %{'TT1', 'TT11', 'TT12'};
temp=split_string(pwd,'\');
% expects cluster files to be stored at <drive>\<user>\<animal>\Clusters\<day> with only one file per-tetrode
clustDir = pwd;%[temp{1},dir_delim,temp{2},dir_delim,temp{3},dir_delim,'Clusters',dir_delim, temp{end}];

if ~isdir(clustDir)
	cd ..;
	return;
end

filelist = dir(clustDir);
nn=0;
clusterList=[];
for aa=1:length(filelist)
	temp=filelist(aa).name;
	if(~isempty(strfind(temp,'.clusters')))
		nn=nn+1;
		clusterList{nn} = [clustDir, dir_delim, temp];
		temp2 = split_string(temp,'.');
		ttList{nn} = temp2{1};
		
	end
end

%clusterList = {[clustDir,dir_delim,'TT1_20131019_cleaned.clusters']};
if(isempty(clusterList))
	cd ..;
	return;
end

%keyboard
writeBatchTFiles(goodList, clusterList, ttList)
cd ..;
fprintf('\n\n+++++++++++++++ T-Files Written! +++++++++++++++\n\n');