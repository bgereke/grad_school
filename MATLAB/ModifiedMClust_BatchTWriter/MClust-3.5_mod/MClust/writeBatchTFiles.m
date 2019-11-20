function [ttfile] = writeBatchTFiles(dirlist, clusterlist, ttlist)

%% --------writeBatchTFiles-------
% Batch-apply cluster files to .ntt files when given a
% generalized list of .ntt files, convex hull .cluster files, and a list of directories
% in which to find them.
%% -------------------------------

global MClust_Directory;
global MClust_NeuralLoadingFunction;
global MClust_Directory     % directory where the MClust main function and features reside
global MClust_FeatureData  % feature data
global MClust_TTfn         % file name for tt file
global MClust_NeuralLoadingFunction % Loading Engine 
global MClust_TTdn;    % directory name for tt file
global MClust_FDdn;
global MClust_Clusters

global MCLUST_VERSION



% check that the dirlist and ttlist inputs are valid

for aa=1:length(dirlist)
	% Move to each directory...
	eval(['cd ',dirlist{aa}]);

	ttfile=[];	
	for bb=1:length(ttlist)
		% load each .ntt file
		MClustResetGlobals('initialize');
		MClust_NeuralLoadingFunction = 'LoadTT_NeuralynxNT';
		MClust_TTdn = pwd;      % directory name for tt file
		MClust_FDdn =pwd;  
		MClust_TTfn = [ttlist{bb}];
		if(~exist([ttlist{bb},'.ntt'], 'file'))
			continue;
		end
		[t] = MClust_LoadNeuralData([ttlist{bb},'.ntt']);
		if(~isempty(t))
			CalculateFeatures(ttlist{bb}, {'Energy','Peak','PeakValleyDiff','Valley','WavePC2','WavePC3','area','waveFFT','wavePC1'});
			
			% apply the appropriate convex file
			ApplyConvexHulls([clusterlist{bb}])

			% write the .t files
			[~,ttfile_in]=WriteTFiles_mod1([]);
			if isempty(ttfile)
				ttfile=ttfile_in;
			else
				ttfile=horzcat(ttfile,ttfile_in);
			end
		end
	end
end