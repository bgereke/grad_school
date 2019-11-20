function [redraw, rekey, undoable] = FindSpikesWithoutRP(iClust)

% f = FindSpikesWithoutRP(iClust, varargin)
%
% INPUTS
%     iClust - index into MClust_Clusters
% REQUIRES
%  figHandle from parent
%
% OUTPUTS
%     MCC - The updated cluster
%
% 
% Added by Ernie Hwaun 5-13-2015
%

global MClust_Clusters
global MClust_Hide MClust_AvailableClusterTypes
global MClust_FeatureTimestamps

redraw= true; rekey = true; undoable = true;

[f MClust_Clusters{iClust}] = FindInCluster(MClust_Clusters{iClust});
timediff = diff(MClust_FeatureTimestamps(f));
WithoutRPIndex = find(timediff < 10);
SpikeWithoutRPIndex = sort([f(WithoutRPIndex); f(WithoutRPIndex+1)]);

ClustInd = repmat(0,size(MClust_FeatureTimestamps));
ClustInd(SpikeWithoutRPIndex) = 1;

%Create mccluster for displaying the spikes without RP
    clusterType = MClust_AvailableClusterTypes{get(findobj('Tag', 'AddAsType'), 'value')};
    figHandle = findobj('Type','figure','Tag', 'ClusterCutWindow');
    
    if isempty(MClust_Clusters)
        MClust_Clusters{1} = feval(clusterType, 'Cluster 01');
        MClust_Hide(2) = 0;
    else
        MClust_Clusters{end+1} = feval(clusterType, sprintf('Cluster %02d', length(MClust_Clusters)+1));
        MClust_Hide(length(MClust_Clusters) + 1) = 0;
    end
    MClustCutterClearClusterKeys(figHandle);
    MClustCutterRedrawClusterKeys(figHandle, max(0,length(MClust_Clusters)-16));


[newCluster] = ReplaceInCluster(ClustInd, MClust_Clusters{end});

newCluster = SetName(newCluster, [' WithoutRP ' GetName(newCluster)]);
MClust_Clusters{end} = newCluster;
MClust_Hide(end) = 0;
	
warndlg(['Cluster ' num2str(iClust) ' Spikes without RP are copied to cluster ' ...
	num2str(length(MClust_Clusters)) '.'], 'Copy successful');

[MDistance,FeaturesCombination] = ClusterSeparation_EH(length(MClust_Clusters),iClust);

    if size(MDistance,2) <= 10
        [~,sortInd] = sort(max(MDistance,[],2),'descend');
    elseif size(MDistance,2) > 10
        [~,sortInd] = sort(mean(MDistance,2),'descend');
    end
    
for rr = 1:15 %include the top 15 projection with the furtherest Mahalanobis distance
    rank{rr} = FeaturesCombination{sortInd(rr)};
end

msg{1} = GetName(MClust_Clusters{end});
msg{2} = ['class(' class(MClust_Clusters{end}) ')'];
msg{3} = '-----';
msg = cat(2,msg,rank);
msgbox(msg, ['Cluster ', num2str(length(MClust_Clusters))]);

