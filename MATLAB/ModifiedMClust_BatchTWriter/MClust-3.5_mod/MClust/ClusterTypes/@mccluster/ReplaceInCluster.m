function [MCC] = ReplaceInCluster(f, MCC)

% this function is used in FindSpikesWithoutRP to replace the points in 
% mccluster
% Added by Ernie Hwaun 5-13-2015

MCC.myOrigPoints = find(f);