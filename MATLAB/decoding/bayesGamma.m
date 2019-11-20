function [scores] =  bayesGamma(scores,listfile)

f1 = 60;
f2 = 100;
s1 = 25;
s2 = 55;

tetnums = gettetnumbers(listfile);

%load the tetrodes of interest and get their TFRs
disp('Calculating power for the following Tetrodes:');
for i = 1:size(tetnums,1)
    disp(strcat('CSC',int2str(tetnums{i,1}),'.ncs'));
    temptet = strcat('CSC',int2str(tetnums{i,1}),'.ncs');
    [tetnums{i,2},~,tetnums{i,3}] = loadEeg8(temptet);
    tetnums{i,4} = TFR_frequency_band(tetnums{i,2},2000,5,f1,f2);
    tetnums{i,5} = TFR_frequency_band(tetnums{i,2},2000,5,s1,s2);
    tetnums{i,6} = zscore(tetnums{i,4});
    tetnums{i,7} = zscore(tetnums{i,5});
end

%go through each window listed and get the power for fast and slow
for i = 1:size(scores,1)
    startEegInd = spikeTStoEegIndex(scores{i,2}(1,1),tetnums{1,3},2000);
    stopEegInd =  spikeTStoEegIndex(scores{i,2}(1,2),tetnums{1,3},2000);
    %get mean slow/fast power for the window from each tet
    meanF = 0;meanFz = 0;
    meanS = 0;meanSz = 0;
    for j = 1:size(tetnums,1)
        meanF = meanF + mean(tetnums{j,4}(startEegInd:stopEegInd));
        meanS = meanS + mean(tetnums{j,5}(startEegInd:stopEegInd));
        meanFz = meanFz + mean(tetnums{j,6}(startEegInd:stopEegInd));
        meanSz = meanSz + mean(tetnums{j,7}(startEegInd:stopEegInd));
    end
    scores{i,23} = meanF/j;
    scores{i,24} = meanS/j;
    scores{i,25} = meanFz/j;
    scores{i,26} = meanSz/j;
end


%figure out which tetrodes to use
function tets = gettetnumbers(file)

tets = nan(numel(textread(file,'%1c%*[^\n]')),1); %#ok<REMFF1>

fid = fopen(file,'r');
if fid == -1
    msgbox('Could not open the input file','ERROR');
end

for i = 1:size(tets,1)
  tline = fgetl(fid);
  if ~ischar(tline) 
      break 
  else
      if tline(4)=='_'
        tets(i) = str2double(tline(3));
      elseif tline(5) == '_'
        tets(i) = nan;
      end
  end          
end
       
fclose(fid);

tets(isnan(tets)) = [];
tets(tets>6) = [];
tets = unique(tets);
tets = num2cell(tets);
