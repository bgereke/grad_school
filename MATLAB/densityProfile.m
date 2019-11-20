function [densityProfile] = densityProfile(data,subLen)
%% Initiate Variables

excZoneLen = round(subLen * 0.5);
dataLen = length(data);
proLen = dataLen - subLen + 1;
densityProfile = zeros(1,proLen);
if dataLen == size(data, 2)
    data = data';
end

%% Preprocess inputs to mass

[dataFreq, dataMu, dataSig] = massPre(data, dataLen, subLen);

%% Compute complexity vector

CV = zeros(proLen,1);  
for j = 1: proLen
    CV(j) = sqrt(sum(diff(data(j:j+subLen-1)).^2));
end

%% Main loop
completed = []; it = floor(proLen/2); delta = inf;
while delta > 1e-6 
    tmpProfile = densityProfile;
    for idx = 1:it:proLen    
        if sum(idx == completed) == 1
            continue
        end
        
        % compute the distance profile and update density profile
        query = data(idx:idx+subLen-1);
        distProfile = mass(dataFreq, query, dataLen, subLen, ...
            dataMu, dataSig, dataMu(idx), dataSig(idx));
        
        %complexity correction
        distProfile = sqrt(abs(distProfile)).*max([CV CV(idx)*ones(proLen,1)],[],2)./...
            min([CV CV(idx)*ones(proLen,1)],[],2);
        
        [~,locs] = findpeaks(-distProfile,'MinPeakDist',2*excZoneLen);
        
        excZoneStart = max(1, idx - excZoneLen);
        excZoneEnd = min(proLen, idx + excZoneLen);
        locs(locs>=excZoneStart & locs<=excZoneEnd) = [];
        dSort = sort(distProfile(locs));
        
        densityProfile(idx) = mean(dSort(1:round(0.01*length(dSort))));
        
        completed = [completed idx];
    end
    it = floor(it/2);
    densityProfile = pchip(completed,densityProfile(completed),1:proLen);
    delta = sqrt(sum((densityProfile-tmpProfile).^2))/sqrt(sum(densityProfile.^2));    
    fprintf('%s%s\n','Completed: ',num2str(length(completed)/proLen*100));
    fprintf('%s%s\n','Delta: ',num2str(delta));keyboard
end

%%
% The following two functions are modified from the code provided in the 
% following URL
% http://www.cs.unm.edu/~mueen/FastestSimilaritySearch.html
function [dataFreq, dataMu, dataSig] = massPre(data, dataLen, subLen)
data(dataLen + 1:(subLen + dataLen)) = 0;
dataFreq = fft(data);
dataCumsum = cumsum(data);
data2Cumsum =  cumsum(data .^ 2);
data2Sum = data2Cumsum(subLen:dataLen) - ...
    [0; data2Cumsum(1:dataLen - subLen)];
dataSum = dataCumsum(subLen:dataLen) - ...
    [0; dataCumsum(1:dataLen - subLen)];
dataMu = dataSum ./ subLen;
data2Sig = (data2Sum ./ subLen) - (dataMu .^ 2);
dataSig = sqrt(data2Sig);

%%
function distProfile = mass(dataFreq, query, ...
    dataLen, subLen, dataMu, dataSig, queryMu, querySig)
query = query(end:-1:1);
query(subLen+1:(subLen+dataLen)) = 0;
queryFreq = fft(query);
productFreq = dataFreq .* queryFreq;
product = ifft(productFreq);
distProfile = 2 * (subLen - ...
    (product(subLen:dataLen) - subLen * dataMu * queryMu) ./ ...
    (dataSig * querySig));