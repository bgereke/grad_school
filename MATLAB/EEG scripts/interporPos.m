function [x,y] = interporPos(x,y,timeTreshold,sampRate)

% Turn of warnings
warning off;

% Number of samples that corresponds to the time treshold.
sampTreshold = floor(timeTreshold * sampRate);

% number of samples
numSamp = length(x);
% Find the indexes to the missing samples
temp1 = 1./x;
temp2 = 1./y;
indt1 = isinf(temp1);
indt2 = isinf(temp2);
ind = indt1 .* indt2;
ind2 = find(ind==1);
% Number of missing samples
N = length(ind2);

if N == 0
    % No samples missing, and we return
    return
end

change = 0;

% Remove NaN in the start of the path
if ind2(1) == 1
    change = 1;
    count = 0;
    while 1
        count = count + 1;
        if ind(count)==0
            break
        end
    end
    x(1:count) = x(count);
    y(1:count) = y(count);
%     ind(1:count) = 0;
%     ind2(1:count) = [];
end

% Remove NaN in the end of the path
if ind2(end) == numSamp
    change = 1;
    count = length(x);
    while 1
       count = count - 1;
       if ind(count)==0
           break
       end
    end
    x(count:numSamp) = x(count);
    y(count:numSamp) = y(count);
end

if change
    % Recalculate the missing samples
    temp1 = 1./x;
    temp2 = 1./y;
    indt1 = isinf(temp1);
    indt2 = isinf(temp2);
    % Missing samples are where both x and y are equal to zero
    ind = indt1 .* indt2;
    ind2 = find(ind==1);
    % Number of samples missing
    N = length(ind2);
end

for ii = 1:N
    % Start of missing segment (may consist of only one sample)
    start = ind2(ii);
    % Find the number of samples missing in a row
    count = 0;
    while 1
        count = count + 1;
        if ind(start+count)==0
            break
        end
    end
    % Index to the next good sample
    stop = start+count;
    if start == stop
        % Only one sample missing. Setting it to the last known good
        % sample
        x(start) = x(start-1);
        y(start) = y(start-1);
    else
        if count < sampTreshold
            % Last good position before lack of tracking
            x1 = x(start-1);
            y1 = y(start-1);
            % Next good position after lack of tracking
            x2 = x(stop);
            y2 = y(stop);
            % Calculate the interpolated positions
            X = interp1([1,2],[x1,x2],1:1/count:2);
            Y = interp1([1,2],[y1,y2],1:1/count:2);
            % Switch the lacking positions with the estimated positions
            x(start:stop) = X;
            y(start:stop) = Y;

            % Increment the counter (avoid estimating allready estimated
            % samples)
            ii = ii+count;
        else
            % To many samples missing in a row and they are left as missing
            ii = ii+count;
        end
    end
end