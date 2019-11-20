function [points segments] = getTextureSegments(texture)

shapeArray = texture.shapes;
for ndx = length(shapeArray):-1:1
    if strcmp(class(shapeArray{ndx}),'shapeColor')
        shapeArray(ndx) = [];
    end
end

% Get point coordinates
points = zeros(0,2);
segments = zeros(0,2);
for ndx = 1:length(shapeArray)
    [x y] = shapeArray{ndx}.coords2D;
    
    for loc = 1:length(shapeArray{ndx}.x)
        if isempty(find(x == shapeArray{ndx}.x(loc) & y == shapeArray{ndx}.y(loc),1))
            x = [x; NaN; shapeArray{ndx}.x(loc)]; %#ok<AGROW>
            y = [y; NaN; shapeArray{ndx}.y(loc)]; %#ok<AGROW>
        end
    end
    
    f = find(isnan(x) | isnan(y));
    f = [0; f; length(x)+1]; %#ok<AGROW>
    for fndx = 1:length(f)-1
        lst = f(fndx)+1:f(fndx+1)-1;
        if length(lst)>1
            segments = [segments; [(1:length(lst)-1)' (2:length(lst))']+size(points,1)]; %#ok<AGROW>
        elseif length(lst)==1
            segments = [segments; [size(points,1)+1 size(points,1)+1]]; %#ok<AGROW>
        end
        points = [points; x(lst) y(lst)]; %#ok<AGROW>
    end
end

% Find matching points
tolerance = sqrt(sum(range(points,1).^2))/100;
for ndx = 1:size(points,1)
    dist = sqrt(sum(bsxfun(@minus,points,points(ndx,:)).^2,2));
    f = setdiff(find(dist < tolerance),ndx);
    for fndx = 1:length(f)
        points(f(fndx),:) = NaN(1,2);
        segments(segments==f(fndx)) = ndx;
    end
end

% Remove duplicate points
[dummy dummy segments] = unique(segments(:)); %#ok<ASGLU>
segments = reshape(segments,numel(segments)/2,2);
points(isnan(points(:,1)),:) = [];

% Remove single-point segments
segments(segments(:,1)==segments(:,2),:) = [];

% Find point-segment intersections
hasChanged = 1;
while hasChanged == 1
    hasChanged = 0;

%     figure
%     for ndx = 1:size(segments,1)
%         plot(points(segments(ndx,:),1),points(segments(ndx,:),2),'color',rand(1,3))
%         hold on
%     end
    
    for ndx = 1:size(points,1)
        x2 = points(segments(:,1),:);
        x3 = points(segments(:,2),:);
        x1 = ones(size(x2,1),2);
        x1 = bsxfun(@times,x1,points(ndx,:));
        
        edgeLength_sq = sum((x3-x2).^2,2);
        alpha = sum((x1-x2).*(x3-x2),2) ./ edgeLength_sq;
        x_alpha = x2 + [alpha alpha].*(x3-x2);
        n = x_alpha - x1;
        norm_n = sqrt(sum(n.^2,2));
        f = find(norm_n<tolerance & alpha>0 & alpha<1);
        if ~isempty(f)
            f = f(1);
            if norm(x1(f,:)-x2(f,:)) > tolerance && norm(x1(f,:)-x3(f,:)) > tolerance
                segments(end+1,:) = segments(f,:); %#ok<AGROW>
                segments(f,2) = ndx;
                segments(end,1) = ndx;
                points(ndx,:) = x_alpha(f,:);
                points(ndx,:);
                hasChanged = 1;
                break
            end
        end
    end
end

%
segments = sort(segments')'; %#ok<TRSRT>
segments = unique(segments,'rows');