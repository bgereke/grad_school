function obj = computeTexture(obj)

% Refining parameters
rv = obj.height*obj.refining(1);
rh = obj.width*obj.refining(2);
rvf = obj.height*obj.grid(1);
rhf = obj.width*obj.grid(2);

% Obtain sketch points and segments
[objPoints objSegments] = obj.getSegments;

% Eliminate lonely points
objSegments(objSegments(:,1)==objSegments(:,2),:) = [];
f = 1:size(objPoints,1);
f(objSegments(:)) = [];
for ndx = length(f):-1:1
    objPoints(f(ndx),:) = [];
    g = objSegments(:)>f(ndx);
    objSegments(g) = objSegments(g)-1;
end
objSegments = unique(objSegments,'rows');

allPoints = objPoints;
sg = objSegments;
if obj.tilable(2)==1
    tolerance = 0.01*range(objPoints(:,1));
    y1 = allPoints(allPoints(:,1)==0,2);
    y2 = allPoints(allPoints(:,1)==obj.width,2);
    allY = unique([y1; y2]);
    newy = unique(y2);
    for ndx = length(newy):-1:1
        if ~isempty(find(abs(newy(ndx)-y1)<tolerance, 1))
            newy(ndx) = [];
        end
    end
    for ndx = 1:length(newy)
        f = find(allPoints(sg(:,1),1)==0 & allPoints(sg(:,2),1)==0 & ...
            ((allPoints(sg(:,1),2) < newy(ndx) & allPoints(sg(:,2),2) > newy(ndx)) | ...
            (allPoints(sg(:,1),2) > newy(ndx) & allPoints(sg(:,2),2) < newy(ndx))));
        allPoints(end+1,:) = [0 newy(ndx)]; %#ok<AGROW>
        sg(end+1,:) = [size(allPoints,1) sg(f,2)]; %#ok<AGROW>
        sg(f,2) = size(allPoints,1);
    end
    
    
    newy = setdiff(allY,y2);
    for ndx = length(newy):-1:1
        if ~isempty(find(abs(newy(ndx)-y2)<tolerance, 1))
            newy(ndx) = [];
        end
    end
    for ndx = 1:length(newy)
        f = find(allPoints(sg(:,1),1)==obj.width & allPoints(sg(:,2),1)==obj.width & ...
            ((allPoints(sg(:,1),2) < newy(ndx) & allPoints(sg(:,2),2) > newy(ndx)) | ...
            (allPoints(sg(:,1),2) > newy(ndx) & allPoints(sg(:,2),2) < newy(ndx))));
        allPoints(end+1,:) = [obj.width newy(ndx)]; %#ok<AGROW>
        sg(end+1,:) = [size(allPoints,1) sg(f,2)]; %#ok<AGROW>
        sg(f,2) = size(allPoints,1);
    end
end
if obj.tilable(1)==1
    tolerance = 0.01*range(objPoints(:,2));
    x1 = allPoints(allPoints(:,2)==0,1);
    x2 = allPoints(allPoints(:,2)==obj.height,1);
    allX = unique([x1; x2]);
    newx = setdiff(allX,x1);
    for ndx = length(newx):-1:1
        if ~isempty(find(abs(newx(ndx)-x1)<tolerance, 1))
            newx(ndx) = [];
        end
    end
    for ndx = 1:length(newx)
        f = find(allPoints(sg(:,1),2)==0 & allPoints(sg(:,2),2)==0 & ...
            ((allPoints(sg(:,1),1) < newx(ndx) & allPoints(sg(:,2),1) > newx(ndx)) | ...
            (allPoints(sg(:,1),1) > newx(ndx) & allPoints(sg(:,2),1) < newx(ndx))));
        allPoints(end+1,:) = [newx(ndx) 0]; %#ok<AGROW>
        sg(end+1,:) = [size(allPoints,1) sg(f,2)]; %#ok<AGROW>
        sg(f,2) = size(allPoints,1);
    end
    
    newx = setdiff(allX,x2);
    for ndx = length(newx):-1:1
        if ~isempty(find(abs(newx(ndx)-x2)<tolerance, 1))
            newx(ndx) = [];
        end
    end
    for ndx = 1:length(newx)
        f = find(allPoints(sg(:,1),2)==obj.height & allPoints(sg(:,2),2)==obj.height & ...
            ((allPoints(sg(:,1),1) < newx(ndx) & allPoints(sg(:,2),1) > newx(ndx)) | ...
            (allPoints(sg(:,1),1) > newx(ndx) & allPoints(sg(:,2),1) < newx(ndx))));
        allPoints(end+1,:) = [newx(ndx) obj.height]; %#ok<AGROW>
        sg(end+1,:) = [size(allPoints,1) unique(sg(f,2))]; %#ok<AGROW>
        sg(f,2) = size(allPoints,1);
    end
end

% Grid points
xgrid = (0:obj.width/rhf:obj.width)';
ygrid = (0:obj.height/rvf:obj.height)';

% Refine edges
w = abs(allPoints(sg(:,1),1)-allPoints(sg(:,2),1));
h = abs(allPoints(sg(:,1),2)-allPoints(sg(:,2),2));
sg2 = zeros(0,2);
for ndx = 1:size(sg,1)
    numseg = max(ceil([w(ndx)/rh h(ndx)/rv]));
    if numseg==0
        numseg = 1;
    end
    alpha = linspace(0,1,numseg+1)';
    alpha = alpha(2:end-1);
    
    % Place grid points on segments
    f = (xgrid > allPoints(sg(ndx,1),1) & xgrid < allPoints(sg(ndx,2),1)) | (xgrid < allPoints(sg(ndx,1),1) & xgrid > allPoints(sg(ndx,2),1));
    alpha = [alpha; (xgrid(f)-allPoints(sg(ndx,1),1))./(allPoints(sg(ndx,2),1)-allPoints(sg(ndx,1),1))]; %#ok<AGROW>
    f = (ygrid > allPoints(sg(ndx,1),2) & ygrid < allPoints(sg(ndx,2),2)) | (ygrid < allPoints(sg(ndx,1),2) & ygrid > allPoints(sg(ndx,2),2));
    alpha = [alpha; (ygrid(f)-allPoints(sg(ndx,1),2))./(allPoints(sg(ndx,2),2)-allPoints(sg(ndx,1),2))]; %#ok<AGROW>

    alpha = unique(alpha);
    
    x = allPoints(sg(ndx,1),1) + alpha*(allPoints(sg(ndx,2),1)-allPoints(sg(ndx,1),1));
    y = allPoints(sg(ndx,1),2) + alpha*(allPoints(sg(ndx,2),2)-allPoints(sg(ndx,1),2));
    
    indx = [sg(ndx,1) size(allPoints,1)+(1:length(x)) sg(ndx,2)]';
    allPoints = [allPoints; x y]; %#ok<AGROW>
    sg2 = [sg2; indx(1:end-1) indx(2:end)]; %#ok<AGROW>
end

% Add grid points
[x y] = meshgrid(xgrid,ygrid);
x = x(2:end-1,2:end-1);
y = y(2:end-1,2:end-1);
x = x(:);
y = y(:);
allPoints = [allPoints; [x y]];

% Add segments connecting grid points
sgnew = zeros(0,2);
for ndx = 1:length(xgrid)
    f = find(allPoints(:,1)==xgrid(ndx));
    [~,ord] = sort(allPoints(f,2));
    f = f(ord);
    sgnew = [sgnew; [f(1:end-1) f(2:end)]]; %#ok<AGROW>
end
for ndx = 1:length(ygrid)
    f = find(allPoints(:,2)==ygrid(ndx));
    [~,ord] = sort(allPoints(f,1));
    f = f(ord);
    sgnew = [sgnew; [f(1:end-1) f(2:end)]]; %#ok<AGROW>
end
sgnew = sort(sgnew')'; %#ok<TRSRT>

% Triangulate texture
allPoints = round(allPoints*100000)/100000;
[allPoints j indx] = unique(allPoints,'rows'); %#ok<ASGLU>
sg2 = indx(sg2);
sgnew = indx(sgnew);
sg2 = sort(sg2')'; %#ok<TRSRT>
sg2 = unique(sg2,'rows');
tri = DelaunayTri(allPoints(:,1),allPoints(:,2));
warning off %#ok<WNOFF>
tri.Constraints = unique(sort([sg2; sgnew]')','rows'); %#ok<TRSRT>
warning on %#ok<WNON>


% Obtain list of colors
colorCoords = [];
colorRGBA = [];
for ndx = 1:length(obj.shapes)
    if strcmp(class(obj.shapes{ndx}),'shapeColor')
        loc = obj.shapes{ndx}.locations;
        colorCoords = [colorCoords; loc]; %#ok<AGROW>
        colorRGBA = [colorRGBA; repmat(obj.shapes{ndx}.RGBA,size(loc,1),1)]; %#ok<AGROW>
    end
end

% Figure out the triangles that contain color locations
cdata = NaN(size(tri.Triangulation,1),4);
if ~isempty(colorCoords)
    ptLoc = pointLocation(tri,colorCoords);
    cdata(ptLoc,:) = colorRGBA;
end

% Color the rest of the triangles
trias = zeros(0,2);
edg = edges(tri);
pairs = edgeAttachments(tri,edg);
for ndx = 1:length(pairs)
    if length(pairs{ndx}) == 1
        continue
    end
        
    % Test if the edge is real
    f = find((sg2(:,1)==edg(ndx,1) & sg2(:,2)==edg(ndx,2)) | ...
        (sg2(:,2)==edg(ndx,1) & sg2(:,1)==edg(ndx,2)),1);
    
    if isempty(f)
        trias(end+1,:) = pairs{ndx}; %#ok<AGROW>
    end
end

hasChanged = 1;
while hasChanged == 1
    hasChanged = 0; 
    for ndx = 1:size(trias,1)
        if isnan(cdata(trias(ndx,1),1)) && ~isnan(cdata(trias(ndx,2),1))
            cdata(trias(ndx,1),:) = cdata(trias(ndx,2),:);
            hasChanged = 1;
        end
        if isnan(cdata(trias(ndx,2),1)) && ~isnan(cdata(trias(ndx,1),1))
            cdata(trias(ndx,2),:) = cdata(trias(ndx,1),:);
            hasChanged = 1;
        end
    end
end

vertices = tri.X;
triangulation = tri.Triangulation;

% Orient triangles away from edges
cdata(isnan(cdata)) = -1;
vatt = vertexAttachments(tri,(1:size(tri.X,1))');
isall = cellfun(@(x)all(cdata(x(1),1)==cdata(x,1)) & all(cdata(x(1),2)==cdata(x,2)) & all(cdata(x(1),3)==cdata(x,3)),vatt);
sg = find(~isall);
oth = {[2 3],[1 3],[1 2]};
edgeTriangles = [];
vertexColor = NaN(size(vertices,1),4);
for ndx = 1:size(triangulation)
    allEdge = true;
    for j = 1:3
        if isempty(find(sg==triangulation(ndx,j),1))
            triangulation(ndx,:) = triangulation(ndx,[j oth{j}]);
            vertexColor(triangulation(ndx,1),:) = cdata(ndx,:);
            allEdge = false;
            break
        end
    end
    if allEdge
        edgeTriangles(end+1) = ndx; %#ok<AGROW>
    end
end

for ndx = 1:length(edgeTriangles)
    tria = triangulation(edgeTriangles(ndx),:);
    vertices(end+1,:) = vertices(tria(1),:); %#ok<AGROW>
    tria(1) = size(vertices,1);
    triangulation(edgeTriangles(ndx),:) = tria;
    vertexColor(end+1,:) = cdata(edgeTriangles(ndx),:); %#ok<AGROW>
end

for ndx = size(vertices,1):-1:1
    if all(isnan(vertexColor(ndx,:)))
        f = find(~all(isnan(vertexColor),2) & vertices(:,1)==vertices(ndx,1) & vertices(:,2)==vertices(ndx,2),1);
        if ~isempty(f)
            vertices(ndx,:) = [];
            vertexColor(ndx,:) = [];
            triangulation(triangulation==ndx) = f;
            triangulation(triangulation>ndx) = triangulation(triangulation>ndx)-1;
        end
    end
end

for ndx = size(vertices,1):-1:1
    if isempty(find(triangulation==ndx,1))
        vertices(ndx,:) = [];
        vertexColor(ndx,:) = [];
        triangulation(triangulation>ndx) = triangulation(triangulation>ndx)-1;
    end
end

obj.triangles.vertices = vertices;
obj.triangles.triangulation = triangulation;
vertexColor(vertexColor < 0) = NaN;
% vertexColor(isnan(vertexColor(:,1:3))) = 0;
obj.triangles.cdata = vertexColor*(1-eps);