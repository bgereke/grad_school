function tiled = tileTexture(textureToTile,tiling)

tiling = round(tiling);

texture = struct;
props = properties(textureToTile);
for ndx = 1:length(props)
    texture.(props{ndx}) = textureToTile.(props{ndx});
end

% Determine tiling
if texture.tilable(1) == 1
    topNdx = find(texture.triangles.vertices(:,2)==0);
    [~,unq] = unique(texture.triangles.vertices(topNdx,1));
    topNdx = topNdx(unq);
    bottomNdx = find(texture.triangles.vertices(:,2)==texture.height);
    [~,unq] = unique(texture.triangles.vertices(bottomNdx,1));
    bottomNdx = bottomNdx(unq);
    [y ord] = sort(texture.triangles.vertices(topNdx,1)); %#ok<*ASGLU>
    topNdx = topNdx(ord);
    [y ord] = sort(texture.triangles.vertices(bottomNdx,1));
    bottomNdx = bottomNdx(ord);
else
    topNdx = [];
    bottomNdx = [];
end
if texture.tilable(2) == 1
    leftNdx = find(texture.triangles.vertices(:,1)==0);
    [~,unq] = unique(texture.triangles.vertices(leftNdx,2));
    leftNdx = leftNdx(unq);
    rightNdx = find(texture.triangles.vertices(:,1)==texture.width);
    [~,unq] = unique(texture.triangles.vertices(rightNdx,2));
    rightNdx = rightNdx(unq);
    [x ord] = sort(texture.triangles.vertices(leftNdx,2));
    leftNdx = leftNdx(ord);
    [x ord] = sort(texture.triangles.vertices(rightNdx,2));
    rightNdx = rightNdx(ord);
else
    leftNdx = [];
    rightNdx = [];
end
for ndx = length(topNdx):-1:1
    if sum(texture.triangles.cdata(topNdx(ndx),:) == texture.triangles.cdata(bottomNdx(ndx),:))<3  && ~isnan(texture.triangles.cdata(topNdx(ndx),1))
        topNdx(ndx) = [];
        bottomNdx(ndx) = [];
    end
end

% Tile
tiled.triangles.triangulation = zeros(0,3);
tiled.triangles.vertices = zeros(0,2);
tiled.triangles.cdata = zeros(0,4);
allTop = zeros(0,1);
allBottom = zeros(0,1);
allLeft = zeros(0,1);
allRight = zeros(0,1);
for y = 0:tiling(1)-1
    if y < tiling(1)-1
        allBottom = [allBottom; size(tiled.triangles.vertices,1)+bottomNdx]; %#ok<*AGROW>
    end
    if y > 0
        allTop = [allTop; size(tiled.triangles.vertices,1)+topNdx];
    end
    allLeft = [allLeft; size(tiled.triangles.vertices,1)+leftNdx];
    allRight = [allRight; size(tiled.triangles.vertices,1)+rightNdx];
    
    tiled.triangles.triangulation = [tiled.triangles.triangulation; texture.triangles.triangulation+size(tiled.triangles.vertices,1)];
    tiled.triangles.vertices = [tiled.triangles.vertices; ...
        [texture.triangles.vertices(:,1) texture.triangles.vertices(:,2)+texture.height*y]];
    tiled.triangles.cdata = [tiled.triangles.cdata; texture.triangles.cdata];
end

for ndx = 1:length(allTop)
    tiled.triangles.triangulation(tiled.triangles.triangulation==allTop(ndx)) = -allBottom(ndx);
    allLeft(allLeft==allTop(ndx)) = -allBottom(ndx);
    allRight(allRight==allTop(ndx)) = -allBottom(ndx);
end

tiled.triangles.triangulation = abs(tiled.triangles.triangulation);
allLeft = abs(allLeft);
allRight = abs(allRight);
f = find(allLeft(2:end)==allLeft(1:end-1));
allLeft(f) = [];
allRight(f) = [];

toDelete = zeros(max(tiled.triangles.triangulation(:)),1);
toDelete(allTop) = 1;
toDelete = cumsum(toDelete);
tiled.triangles.triangulation = tiled.triangles.triangulation-toDelete(tiled.triangles.triangulation);
allLeft = allLeft-toDelete(allLeft);
allRight = allRight-toDelete(allRight);
tiled.triangles.vertices(allTop,:) = [];
tiled.triangles.cdata(allTop,:) = [];

% Horizontal tiling
leftNdx = allLeft;
allLeft = zeros(0,1);
rightNdx = allRight;
allRight = zeros(0,1);
texture.triangles.triangulation = tiled.triangles.triangulation;
tiled.triangles.triangulation = zeros(0,3);
texture.triangles.vertices = tiled.triangles.vertices;
tiled.triangles.vertices = zeros(0,2);
texture.triangles.cdata = tiled.triangles.cdata;
tiled.triangles.cdata = zeros(0,4);

for ndx = length(leftNdx):-1:1
    if sum(texture.triangles.cdata(leftNdx(ndx),:) == texture.triangles.cdata(rightNdx(ndx),:))<3 && ~isnan(texture.triangles.cdata(leftNdx(ndx),1))
        leftNdx(ndx) = [];
        rightNdx(ndx) = [];
    end
end

for x = 0:tiling(2)-1
    if x < tiling(2)-1
        allRight = [allRight; size(tiled.triangles.vertices,1)+rightNdx]; %#ok<*AGROW>
    end
    if x > 0
        allLeft = [allLeft; size(tiled.triangles.vertices,1)+leftNdx];
    end
    
    tiled.triangles.triangulation = [tiled.triangles.triangulation; texture.triangles.triangulation+size(tiled.triangles.vertices,1)];
    tiled.triangles.vertices = [tiled.triangles.vertices; ...
        [texture.triangles.vertices(:,1)+texture.width*x texture.triangles.vertices(:,2)]];
    tiled.triangles.cdata = [tiled.triangles.cdata; texture.triangles.cdata];
end

if ~isempty(allLeft)
    tiled.triangles.triangulation = arrayReplace(tiled.triangles.triangulation,allLeft,allRight);
end
toDelete = zeros(max(tiled.triangles.triangulation(:)),1);
toDelete(allLeft) = 1;
toDelete = cumsum(toDelete);


tiled.triangles.triangulation = tiled.triangles.triangulation-toDelete(tiled.triangles.triangulation);
tiled.triangles.vertices(allLeft,:) = [];
tiled.triangles.cdata(allLeft,:) = [];