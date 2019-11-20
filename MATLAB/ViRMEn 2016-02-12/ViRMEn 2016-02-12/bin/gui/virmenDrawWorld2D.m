function [surfaceHandles, edgeHandles, playerHandle] = virmenDrawWorld2D(worldData)

% Plot surfaces
surfaceHandles = zeros(1,size(worldData.objectVertices,1));
for ndx = 1:size(worldData.objectVertices,1)
    lst1 = worldData.objectTriangles(ndx,1);
    lst2 = worldData.objectTriangles(ndx,2);
    cons = worldData.triangulation(lst1:lst2,:) - worldData.objectVertices(ndx,1) + 1;
    lst1 = worldData.objectVertices(ndx,1);
    lst2 = worldData.objectVertices(ndx,2);
    pts = worldData.points(lst1:lst2,:);
    hold on
    x = reshape(pts(cons,1),size(cons))';
    y = reshape(pts(cons,2),size(cons))';
    x(end+1,:) = NaN; %#ok<AGROW>
    y(end+1,:) = NaN; %#ok<AGROW>

    surfaceHandles(ndx) = plot(x(:),y(:),'k');
end

% Plot edges
edgeHandles = zeros(1,size(worldData.edgeGroups,1));
for eg = 1:size(worldData.edgeGroups,1)
    for ndx = worldData.edgeGroups(eg,1):worldData.edgeGroups(eg,2)
        r = worldData.radius(ndx);
        x = worldData.endpoints(ndx,1)+r*cos(linspace(0,2*pi,100));
        y = worldData.endpoints(ndx,2)+r*sin(linspace(0,2*pi,100));
        x = [x NaN worldData.endpoints(ndx,3)+r*cos(linspace(0,2*pi,100))]; %#ok<AGROW>
        y = [y NaN worldData.endpoints(ndx,4)+r*sin(linspace(0,2*pi,100))]; %#ok<AGROW>
        ang = atan2(worldData.endpoints(ndx,4)-worldData.endpoints(ndx,2), ...
            worldData.endpoints(ndx,3)-worldData.endpoints(ndx,1));
        ang = ang+pi/2;
        x = [x NaN worldData.endpoints(ndx,[1 3])+r*cos(ang)]; %#ok<AGROW>
        y = [y NaN worldData.endpoints(ndx,[2 4])+r*sin(ang)]; %#ok<AGROW>
        x = [x NaN worldData.endpoints(ndx,[1 3])-r*cos(ang)]; %#ok<AGROW>
        y = [y NaN worldData.endpoints(ndx,[2 4])-r*sin(ang)]; %#ok<AGROW>
        plot(x,y,':k');
    end
end

% Plot start location
loc = [0 1; 1 1; 1 -1; -1 -1; -1 1; 0 1; 0 4; 1 2; NaN NaN; -1 2; 0 4];
rot = [cos(worldData.defaultPosition(4)) -sin(worldData.defaultPosition(4)); ...
    sin(worldData.defaultPosition(4)) cos(worldData.defaultPosition(4))];
loc = rot*loc';
playerHandle = plot(loc(1,:)+worldData.defaultPosition(1),loc(2,:)+worldData.defaultPosition(2),'k','linewidth',2);