classdef objectFloor < virmenObject
    properties (SetObservable)
        width = 100;
        height = 100;
        elevation = 0;
        rotation = 0;
    end
    methods
        function obj = objectFloor
            obj.iconLocations = [0 0];
            obj.helpString = 'Click and drag a floor rectangle';
        end
        function obj = getPoints(obj)
            rect = getrect(gcf);
            obj.x = rect(1)+rect(3)/2;
            obj.y = rect(2)+rect(4)/2;
            obj.width = rect(3);
            obj.height = rect(4);
        end
        function [x y z] = coords2D(obj)
            rot = [cosd(obj.rotation) sind(obj.rotation); -sind(obj.rotation) cosd(obj.rotation)];
            x = [-obj.width/2; -obj.width/2; obj.width/2; obj.width/2; -obj.width/2; obj.width/2; obj.width/2; -obj.width/2];
            y = [-obj.height/2; obj.height/2; obj.height/2; -obj.height/2; -obj.height/2; obj.height/2; -obj.height/2; obj.height/2];
            xy = [x y]*rot;
            eachx = xy(:,1);
            eachy = xy(:,2);
            
            loc = obj.locations;
            x = zeros(0,1);
            y = zeros(0,1);
            for ndx = 1:size(loc,1)
                x = [x; NaN; eachx+loc(ndx,1)]; %#ok<AGROW>
                y = [y; NaN; eachy+loc(ndx,2)]; %#ok<AGROW>
            end
            x(1) = [];
            y(1) = [];
            z = obj.elevation*ones(size(x));
        end
        function objSurface = coords3D(obj)
            texture = tile(obj.texture,obj.tiling);
            
            x = texture.triangles.vertices(:,1);
            y = texture.triangles.vertices(:,2);
            x = ((x-min(x(:)))/range(x(:))-0.5)*obj.width;
            y = ((y-min(y(:)))/range(y(:))-0.5)*obj.height;
            rot = [cosd(obj.rotation) sind(obj.rotation); -sind(obj.rotation) cosd(obj.rotation)];
            loc = [x y]*rot;
            x = loc(:,1);
            y = loc(:,2);
            objSurface.vertices = zeros(0,3);
            objSurface.triangulation = zeros(0,3);
            loc = obj.locations;
            for ndx = 1:size(loc,1)
                objSurface.triangulation = [objSurface.triangulation; texture.triangles.triangulation+size(objSurface.vertices,1)];
                objSurface.vertices = [objSurface.vertices; bsxfun(@plus,[x y],loc(ndx,:)) repmat(obj.elevation,size(x,1),1)];
            end
            objSurface.cdata = repmat(texture.triangles.cdata,size(loc,1),1);
        end
        function edges = edges(obj)
            edges = zeros(0,4);
            loc = obj.locations;
            rot = [cosd(obj.rotation) sind(obj.rotation); -sind(obj.rotation) cosd(obj.rotation)];
            for ndx = 1:size(loc,1)
                x = [loc(ndx,1)-obj.width/2; loc(ndx,1)-obj.width/2; loc(ndx,1)+obj.width/2; loc(ndx,1)+obj.width/2; loc(ndx,1)-obj.width/2];
                y = [loc(ndx,2)-obj.height/2; loc(ndx,2)+obj.height/2; loc(ndx,2)+obj.height/2; loc(ndx,2)-obj.height/2; loc(ndx,2)-obj.height/2];
                xy = [x y]*rot;
                x = xy(:,1);
                y = xy(:,2);
                edges = [edges; x(1:end-1) y(1:end-1) x(2:end) y(2:end)]; %#ok<AGROW>
            end
        end
    end
end