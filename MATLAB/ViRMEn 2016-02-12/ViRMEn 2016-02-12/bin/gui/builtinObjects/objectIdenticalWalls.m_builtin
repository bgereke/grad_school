classdef objectIdenticalWalls < virmenObject
    properties (SetObservable)
        bottom = 0;
        top = 10;
        rotation = 90;
        width = 30;
    end
    methods
        function obj = objectIdenticalWalls
            obj.iconLocations = [-10 0; 10 0];
            obj.helpString = 'Click wall centers, then press Enter';
        end
        function obj = getPoints(obj)
            [x y] = getpts(gcf);
            obj.x = x;
            obj.y = y;
        end
        function [x y z] = coords2D(obj)
            sz = obj.width*[cosd(obj.rotation) sind(obj.rotation)]/2;
            x = zeros(0,1);
            y = zeros(0,1);
            z = zeros(0,1);
            loc = obj.locations;
            for ndx = 1:size(loc,1)
                x = [x; loc(ndx,1)+sz(1)*[-1 1 1 -1 -1 1 NaN 1 -1 NaN]']; %#ok<AGROW>
                y = [y; loc(ndx,2)+sz(2)*[-1 1 1 -1 -1 1 NaN 1 -1 NaN]']; %#ok<AGROW>
                z = [z; obj.bottom; obj.bottom; obj.top; obj.top; obj.bottom; obj.top; NaN; obj.bottom; obj.top; NaN]; %#ok<AGROW>
            end
            x(end) = [];
            y(end) = [];
            z(end) = [];
        end
        function objSurface = coords3D(obj)
            texture = tile(obj.texture,obj.tiling);
            
            x = texture.triangles.vertices(:,1);
            y = texture.triangles.vertices(:,2);
            x = ((x-min(x(:)))/range(x(:))-0.5)*obj.width;
            z = (y-min(y(:)))/range(y(:))*(obj.top-obj.bottom)+obj.bottom;
            y = zeros(size(x));
            rot = [cosd(obj.rotation) sind(obj.rotation); -sind(obj.rotation) cosd(obj.rotation)];
            loc = [x y]*rot;
            x = loc(:,1);
            y = loc(:,2);
            objSurface.vertices = zeros(0,3);
            objSurface.triangulation = zeros(0,3);
            loc = obj.locations;
            for ndx = 1:size(loc,1)
                objSurface.triangulation = [objSurface.triangulation; texture.triangles.triangulation+size(objSurface.vertices,1)];
                objSurface.vertices = [objSurface.vertices; bsxfun(@plus,[x y z],[loc(ndx,:) 0])];
            end
            objSurface.cdata = repmat(texture.triangles.cdata,size(loc,1),1);
        end
        function edges = edges(obj)
            edges = zeros(0,4);
            sz = obj.width*[cosd(obj.rotation) sind(obj.rotation)]/2;
            loc = obj.locations;
            for ndx = 1:size(loc,1)
                edges = [edges; loc(ndx,:)-sz loc(ndx,:)+sz]; %#ok<AGROW>
            end
        end
    end
end