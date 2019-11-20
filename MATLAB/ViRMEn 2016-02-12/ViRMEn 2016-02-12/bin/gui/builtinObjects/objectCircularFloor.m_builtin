classdef objectCircularFloor < virmenObject
    properties (SetObservable)
        radius = 100;
        elevation = 0;
        polarTiling = 0;
    end
    methods
        function obj = objectCircularFloor
            obj.iconLocations = [0 0];
            obj.helpString = 'Click floor centers, then press Enter';
        end
        function obj = getPoints(obj)
            [x y] = getpts(gcf);
            obj.x = x;
            obj.y = y;
        end
        function [x y z] = coords2D(obj)
            theta = linspace(0,2*pi,1000);
            eachx = obj.radius*[cos(theta) NaN cos(pi/4) cos(5*pi/4) NaN cos(3*pi/4) cos(7*pi/4)]';
            eachy = obj.radius*[sin(theta) NaN sin(pi/4) sin(5*pi/4) NaN sin(3*pi/4) sin(7*pi/4)]';
            
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
            
            if ~obj.polarTiling
                isbad = sqrt(sum(bsxfun(@minus,bsxfun(@rdivide,texture.triangles.vertices,max(texture.triangles.vertices,[],1)),[.5 .5]).^2,2))>.5;
                badtri = all(isbad(texture.triangles.triangulation),2);
                texture.triangles.triangulation(badtri,:) = [];
                
                notused = setdiff((1:size(texture.triangles.vertices,1))',texture.triangles.triangulation(:));
                numbef = zeros(size(texture.triangles.vertices,1),1);
                numbef(notused) = 1;
                numbef = cumsum(numbef);
                texture.triangles.vertices(notused,:) = [];
                texture.triangles.cdata(notused,:) = [];
                texture.triangles.triangulation = texture.triangles.triangulation-numbef(texture.triangles.triangulation);
                x = texture.triangles.vertices(:,1);
                y = texture.triangles.vertices(:,2);
                x = ((x-min(x(:)))/range(x(:))-0.5)*2*obj.radius;
                y = ((y-min(y(:)))/range(y(:))-0.5)*2*obj.radius;
            else
                theta = 2*pi*(texture.triangles.vertices(:,1)-min(texture.triangles.vertices(:,1)))/range(texture.triangles.vertices(:,1));
                r = obj.radius*(texture.triangles.vertices(:,2)-min(texture.triangles.vertices(:,2)))/range(texture.triangles.vertices(:,2));
                x = r.*cos(theta);
                y = r.*sin(theta);
            end
            
            
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
            loc = obj.locations;
            edges = [loc loc];
        end
    end
end