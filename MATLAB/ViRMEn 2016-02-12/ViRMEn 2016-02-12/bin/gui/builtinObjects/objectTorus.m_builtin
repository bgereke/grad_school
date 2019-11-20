classdef objectTorus < virmenObject
    properties (SetObservable)
        minorRadius = 2;
        majorRadius = 10;
        elevation = 0;
        tilt = 0;
        zRotation = 0;
    end
    methods
        function obj = objectTorus
            obj.iconLocations = [0 0];
            obj.helpString = 'Click torus centers, then press Enter';
        end
        function obj = getPoints(obj)
            [obj.x obj.y] = getpts(gcf);
        end
        function [x y z] = coords2D(obj)
            loc = obj.locations;
            
            [theta phi] = meshgrid(linspace(0,2*pi,20),linspace(0,2*pi,10));
            
            xgrid = (obj.majorRadius+obj.minorRadius*cos(phi)).*cos(theta);
            ygrid = (obj.majorRadius+obj.minorRadius*cos(phi)).*sin(theta);
            zgrid = obj.minorRadius*sin(phi);
            
            pt = zeros(0,3);
            for ndx = 1:size(xgrid,1)
                pt = [pt; [xgrid(ndx,:)' ygrid(ndx,:)' zgrid(ndx,:)']; [NaN NaN NaN]]; %#ok<AGROW>
            end
            for ndx = 1:size(xgrid,2)
                pt = [pt; [xgrid(:,ndx) ygrid(:,ndx) zgrid(:,ndx)]; [NaN NaN NaN]]; %#ok<AGROW>
            end
            
            rotTilt = [cosd(obj.tilt) -sind(obj.tilt); sind(obj.tilt) cosd(obj.tilt)];
            pt(:,[1 3]) = (rotTilt*pt(:,[1 3])')';
            
            rotZ = [cosd(obj.zRotation) -sind(obj.zRotation); sind(obj.zRotation) cosd(obj.zRotation)];
            pt(:,[1 2]) = (rotZ*pt(:,[1 2])')';
            
            elev = obj.elevation;
            if length(elev) == 1
                elev = repmat(elev,size(loc,1),1);
            end
            elev = [elev; zeros(size(loc,1)-length(elev),1)];
            
            x = zeros(0,1);
            y = zeros(0,1);
            z = zeros(0,1);
            for ndx = 1:size(loc,1)
                x = [x; pt(:,1)+loc(ndx,1)]; %#ok<AGROW>
                y = [y; pt(:,2)+loc(ndx,2)]; %#ok<AGROW>
                z = [z; pt(:,3)+elev(ndx)]; %#ok<AGROW>
            end
            
            x(end) = [];
            y(end) = [];
            z(end) = [];
        end
        function objSurface = coords3D(obj)
            texture = tile(obj.texture,obj.tiling);
            loc = obj.locations;
            
            theta = 2*pi*(texture.triangles.vertices(:,1)-min(texture.triangles.vertices(:,1)))/range(texture.triangles.vertices(:,1));
            phi = 2*pi*(texture.triangles.vertices(:,2)-min(texture.triangles.vertices(:,2)))/range(texture.triangles.vertices(:,2));
            pt(:,1) = (obj.majorRadius+obj.minorRadius*cos(phi)).*cos(theta);
            pt(:,2) = (obj.majorRadius+obj.minorRadius*cos(phi)).*sin(theta);
            pt(:,3) = obj.minorRadius*sin(phi);
            
            rotTilt = [cosd(obj.tilt) -sind(obj.tilt); sind(obj.tilt) cosd(obj.tilt)];
            pt(:,[1 3]) = (rotTilt*pt(:,[1 3])')';
            
            rotZ = [cosd(obj.zRotation) -sind(obj.zRotation); sind(obj.zRotation) cosd(obj.zRotation)];
            pt(:,[1 2]) = (rotZ*pt(:,[1 2])')';
            
            elev = obj.elevation;
            if length(elev) == 1
                elev = repmat(elev,size(loc,1),1);
            end
            elev = [elev; zeros(size(loc,1)-length(elev),1)];
            
            objSurface.vertices = zeros(0,3);
            objSurface.triangulation = zeros(0,3);
            for ndx = 1:size(loc,1)
                objSurface.triangulation = [objSurface.triangulation; texture.triangles.triangulation+size(objSurface.vertices,1)];
                objSurface.vertices = [objSurface.vertices; bsxfun(@plus,pt,[loc(ndx,:) elev(ndx)])];
            end
            objSurface.cdata = repmat(texture.triangles.cdata,size(loc,1),1);
        end
        function edges = edges(obj)
            loc = obj.locations;
            edges = [loc loc];
        end
    end
end