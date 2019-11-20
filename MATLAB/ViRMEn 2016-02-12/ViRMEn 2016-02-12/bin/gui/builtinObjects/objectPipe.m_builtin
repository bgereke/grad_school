classdef objectPipe < virmenObject
    properties (SetObservable)
        radius = 2;
        elevation1 = 0;
        elevation2 = 0;
    end
    methods
        function obj = objectPipe
            obj.iconLocations = [0 -5; 10 -10; 15 -20];
            obj.helpString = 'Click pipe endpoints, then press Enter';
        end
        function obj = getPoints(obj)
            [x y] = getline(gcf);
            if length(x)==1
                obj.x = [];
                obj.y = [];
            else
                obj.x = x;
                obj.y = y;
            end
        end
        function [x y z] = coords2D(obj)
            loc = obj.locations;
            
            [theta zgrid] = meshgrid(linspace(0,2*pi,10),linspace(0,1,10));
            xgrid = obj.radius*cos(theta);
            ygrid = obj.radius*sin(theta);
            
            ptb = zeros(0,3);
            for ndx = 1:size(xgrid,1)
                ptb = [ptb; [xgrid(ndx,:)' ygrid(ndx,:)' zgrid(ndx,:)']; [NaN NaN NaN]]; %#ok<AGROW>
            end
            for ndx = 1:size(xgrid,2)
                ptb = [ptb; [xgrid(:,ndx) ygrid(:,ndx) zgrid(:,ndx)]; [NaN NaN NaN]]; %#ok<AGROW>
            end
            
            elev = linspace(obj.elevation1,obj.elevation2,size(loc,1));
            
            x = zeros(0,1);
            y = zeros(0,1);
            z = zeros(0,1);
            for ndx = 1:size(loc,1)-1
                pt = ptb;
                pt(:,3) = pt(:,3)*sqrt(sum((loc(ndx+1,:)-loc(ndx,:)).^2)+(elev(ndx+1)-elev(ndx))^2);
                
                ang = atan2(elev(ndx+1)-elev(ndx),norm(loc(ndx+1,:)-loc(ndx,:)))-pi/2;
                rotTilt = [cos(ang) -sin(ang); sin(ang) cos(ang)];
                pt(:,[1 3]) = (rotTilt*pt(:,[1 3])')';
                
                ang = atan2(loc(ndx+1,2)-loc(ndx,2),loc(ndx+1,1)-loc(ndx,1));
                rotZ = [cos(ang) -sin(ang); sin(ang) cos(ang)];
                pt(:,[1 2]) = (rotZ*pt(:,[1 2])')';
                
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
            
            xs = obj.radius*cos(theta);
            ys = obj.radius*sin(theta);
            zs = (texture.triangles.vertices(:,2)-min(texture.triangles.vertices(:,2)))/range(texture.triangles.vertices(:,2));
            
            elev = linspace(obj.elevation1,obj.elevation2,size(loc,1));
            
            objSurface.vertices = zeros(0,3);
            objSurface.triangulation = zeros(0,3);
            for ndx = 1:size(loc,1)-1
                pt(:,1) = xs;
                pt(:,2) = ys;
                pt(:,3) = zs*sqrt(sum((loc(ndx+1,:)-loc(ndx,:)).^2)+(elev(ndx+1)-elev(ndx))^2);
                
                ang = atan2(elev(ndx+1)-elev(ndx),norm(loc(ndx+1,:)-loc(ndx,:)))-pi/2;
                rotTilt = [cos(ang) -sin(ang); sin(ang) cos(ang)];
                pt(:,[1 3]) = (rotTilt*pt(:,[1 3])')';
                
                ang = atan2(loc(ndx+1,2)-loc(ndx,2),loc(ndx+1,1)-loc(ndx,1));
                rotZ = [cos(ang) -sin(ang); sin(ang) cos(ang)];
                pt(:,[1 2]) = (rotZ*pt(:,[1 2])')';
                
                objSurface.triangulation = [objSurface.triangulation; texture.triangles.triangulation+size(objSurface.vertices,1)];
                objSurface.vertices = [objSurface.vertices; bsxfun(@plus,pt,[loc(ndx,:) elev(ndx)])];
            end
            objSurface.cdata = repmat(texture.triangles.cdata,size(loc,1)-1,1);
        end
        function edges = edges(obj)
            loc = obj.locations;
            edges = [loc(1:end-1,:) loc(2:end,:)];
        end
    end
end