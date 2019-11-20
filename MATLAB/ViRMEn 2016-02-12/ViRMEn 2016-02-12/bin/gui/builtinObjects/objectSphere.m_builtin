classdef objectSphere < virmenObject
    properties (SetObservable)
        radius = 10;
        elevation = 0;
    end
    methods
        function obj = objectSphere
            obj.iconLocations = [0 0];
            obj.helpString = 'Click sphere centers, then press Enter';
        end
        function obj = getPoints(obj)
            [obj.x obj.y] = getpts(gcf);
        end
        function [x y z] = coords2D(obj)
            loc = obj.locations;
            theta = linspace(0,2*pi,20)';
            x = zeros(0,1);
            y = zeros(0,1);
            z = zeros(0,1);
            elev = linspace(-1,1,7);
            elev([1 end])= [];
            thlist = linspace(0,2*pi,8);
            thlist(end) = [];
            for ndx = 1:size(loc,1)
                for el = elev
                    x = [x; loc(ndx,1)+obj.radius*sqrt(1-el^2)*cos(theta); NaN]; %#ok<AGROW>
                    y = [y; loc(ndx,2)+obj.radius*sqrt(1-el^2)*sin(theta); NaN]; %#ok<AGROW>
                    z = [z; obj.elevation+obj.radius*el*ones(size(theta)); NaN]; %#ok<AGROW>
                end
                elev2 = linspace(-1,1,length(theta))';
                for th = thlist
                    x = [x; obj.radius*sqrt(1-elev2.^2).*cos(th)+loc(ndx,1)]; %#ok<AGROW>
                    y = [y; obj.radius*sqrt(1-elev2.^2).*sin(th)+loc(ndx,2)]; %#ok<AGROW>
                    z = [z; obj.elevation+elev2*obj.radius]; %#ok<AGROW>
                end
                x(end+1) = NaN; %#ok<AGROW>
                y(end+1) = NaN; %#ok<AGROW>
                z(end+1) = NaN; %#ok<AGROW>
            end
        end
        function objSurface = coords3D(obj)
            texture = tile(obj.texture,obj.tiling);
            
            x = texture.triangles.vertices(:,1);
            y = texture.triangles.vertices(:,2);
            ynorm = (y-min(y(:)))/range(y(:));
            pt(:,1) = obj.radius*sqrt(1-(2*ynorm-1).^2).*cosd(360*(x-min(x(:)))./range(x(:)));
            pt(:,2) = obj.radius*sqrt(1-(2*ynorm-1).^2).*sind(360*(x-min(x(:)))./range(x(:)));
            
            objSurface.vertices = zeros(0,2);
            objSurface.triangulation = zeros(0,3);
            loc = obj.locations;
            for ndx = 1:size(loc,1)
                objSurface.triangulation = [objSurface.triangulation; texture.triangles.triangulation+size(objSurface.vertices,1)];
                objSurface.vertices = [objSurface.vertices; bsxfun(@plus,pt,loc(ndx,:))];
            end
            
            objSurface.vertices(:,3) = repmat((2*ynorm-1)*obj.radius+obj.elevation,size(loc,1),1);
            objSurface.cdata = repmat(texture.triangles.cdata,size(loc,1),1);
        end
        function edges = edges(obj)
            loc = obj.locations;
            edges = [loc loc];
        end
    end
end