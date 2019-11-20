classdef objectVerticalWall < virmenObject
    properties (SetObservable)
        bottom = 0;
        top = 10;
    end
    methods
        function obj = objectVerticalWall
            obj.iconLocations = [0 -5; 10 -10; 15 -20];
            obj.helpString = 'Click wall endpoints, then press Enter';
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
            x = [loc(:,1); loc(end:-1:1,1); NaN];
            y = [loc(:,2); loc(end:-1:1,2); NaN];
            z = [obj.bottom*ones(size(loc,1),1); obj.top*ones(size(loc,1),1); NaN];
            for ndx = 1:size(loc,1)-1
                x = [x; loc([ndx+1; ndx; ndx; ndx+1],1); NaN]; %#ok<AGROW>
                y = [y; loc([ndx+1; ndx; ndx; ndx+1],2); NaN]; %#ok<AGROW>
                z = [z; obj.bottom; obj.top; obj.bottom; obj.top; NaN]; %#ok<AGROW>
            end
            x(end) = [];
            y(end) = [];
            z(end) = [];
        end
        function objSurface = coords3D(obj)
            texture = tile(obj.texture,obj.tiling);
            
            x = texture.triangles.vertices(:,1);
            y = texture.triangles.vertices(:,2);
            x_norm = (x-min(x(:)))/range(x(:));
            loc = obj.locations;
            dst = [0; sqrt(sum((loc(2:end,:)-loc(1:end-1,:)).^2,2))];
            pos = cumsum(dst);
            pos = pos/pos(end);
            X = zeros(size(x));
            Y = zeros(size(x));
            loc = obj.locations;
            for ndx = 1:length(pos)-1
                f = find(x_norm(:)>=pos(ndx));
                X(f) = loc(ndx,1)+(x_norm(f)-pos(ndx))/(pos(ndx+1)-pos(ndx))*(loc(ndx+1,1)-loc(ndx,1));
                Y(f) = loc(ndx,2)+(x_norm(f)-pos(ndx))/(pos(ndx+1)-pos(ndx))*(loc(ndx+1,2)-loc(ndx,2));
            end
            Z = (y-min(y(:)))/range(y(:))*(obj.top-obj.bottom)+obj.bottom;
            objSurface.vertices = [X Y Z];
            objSurface.cdata = texture.triangles.cdata;
            objSurface.triangulation = texture.triangles.triangulation;
        end
        function edges = edges(obj)
            loc = obj.locations;
            edges = [loc(1:end-1,:) loc(2:end,:)];
        end
    end
end