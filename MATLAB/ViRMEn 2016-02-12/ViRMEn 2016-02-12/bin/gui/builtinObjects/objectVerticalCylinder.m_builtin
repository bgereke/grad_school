classdef objectVerticalCylinder < virmenObject
    properties (SetObservable)
        radius = 5;
        bottom = 0;
        top = 40;
        startAngle = 0;
        stopAngle = 360;
    end
    methods
        function obj = objectVerticalCylinder
            obj.iconLocations = [0 0];
            obj.helpString = 'Click cylinder centers, then press Enter';
        end
        function obj = getPoints(obj)
            [obj.x, obj.y] = getpts(gcf);
        end
        function [x, y, z] = coords2D(obj)
            st = obj.stopAngle;
            if st < obj.startAngle
                st = st+360;
            end
            theta = linspace(obj.startAngle,st,20)';
            xs = obj.radius*cosd(theta);
            ys = obj.radius*sind(theta);
            x = zeros(0,1);
            y = zeros(0,1);
            z = zeros(0,1);
            zlist = linspace(obj.bottom,obj.top,5);
            thlist = linspace(obj.startAngle,st,8);
            loc = obj.locations;
            for zs = zlist
                for ndx = 1:size(loc,1)
                    x = [x; xs+loc(ndx,1); NaN]; %#ok<AGROW>
                    y = [y; ys+loc(ndx,2); NaN]; %#ok<AGROW>
                    z = [z; zs*ones(size(xs)); NaN]; %#ok<AGROW>
                    for th = thlist
                        x = [x; obj.radius*cosd([th; th])+loc(ndx,1); NaN]; %#ok<AGROW>
                        y = [y; obj.radius*sind([th; th])+loc(ndx,2); NaN]; %#ok<AGROW>
                        z = [z; obj.bottom; obj.top; NaN]; %#ok<AGROW>
                    end
                end
            end
            x(end) = [];
            y(end) = [];
            z(end) = [];
        end
        function objSurface = coords3D(obj)
            texture = tile(obj.texture,obj.tiling);
            
            x = texture.triangles.vertices(:,1);
            angle_range = mod(obj.stopAngle-obj.startAngle,360);
            if angle_range == 0
                angle_range = 360;
            end
            pt(:,1) = obj.radius*cosd(angle_range*(x-min(x(:)))./range(x(:))+obj.startAngle);
            pt(:,2) = obj.radius*sind(angle_range*(x-min(x(:)))./range(x(:))+obj.startAngle);
            
            objSurface.vertices = zeros(0,2);
            objSurface.triangulation = zeros(0,3);
            loc = obj.locations;
            for ndx = 1:size(loc,1)
                objSurface.triangulation = [objSurface.triangulation; texture.triangles.triangulation+size(objSurface.vertices,1)];
                objSurface.vertices = [objSurface.vertices; bsxfun(@plus,pt,loc(ndx,:))];
            end
            
            y = texture.triangles.vertices(:,2);
            objSurface.vertices(:,3) = repmat((y-min(y(:)))/range(y(:))*(obj.top-obj.bottom)+obj.bottom,size(obj.locations,1),1);
            objSurface.cdata = repmat(texture.triangles.cdata,size(obj.locations,1),1);
        end
        function edges = edges(obj)
            edges = [obj.locations obj.locations];
        end
    end
end