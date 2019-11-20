classdef shapeRegularPolygon < virmenShape
    properties (SetObservable)
        points = 8;
        radius = .1;
        rotation = 0;
    end
    methods
        function obj = shapeRegularPolygon
            obj.iconLocations = [0 0; .3 .3];
            obj.helpString = 'Click polygon centers, then press Enter';
        end
        function obj = getPoints(obj)
            [x y] = getpts(gcf);
            obj.x = x;
            obj.y = y;
        end
        function [x y] = coords2D(obj)
            theta = linspace(0,360,obj.points+1)+obj.rotation;
            
            xc = obj.radius*cosd(theta(:));
            yc = obj.radius*sind(theta(:));
            x = zeros(0,1);
            y = zeros(0,1);
            loc = obj.locations;
            for ndx = 1:size(obj.locations,1)
                x = [x; xc+loc(ndx,1); NaN]; %#ok<AGROW>
                y = [y; yc+loc(ndx,2); NaN]; %#ok<AGROW>
            end
            x(end) = [];
            y(end) = [];
        end
    end
end