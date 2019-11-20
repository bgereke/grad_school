classdef shapeImage < virmenShape
    properties (SetObservable)
        height
        width
        interpolate = 1;
        brightness = 0;
        contrast = 0;
        saturation = 0;
        Alpha = NaN;
    end
    properties (Hidden)
        image = [];
    end
    methods
        function obj = shapeImage(obj) %#ok<INUSD>
            obj.x = 0.5;
            obj.y = 0.5;
        end
        function [x, y] = coords2D(obj) %#ok<MANU>
        crd = [0.175 0.197; 0.359 0.545; 0.438 0.401; 0.582 0.655; 0.825 0.196; 0.175 0.197];
        th = linspace(0,2*pi,100)';
        crd(end+1,:) = NaN;
        crd = [crd; .17+.08*cos(th) .72-.103*sin(th)];
        crd(end+1,:) = NaN;
        crd(:,1) = crd(:,1)+0.05;
        crd = [crd; 0 0; 0 1; 1 1; 1 0; 0 0];
        crd(:,2) = crd(:,2)*0.776;
        
        crd = (crd-0.5)*0.25 + 0.5;
        
        x = crd(:,1);
        y = crd(:,2);
        end
    end
end