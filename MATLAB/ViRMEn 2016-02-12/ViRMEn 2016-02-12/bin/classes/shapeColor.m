classdef shapeColor < virmenShape
    properties (SetObservable)
        R = [];
        G = [];
        B = [];
        Alpha = NaN;
    end
    methods
        function obj = shapeColor
            obj.helpString = 'Click color locations, then press Enter';
        end
        function obj = getPoints(obj)
            [obj.x obj.y] = getpts(gcf);
            if length(obj.RGBA)==1
                RGB = uisetcolor([1 1 1]);
                obj.R = RGB(1);
                obj.G = RGB(2);
                obj.B = RGB(3);
            end
        end
        function [x y] = coords2D(obj)
            loc = obj.locations;
            x = loc(:,1);
            x = [x'; NaN(size(x'))];
            x = x(:);
            y = loc(:,2);
            y = [y'; NaN(size(y'))];
            y = y(:);
        end
        function col = RGBA(varargin)
            obj = varargin{1};
            col = [obj.R obj.G obj.B obj.Alpha];
            if nargin>1
                col = col(varargin{2});
            end
        end
    end
end