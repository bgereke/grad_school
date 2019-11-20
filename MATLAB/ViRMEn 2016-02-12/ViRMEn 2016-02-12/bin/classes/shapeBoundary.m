classdef shapeBoundary < virmenShape
    properties
    end
    methods
        function [x y] = coords2D(obj)
            loc = obj.locations;
            x = loc([1 1 2 2 1]',1);
            y = loc([1 2 2 1 1]',2);
        end
    end
end