classdef virmenWorldPrimitives
    properties
        surface = struct;
        objects = struct;
        edges = struct;
        backgroundColor = [0 0 0]
        startLocation = [0 0 0 0]
    end
    properties (Hidden = true)
        walls = struct;
    end
end