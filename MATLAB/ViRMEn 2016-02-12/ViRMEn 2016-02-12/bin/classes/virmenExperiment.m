classdef virmenExperiment < virmenClass
    properties (SetObservable)
        windows = {virmenWindow};
        movementFunction = @undefined;
        transformationFunction = @undefined;
        experimentCode = @undefined;
        codeText = {};
        worlds = {};
    end
    methods
        function obj = virmenExperiment
            w = virmenWorld;
            w.parent = obj;
            obj.worlds = {w};
        end
        function ch = children(obj)
            ch = obj.worlds;
        end
        function exper = addWorld(exper,worldToAdd,varargin)
            if nargin == 2
                mode = 'copy';
            else
                mode = varargin{1};
            end
            if strcmp(mode,'copy')
                world = copyVirmenObject(worldToAdd);
            else
                world = worldToAdd;
            end
            world.parent = exper;
            exper.worlds{end+1} = world;
        end
        function obj = updateCodeText(obj)
            obj = updateCode(obj);
        end
        function varargout = run(obj)
            err = virmenEngine(obj);
            if nargout == 1
                varargout{1} = err;
            end
        end
    end
end