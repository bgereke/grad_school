classdef virmenShape < virmenClass
    properties (SetObservable)
        x = zeros(0,1);
        y = zeros(0,1);
    end
    properties
        iconLocations = [0 0; 100 100];
        helpString = '';
    end
    methods
        function ch = children(obj) %#ok<MANU>
            ch = {};
        end
        function loc = locations(varargin)
            obj = varargin{1};
            if length(obj.x)==1 && length(obj.y)>1
                loc = [repmat(obj.x(:),length(obj.y),1) obj.y(:)];
            elseif length(obj.x)>1 && length(obj.y)==1
                loc = [obj.x(:) repmat(obj.y(:),length(obj.x),1)];
            else
                loc = [obj.x(:) obj.y(:)];
            end
            if nargin == 2
                loc = loc(varargin{2});
            elseif nargin == 3
                loc = loc(varargin{2},varargin{3});
            end
        end
    end
end