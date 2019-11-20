classdef virmenObject < virmenClass
    properties (SetObservable)
        x = zeros(0,1);
        y = zeros(0,1);
        tiling = [5 5];
        edgeRadius = NaN;
        texture = virmenTexture;
    end
    properties
        iconLocations = [0 0; 100 100];
        helpString = '';
    end
    methods
        function obj = virmenObject
            t = virmenTexture;
            obj = setTexture(obj,t,'new');
        end
        function ch = children(obj)
            ch = {obj.texture};
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
        function h = draw3D(obj)
            objSurface = obj.coords3D;
            h = trisurf(objSurface.triangulation,objSurface.vertices(:,1), ...
                objSurface.vertices(:,2),objSurface.vertices(:,3), ...
                zeros(size(objSurface.vertices,1),1));
            set(h,'facevertexcdata',objSurface.cdata(:,1:3),'facecolor','flat');
            set(h,'edgecolor','none');
            col = objSurface.cdata(:,4);
            col(isnan(col)) = 1-eps;
            set(h,'alphadatamapping','none','facealpha','flat','facevertexalphadata',col);
        end
        function [h he] = draw2D(obj)
            [x, y, z] = obj.coords2D; %#ok<*PROP>
            if isempty(x)
                h = plot(NaN,NaN,'k');
            else
                h = plot3(x,y,z,'k');
            end
            hold on
            edges = obj.edges;
            xs = zeros(1,0);
            ys = zeros(1,0);
            for ndx = 1:size(edges,1)
                x = edges(ndx,1)+obj.edgeRadius*cos(linspace(0,2*pi,100));
                y = edges(ndx,2)+obj.edgeRadius*sin(linspace(0,2*pi,100));
                x = [x NaN edges(ndx,3)+obj.edgeRadius*cos(linspace(0,2*pi,100))]; %#ok<AGROW>
                y = [y NaN edges(ndx,4)+obj.edgeRadius*sin(linspace(0,2*pi,100))]; %#ok<AGROW>
                ang = atan2(edges(ndx,4)-edges(ndx,2),edges(ndx,3)-edges(ndx,1));
                ang = ang+pi/2;
                x = [x NaN edges(ndx,[1 3])+obj.edgeRadius*cos(ang)]; %#ok<AGROW>
                y = [y NaN edges(ndx,[2 4])+obj.edgeRadius*sin(ang)]; %#ok<AGROW>
                x = [x NaN edges(ndx,[1 3])-obj.edgeRadius*cos(ang)]; %#ok<AGROW>
                y = [y NaN edges(ndx,[2 4])-obj.edgeRadius*sin(ang)]; %#ok<AGROW>
                xs = [xs x NaN]; %#ok<AGROW>
                ys = [ys y NaN]; %#ok<AGROW>
            end
            if ~isempty(xs)
                xs(end) = [];
                ys(end) = [];
            end
            if isempty(xs)
                he = plot(NaN,NaN,':k');
            else
                he = plot(xs,ys,':k');
            end
        end
        function obj = setTexture(obj,textureToAdd,varargin)
            if nargin == 2
                mode = 'copy';
            else
                mode = varargin{1};
            end
            if strcmp(mode,'copy')
                newTexture = copyVirmenObject(textureToAdd);
            else
                newTexture = textureToAdd;
            end
            newTexture.parent = obj;
            obj.texture = newTexture;
        end
    end
end