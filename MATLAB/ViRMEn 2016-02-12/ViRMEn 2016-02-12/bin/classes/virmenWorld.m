classdef virmenWorld < virmenClass
    properties (SetObservable)
        backgroundColor = [0 0 0];
        startLocation = [0 0 0 0];
        transparency = 0;
        objects = {};
    end
    methods (Hidden = true)
        function [objSurface, objVertices, objTriangles] = coords3D(obj)
            objSurface.vertices = zeros(0,3);
            objSurface.cdata = zeros(0,4);
            objSurface.triangulation = zeros(0,3);
            objTriangles = zeros(length(obj.objects),1);
            objVertices = zeros(length(obj.objects),1);
            for ndx = 1:length(obj.objects)
                os = obj.objects{ndx}.coords3D;
                
                os.triangulation(all(isnan(os.cdata(os.triangulation(:,1),:)),2),:) = [];
                f = unique(os.triangulation(:));                
                vlist = zeros(size(os.vertices,1),1);
                vlist(f) = 1;
                numbefore = cumsum(vlist);
                nummissing = (1:size(os.vertices,1))'-numbefore;
                os.triangulation = os.triangulation - nummissing(os.triangulation);
                os.vertices = os.vertices(f,:);
                os.cdata = os.cdata(f,:);
                
                objTriangles(ndx,1) = size(objSurface.triangulation,1)+1;
                objSurface.triangulation = [objSurface.triangulation; os.triangulation+size(objSurface.vertices,1)];
                objTriangles(ndx,2) = size(objSurface.triangulation,1);
                objVertices(ndx,1) = size(objSurface.vertices,1)+1;
                objSurface.vertices = [objSurface.vertices; os.vertices];
                objVertices(ndx,2) = size(objSurface.vertices,1);
                objSurface.cdata = [objSurface.cdata; os.cdata];
            end
        end
    end
    methods
        function ch = children(obj)
            ch = obj.objects;
        end
        function h = draw3D(obj)
            [objSurface, objVertices, objTriangles] = obj.coords3D; %#ok<ASGLU>
            if isempty(objSurface.vertices)
                h = [];
            else
                h = trisurf(objSurface.triangulation,objSurface.vertices(:,1), ...
                    objSurface.vertices(:,2),objSurface.vertices(:,3), ...
                    zeros(size(objSurface.vertices,1),1));
                set(h,'facevertexcdata',objSurface.cdata(:,1:3),'facecolor','flat');
                set(h,'edgecolor','none');
                col = objSurface.cdata(:,4);
                if obj.transparency == 0
                    col(:) = 1-eps;
                else
                    col(isnan(col)) = 1-eps;
                end
                set(h,'alphadatamapping','none','facealpha','flat','facevertexalphadata',col);
            end
        end
        function [h, he, hp] = draw2D(obj)
            h = [];
            he = [];
            for ndx = 1:length(obj.objects)
                [h(ndx), he(ndx)] = obj.objects{ndx}.draw2D; %#ok<AGROW>
                hold on
            end
            loc = [0 1; 1 1; 1 -1; -1 -1; -1 1; 0 1; 0 4; 1 2; NaN NaN; -1 2; 0 4];
            rot = [cos(obj.startLocation(4)) -sin(obj.startLocation(4)); ...
                sin(obj.startLocation(4)) cos(obj.startLocation(4))];
            loc = rot*loc';
            hp = plot(loc(1,:)+obj.startLocation(1),loc(2,:)+obj.startLocation(2),'k','linewidth',2);
            hold on
        end
        function world = addObject(world,objectToAdd,varargin)
            if nargin == 2
                mode = 'copy';
            else
                mode = varargin{1};
            end
            if strcmp(mode,'copy')
                object = copyVirmenObject(objectToAdd);
            else
                object = objectToAdd;
            end
            object.parent = world;
            world.objects{end+1} = object;
        end
        function vrWorld = triangulate(world)
            vrWorld = loadVirmenWorld(world);
        end
    end
end