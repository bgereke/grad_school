classdef virmenTexture < virmenClass
    properties (SetObservable)
        width = 1;
        height = 1;
        tilable = [1 1];
        refining = [1 1];
        grid = [1 1];
        shapes = {};
        triangles = struct;
    end
    methods (Hidden = true)
        function [points, segments] = getSegments(obj)
            [points, segments] = getTextureSegments(obj);
        end
        function tiled = tile(obj,tiling)
            tiled = tileTexture(obj,tiling);
        end
    end
    methods
        function obj = virmenTexture
            obj.triangles.vertices = [0 0; 0 1; 1 0; 1 1];
            obj.triangles.triangulation = [4 2 3; 2 1 3];
            obj.triangles.cdata = NaN(4,4);
            rect = shapeBoundary;
            rect.parent = obj;
            rect.x = [0 obj.width];
            rect.y = [0 obj.height];
            obj = addShape(obj,rect,'new');
        end
        function ch = children(obj)
            ch = obj.shapes;
        end
        function h = draw(obj)
            h = trisurf(obj.triangles.triangulation,obj.triangles.vertices(:,1),obj.triangles.vertices(:,2),zeros(size(obj.triangles.vertices,1),1));
            set(h,'facevertexcdata',obj.triangles.cdata(:,1:3),'facecolor','flat');
            col = obj.triangles.cdata(:,4);
            col(isnan(col)) = 1-eps;
            set(h,'alphadatamapping','none','facealpha','flat','facevertexalphadata',col);
        end
        function h = sketch(obj)
            h = zeros(1,length(obj.shapes));
            for ndx = 1:length(obj.shapes)
                hold on
                [x, y] = obj.shapes{ndx}.coords2D;
                if strcmp(class(obj.shapes{ndx}),'shapeColor')
                    h(ndx) = plot(x,y,'o','markerfacecolor',obj.shapes{ndx}.RGBA(1:3));
                    set(h(ndx),'color',[.5 .5 .5]);
                else
                    h(ndx) = plot(x,y);
                    set(h(ndx),'color','k','linewidth',1.5);
                end
            end
        end
        function obj = compute(obj)
            if length(obj.shapes)==2 & strcmp(class(obj.shapes{2}),'shapeImage') %#ok<AND2>
                obj = image2texture(obj);
            else
                obj = computeTexture(obj);
            end
        end
        function texture = addShape(texture,shapeToAdd,varargin)
            if nargin == 2
                mode = 'copy';
            else
                mode = varargin{1};
            end
            if strcmp(mode,'copy')
                shape = copyVirmenObject(shapeToAdd);
            else
                shape = shapeToAdd;
            end
            
            shape.parent = texture;
            texture.shapes{end+1} = shape;
            [points, segments] = getTextureSegments(texture); %#ok<NASGU>
            loc = shape.locations;
            if ~strcmp(class(shape),'shapeColor')
                for ndx = 1:size(shape.locations,1)
                    dst = sum(bsxfun(@minus,points,loc(ndx,:)).^2,2);
                    [dummy, f] = min(dst); %#ok<ASGLU>
                    f = f(1);
                    texture.shapes{end}.x(ndx) = points(f,1);
                    texture.shapes{end}.y(ndx) = points(f,2);
                end
            end
        end
        function texture = loadImage(texture,img)
            texture.shapes(2:end) = [];
            texture.addShape(shapeImage);
            texture.shapes{2}.image = img;
            texture.shapes{2}.height = size(img,1);
            texture.shapes{2}.width = size(img,2);
        end
    end
end