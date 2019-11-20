classdef objectVerticalCone < virmenObject
    properties (SetObservable)
        bottomRadius = 5;
        topRadius = 0;
        bottom = 0;
        top = 10;
    end
    methods
        function obj = objectVerticalCone
            obj.iconLocations = [0 0];
            obj.helpString = 'Click cone centers, then press Enter';
        end
        function obj = getPoints(obj)
            [obj.x obj.y] = getpts(gcf);
        end
        function [x y z] = coords2D(obj)
            theta = linspace(0,2*pi,20)';
            x = zeros(0,1);
            y = zeros(0,1);
            z = zeros(0,1);
            zlist = linspace(obj.bottom,obj.top,5);
            thlist = linspace(0,2*pi,8);
            thlist(end) = [];
            loc = obj.locations;
            for zs = zlist
                for ndx = 1:size(loc,1)
                    x = [x; ((zs-obj.bottom)/(obj.top-obj.bottom)*(obj.topRadius-obj.bottomRadius)+obj.bottomRadius)*cos(theta)+loc(ndx,1); NaN]; %#ok<AGROW>
                    y = [y; ((zs-obj.bottom)/(obj.top-obj.bottom)*(obj.topRadius-obj.bottomRadius)+obj.bottomRadius)*sin(theta)+loc(ndx,2); NaN]; %#ok<AGROW>
                    z = [z; zs*ones(size(theta)); NaN]; %#ok<AGROW>
                    for th = thlist
                        x = [x; obj.bottomRadius*cos(th)+loc(ndx,1); obj.topRadius*cos(th)+loc(ndx,1); NaN]; %#ok<AGROW>
                        y = [y; obj.bottomRadius*sin(th)+loc(ndx,2); obj.topRadius*sin(th)+loc(ndx,2); NaN]; %#ok<AGROW>
                        z = [z; obj.bottom; obj.top; NaN]; %#ok<AGROW>
                    end
                end
            end
        end
        function objSurface = coords3D(obj)
            texture = tile(obj.texture,obj.tiling);
            
            x = texture.triangles.vertices(:,1);
            y = texture.triangles.vertices(:,2);
            ynorm = (y-min(y(:)))/range(y(:));
            pt(:,1) = (obj.bottomRadius+(obj.topRadius-obj.bottomRadius)*ynorm).*cosd(360*(x-min(x(:)))./range(x(:)));
            pt(:,2) = (obj.bottomRadius+(obj.topRadius-obj.bottomRadius)*ynorm).*sind(360*(x-min(x(:)))./range(x(:)));
            
            objSurface.vertices = zeros(0,2);
            objSurface.triangulation = zeros(0,3);
            loc = obj.locations;
            for ndx = 1:size(loc,1)
                objSurface.triangulation = [objSurface.triangulation; texture.triangles.triangulation+size(objSurface.vertices,1)];
                objSurface.vertices = [objSurface.vertices; bsxfun(@plus,pt,loc(ndx,:))];
            end
            
            objSurface.vertices(:,3) = repmat(ynorm*(obj.top-obj.bottom)+obj.bottom,size(loc,1),1);
            objSurface.cdata = repmat(texture.triangles.cdata,size(loc,1),1);
        end
        function edges = edges(obj)
            loc = obj.locations;
            edges = [loc loc];
        end
    end
end