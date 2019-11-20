function texture = image2texture(texture)

if texture.shapes{2}.interpolate
    img = imresize(texture.shapes{2}.image,[texture.shapes{2}.height texture.shapes{2}.width]);
else
    img = imresize(texture.shapes{2}.image,[texture.shapes{2}.height texture.shapes{2}.width],'nearest');
end
img = img(end:-1:1,:,:);

c = rgb2hsv(img);

v = c(:,:,3);
baseline = mean(v(:));
v = (v - baseline)*10^(texture.shapes{2}.contrast/100);
v = v + baseline + texture.shapes{2}.brightness/100;
c(:,:,3) = v;

c(:,:,2) = c(:,:,2) + texture.shapes{2}.saturation/100;

c(c<0) = 0;
c(c>1) = 1;

img = hsv2rgb(c);

h = size(img,1);
w = size(img,2);

[x, y] = meshgrid(linspace(0,1,w+1),linspace(0,1,h+1));
t.vertices = [x(:) y(:)];

% First vertex is the bottom-left corner of each triangle
v = (1:(h+1)*w)';
v(h+1:h+1:end)= [];
v = [v v]';
v = v(:);
t.triangulation(:,1) = v;

% Third vertex is the top-right corner of each triangle
t.triangulation(:,3) = v+h+2;

% Second vertex is the top for half of triangles, right for the other half
t.triangulation(1:2:end,2) = v(1:2:end)+1;
t.triangulation(2:2:end,2) = v(2:2:end)+h+1;

% Make all triangles counterclockwise
t.triangulation(1:2:end,[2 3]) = t.triangulation(1:2:end,[3 2]);

% 
t.cdata = NaN(length(x(:)),4);
img(end+1,:,:) = NaN;
img(:,end+1,:) = NaN;
img = img(:);
t.cdata(:,1:3) = reshape(img,(h+1)*(w+1),3);
t.cdata(:,4) = texture.shapes{2}.Alpha;

texture.triangles = t;