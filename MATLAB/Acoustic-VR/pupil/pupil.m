%pupil tracking
function [pupil_c, pupil_d, blinks,F] = pupil(frames)

dt=50; %analyze every dt'th frame
frames = frames(1:dt:end);

numframes = length(frames);
pupil_c = zeros(numframes,2); %pupil center
pupil_d = zeros(numframes,1); %pupil diameter
blinks = zeros(numframes,1);  %blink logical

th = 55;  %max pupil brightness
LB = 1500 ; %lower bound for pupil area

figure
ax = gca;
ax.NextPlot = 'replaceChildren';
F(length(frames)) = struct('cdata',[],'colormap',[]);

for f = 1:length(frames)
    
    frame = double(frames{f}); %get current frame
%     hy = fspecial('sobel');
%             hx = hy';
%             Iy = imfilter(double(frame), hy, 'replicate');
%             Ix = imfilter(double(frame), hx, 'replicate');
%             gradmag = sqrt(Ix.^2 + Iy.^2);
    %     if f == 2653
    imagesc(frame);colormap gray;axis xy;axis square;hold on
    %     end

    %binarize by global threshold
    tmpframe = frame;
    tmpframe(tmpframe>=th) = th;
    tmpframe(tmpframe<th) = 1;
    tmpframe(tmpframe~=1) = 0;
    
    if sum(sum(tmpframe)) == 0
        blinks(f) = 1; %detect blink
        %         imagesc(double(frames{f}));colormap gray;axis xy;axis square;keyboard
        drawnow
        F(f) = getframe(gcf);
    else
        %remove small pixel clusters
        tmpframe = bwareaopen(tmpframe,LB);
        
        if sum(sum(tmpframe)) == 0
            blinks(f) = 1; %detect blink
            %             imagesc(double(frames{f}));colormap gray;axis xy;axis square;pause
            drawnow
            F(f) = getframe(gcf);
        else
            %get cluster outlines
            houtlines = zeros(size(tmpframe));
            voutlines = zeros(size(tmpframe));
            for i=1:size(tmpframe,1)
                houtlines(i,2:end)=abs(diff(tmpframe(i,:)));
            end
            for i=1:size(tmpframe,2)
                voutlines(2:end,i)=abs(diff(tmpframe(:,i)));
            end
            outlines = houtlines | voutlines;
            
            %compute pupil_c, pupil_d
            [y,x] = find(outlines);
            dx = x-x';dy=y-y';
            dist = sqrt(dx.^2+dy.^2);
            pupil_c(f,1) = mean([min(x),max(x)]);
            pupil_c(f,2) = mean([min(y),max(y)]);
            pupil_d(f) = max(dist(:));
%             if f == 2653                
                plot(x,y,'.r');
                plot(pupil_c(f,1),pupil_c(f,2),'+g','MarkerSize',15);
                drawnow
                F(f) = getframe(gcf);
%             end
        end
    end
end