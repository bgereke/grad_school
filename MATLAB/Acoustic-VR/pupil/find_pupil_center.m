function [ cam0_track_left,cam0_track_circle,pupil_data ] = find_pupil_center( cam0,P )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% initialize
warning('off','all');

[frame_width, frame_height]=size(cam0{1});
num_frames=length(cam0);

% intialaize matrices
failed_frames=zeros(1,num_frames);
 c_x_left=zeros(1,num_frames);
 c_y_left=zeros(1,num_frames);
 pupil_r_left=zeros(1,num_frames); 
 CR_x_left=zeros(1,num_frames);
 CR_y_left=zeros(1,num_frames);
 cam0_track_left = cell(1,num_frames);

v2struct(P);

%% Image mean and 
 [cam0_mean,~]=image_mean(cam0);

%%
%average corneal reflection detection
[centers, radii] = imfindcircles(cam0_mean,[5 15],'ObjectPolarity','bright','Sensitivity',.95);
CR_x_left_mean = centers(1,1);
CR_y_left_mean = centers(1,2);


%%
for m=1:num_frames;
    
    disp(m);
    
    
    %filtering PLAY WITH THESE NUMBERS (VIEW A VIDEO OF EACH TEST IMAGE)
    testimage = cam0{m};
    
    testimage2=medfilt2(testimage,[7 7]); %filtering out reflections
    testimage3 = testimage2;
    testimage3(find(testimage2<=min(min(testimage2(range_for_min{1},range_for_min{2})))+17))=0; %making the pupil true BLACK!
    %MAKE SURE the step above makes ALL of the pupil black and hopefully
    %NONE of the shadow. You increase from 25 to include more pixels for the pupil and
    %decrease it, if it includes too many pixels of the shadows. You want
    %to keep the shadow and pupil from ever overlapping.
    
    
    
    %instantaneous pupil size and location detection
    if m > 1
        [centers, radii] = imfindcircles(testimage3,[round(pupil_r_left(m-1)-.25*pupil_r_left(m-1)) round(pupil_r_left(m-1)+.25*pupil_r_left(m-1))],'ObjectPolarity','dark','Sensitivity',sensitivity);
    else
        [centers, radii] = imfindcircles(testimage3,[round(initial_r-.25*initial_r) round(initial_r+.25*initial_r)],'ObjectPolarity','dark','Sensitivity',sensitivity);
    end;
    
    
if (~ismember(m,blink_ind))
    if size(centers,2) > 0 & m > 1;
        %use the similar pupil location as the one detected in the
        %previous image
        [j,k] = min((centers(:,1)-c_x_left(m-1)).^2+(centers(:,2)-c_y_left(m-1)).^2);
        %size change too big = not the pupil: use previous location
        if abs((radii(k) - pupil_r_left(m-1))/pupil_r_left(m-1)) > 0.25 
            %keep trying to find the real pupil
            [centers, radii] = imfindcircles(testimage3,[round(pupil_r_left(m-1)-.25*pupil_r_left(m-1)) round(pupil_r_left(m-1)+.25*pupil_r_left(m-1))],'ObjectPolarity','dark','Sensitivity',sensitivity+.01);
            if size(centers,2) > 0
                [j,k] = min((centers(:,1)-c_x_left(m-1)).^2+(centers(:,2)-c_y_left(m-1)).^2);
                %size change too big = not the pupil: use previous location
                if abs((radii(k) - pupil_r_left(m-1))/pupil_r_left(m-1)) > 0.25
                    %give up
                    c_x_left(m) = c_x_left(m-1);
                    c_y_left(m) = c_y_left(m-1);
                    pupil_r_left(m) = pupil_r_left(m-1);
                else
                    if centers(k,1) > x_bound(1) & centers(k,1) < x_bound(2) & centers(k,2) > y_bound(1) & centers(k,2) < y_bound(2) %based on 250x250 image
                        %yes!
                        c_x_left(m) = centers(k,1);
                        c_y_left(m) = centers(k,2);
                        pupil_r_left(m) = radii(k);
                    else
                        %give up
                        c_x_left(m) = c_x_left(m-1);
                        c_y_left(m) = c_y_left(m-1);
                        pupil_r_left(m) = pupil_r_left(m-1);
                    end
                end;
            else
                %give up
                c_x_left(m) = c_x_left(m-1);
                c_y_left(m) = c_y_left(m-1);
                pupil_r_left(m) = pupil_r_left(m-1);
            end;
        else
            if centers(k,1) > x_bound(1) & centers(k,1) < x_bound(2) & centers(k,2) > y_bound(1) & centers(k,2) < y_bound(2) %based on 250x250 image
                %yes!
                c_x_left(m) = centers(k,1);
                c_y_left(m) = centers(k,2);
                pupil_r_left(m) = radii(k);
            else
                %give up
                c_x_left(m) = c_x_left(m-1);
                c_y_left(m) = c_y_left(m-1);
                pupil_r_left(m) = pupil_r_left(m-1);
            end
        end;
    elseif m > 1;  %no circle found (whisker noise): use previous location
        c_x_left(m) = c_x_left(m-1);
        c_y_left(m) = c_y_left(m-1);
        pupil_r_left(m) = pupil_r_left(m-1);
        failed_frames(m)=1;
    else
        c_x_left(m) = initial_x;
        c_y_left(m) = initial_y;
        pupil_r_left(m) = initial_r;
    end;
    
else
    % blink, use previous frame
    c_x_left(m) = c_x_left(m-1);
        c_y_left(m) = c_y_left(m-1);
        pupil_r_left(m) = pupil_r_left(m-1);
        failed_frames(m)=1;
end
    
    %instantaneous corneal reflection detection
    [centers, radii] = imfindcircles(cam0{m},[5 15],'ObjectPolarity','bright','Sensitivity',.95);
    
    %use the same corneal reflection as the first one detected in the
    %average image
    [j,k] = min((centers(:,1)-CR_x_left_mean).^2+(centers(:,2)-CR_y_left_mean).^2);
    CR_x_left(m) = centers(k,1);
    CR_y_left(m) = centers(k,2);
    
    %generate movie with crosshairs on pupil and corneal reflection
    testimage_track=testimage;
    c_x = round(c_x_left(m));
    c_y = round(c_y_left(m));
    %pupil white crosshairs
    testimage_track(c_y-10:c_y+10,c_x-1:c_x+1)=255;
    testimage_track(c_y-1:c_y+1,c_x-10:c_x+10)=255;
    
     % add circle around pupil - mjd 
    testimage_circle=testimage;
    thick=.1; % circle thickness ( % of radius)
    [rr cc] = meshgrid(1:frame_width,1:frame_height);
    C_outer = sqrt((rr-c_x).^2+(cc-c_y).^2)<=(round(pupil_r_left(m)+pupil_r_left(m)*thick));
    C_inner = sqrt((rr-c_x).^2+(cc-c_y).^2)<=pupil_r_left(m);
    C=C_outer-C_inner;
    [row,col]=find(C);
    for jj=1:length(row);testimage_circle(row(jj),col(jj))=255;end;
    cam0_track_circle{m}=testimage_circle;
    
    
    
    c_x = round(CR_x_left(m));
    c_y = round(CR_y_left(m));
    %corneal reflection black crosshairs
    testimage_track(c_y-5:c_y+5,c_x-1:c_x+1)=0;
    testimage_track(c_y-1:c_y+1,c_x-5:c_x+5)=0;
    cam0_track_left{m}=testimage_track; %WATCH A MOVIE OF THIS!!! THIS IS THE RESULTS/OUTPUT

   
    
end;

c_x_left = c_x_left - mean(c_x_left);   %PUPIL LOCATION (don't have to subtract mean)
c_y_left = c_y_left - mean(c_y_left);
CR_x_left = CR_x_left - CR_x_left_mean; %IR reflection LOCATION (don't have to subtract mean)
CR_y_left = CR_y_left - CR_y_left_mean;


pupil_data=v2struct(c_x_left,c_y_left,CR_x_left,CR_y_left,pupil_r_left,failed_frames);

