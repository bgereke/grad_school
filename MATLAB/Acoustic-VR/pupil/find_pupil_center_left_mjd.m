%% Files method 1: direct read jpegs (slow)

% NOTE: can probably just generate the strings since they are the same for
% every file

%path='C:\\Users\\Matt\\Desktop\\test_eye\\';
path='Z:\eye_vids\frames_214_01262017\';

%%first_frame='0000.jpg';

%file=strcat(path, first_frame);

files = dir(fullfile(path, '*.jpg'));
files = {files.name}';

[a,b,c]=cellfun(@(x) fileparts(x), files, 'UniformOutput', false);
[~,b_sort]=sort(str2double(b));
files=files(b_sort);

%% Files method 2: generate file names 

% path='Z:\eye_vids\frames_md218_01272017\';
num_files=36300;
%num_files=30;

ds=1;
N=1:ds:num_files;
for ii=1:length(N)
    % pad with up to 4 zeros
    if (ii <10000)
        Nstr=num2str(N(ii),'%04d');
        files{ii,1}=strcat(path,Nstr,'.jpg');
    else
        Nstr=num2str(N(ii));
        files{ii,1}=strcat(path,Nstr,'.jpg');
    end
end
    

%%
% Read in images
% timing: 230 secs to load 18150 (60hz-->30hz); 96-151 sec, 60-->10hz)
tic
ds=1;
% ds=1; % 60--> 30hz
select=1:ds:length(files);

for i=1:length(select)
    cam0{i}=imread(files{select(i)});
%     cam0{i}=imread(files{i});
end
cam0_full=cam0;
%%
% select smaller portion

cam0=cam0_full;
toc
%% set to mc image
for i=1:length(Istack_warped)
    im=Istack_warped{i};
    im(1:15,:)=[];
    im(:,1:15)=[];
    im(end-15:end,:)=[];
    im(:,end-15:end)=[];
    img_cl{i}=im;
end
cam0=img_cl;

%% load first frame, estimate center and radius
[frame_width, frame_height]=size(cam0{1});
num_frames=length(cam0);

figure;
imshow(cam0{1});
[P.initial_x,P.initial_y]=ginput(1);
[x,y]=ginput(2);

P.initial_r=sqrt((x(2)-x(1))^2+(y(2)-y(1))^2)/2; %radius

%% Mean
[cam0_mean,whole_frame_mean]=image_mean(cam0);
figure; imshow(uint8(cam0_mean)); % show mean

%% BLINK DETECTION
% NEED TO IMPROVE
%if baseline fluctations, the median may not be a good estimate

wff=imgaussfilt(whole_frame_mean,3);% whole frame filt
wff_nobl=wff; % no blinks
blink_ind=[];

for ii=1:2
    wff_nobl=zscore(wff_nobl);
    med_wf=median(wff_nobl);
    thresh=std(wff_nobl)*3+med_wf;
    ind=find(abs(wff_nobl) > thresh);
    blink_ind=[blink_ind ind];
    wff_nobl(ind)=med_wf;
end


figure;
subplot(211);
plot(wff);
hold on;
plot(f0,'r');

subplot(212);
%wff(blink_ind)=NaN;
plot(wff_nobl);


%%

wff=imgaussfilt(whole_frame_mean,3);% whole frame filt
wff_nobl=wff; % no blinks


imgHz=10;

     wSize= 10; % s, window size
        winFrames= round(100*wSize*1/imgHz); % # frames to grab
    
          f0=[];
        for i=1:length(wff)
            if (i<winFrames+1)
                s=i+1; e=i+winFrames;
                v=wff(1,s:e);
                sv=sort(v);
                f0(i)=mean(sv(1,1:50));

               
            else
                s=i-winFrames; e=i-1;
                v=wff(1,s:e);
                sv=sort(v);
                f0(i)=mean(sv(1,1:50));

               
            end
        end
%%
blink_ind=[];
wff_nobl=wff;
for ii=1:length(wff)
    %wff_nobl=zscore(wff_nobl);
    med_wf=f0(ii);
    thresh=std(wff)*2+med_wf;
    if (abs(wff(ii)) > thresh)
        blink_ind=[blink_ind ii];
        wff_nobl(ii)=med_wf;
    end
end


figure;
subplot(211);
plot(wff);
hold on;
plot(f0,'r');

subplot(212);
%wff(blink_ind)=NaN;
plot(wff_nobl);

%%
pupil_gui(cam0);
%%  MAIN CODE

tic;
% PARAMETERS
P.range_for_min={50:183,20:155}; %range to look for pixels values, (1)=x, (2)=y;  250x250 = {50:200,50:200} 

xb=round(frame_width*.15); 
P.x_bound=[xb frame_width-xb];% bounds on where to look for pupil center, standard approach is 4-5% of dimenesion size. 250x250 ~ [11, 239]; 
yb=round(frame_height*.15); 
P.y_bound=[yb frame_height-yb];

P.blink_ind=blink_ind;

P.sensitivity = 0.95;   %typically use 0.95, but increase if pupil is not dark enough
% END PARAMETERS

[pupil_movie, pupil_movie_circle, PD] = find_pupil_center(cam0,P);

TimeSpent = toc;




%% PLOT PUPIL
figure; 

subplot(121);
plot(PD.pupil_r_left);
title('Pupil Size');
subplot(122);
plot(Pd.c_x_left);
title('Pupil Center - X');
per_failed=length(find(failed_frames))/num_frames;

%%
pupil_gui(pupil_movie_circle)

%% ASSIGN VARS 2 STRUCT
pupil_size=pupil_r_left;
center_x=c_x_left;
center_y=c_y_left;

info=struct;
pupil=struct;

PUPIL=v2struct(pupil_size,center_x,center_y,blink_ind,whole_frame_filt,info);

%% SCATTER PLOT
figure;
scatter(c_x_left,pupil_r_left);


%% Create circle in image
radius=30
thick=.1 %  percent larger to make circle
center_x=100;
center_y=100;

figure;
[rr cc] = meshgrid(1:frame_width,1:frame_height);
C_outer = sqrt((rr-center_x).^2+(cc-center_y).^2)<=(round(radius+radius*thick));
C_inner = sqrt((rr-center_x).^2+(cc-center_y).^2)<=radius;
C=C_outer-C_inner;
imshow(C);

%% Alternative Blink detection
% find peaks of whole frame fluro
blink_width=3; % NOTE: depends on sampling rate
whole_frame_filt=imgaussfilt(whole_frame_mean,3);

%figure;findpeaks(whole_frame_filt,'MinPeakHeight',97,...
%    'MinPeakWidth',blink_width,...
%    'Annotate','extents');
%figure;plot(whole_frame_filt);

% run through twice to get rid of huge fluctations first

%% ************************************************************************
%% MOVIE STUFF
%% ************************************************************************

%% turn images to movie
clear M;
cmap = gray(256);
for ii=1:num_frames
    %M(ii)=im2frame(cam0_track_left{ii},cmap); %movie w/ tracking
    %M(ii)=im2frame(cam0{ii},cmap);
    
    M(ii)=im2frame(pupil_movie_circle{ii},cmap);
end

%figure;movie(M);

%% MOVIE TO AVI
movie_path='C:\\Users\\Matt\\Desktop\\test_eye\\'

v = VideoWriter(strcat(movie_path,'md218_10hz.avi'),'Motion JPEG AVI');
%v.CompressionRatio = 3;
open(v)
writeVideo(v,rM);
close(v);


%% FIND WRONG SIZED FRAMES IN MOVIE
% needed to write the avi
rM=M;
rind=[];
for ii=1:num_frames
    if(size(M(ii).cdata,2) ~=168);display(ii); rind=[rind ii];;
    end;
    
end
rM(rind)=[];

