function [ stack_mean,each_frame_mean ] = image_mean( image_stack)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

num_frames=length(image_stack);
for n=1:num_frames;
    image_all(:,:,n)=image_stack{n};
    each_frame_mean(n)=mean(mean(image_stack{n}));
end;
stack_mean = squeeze(mean(image_all,3));
%clear cam0_all % help memory
cam0_track_left

end

