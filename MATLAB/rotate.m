function [rotated] = rotate(M,deg)

%%% Inputs:
%%%
%%% M - matrix with two rows and any number of columns. The top row should
%%% be your x values and the bottom row should be your y values
%%%
%%% deg - amount of roation in degres
%%%
%%% Returns:
%%%
%%% rotated - matrix the same size as M containing the rotated points
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% convert degrees to radians
rad = deg/360*2*pi;

% make rotation matrix
Q = [cos(-rad), -sin(-rad); sin(-rad) cos(-rad)];

% rotate points
rotated = Q*M;