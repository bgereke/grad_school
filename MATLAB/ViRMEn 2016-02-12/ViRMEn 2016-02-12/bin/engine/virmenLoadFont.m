function [letterGrid letterFont letterAspectRatio] = virmenLoadFont

letterGrid = [ ...
    0.0 1.0 1.0 1.0; ...
    0.0 1.0 0.0 0.5; ...
    0.0 1.0 0.5 0.5; ...
    0.5 1.0 0.5 0.5; ...
    1.0 1.0 0.5 0.5; ...
    1.0 1.0 1.0 0.5; ...
    0.0 0.5 0.5 0.5; ...
    0.5 0.5 1.0 0.5; ...
    0.0 0.0 0.0 0.5; ...
    0.0 0.0 0.5 0.5; ...
    0.5 0.0 0.5 0.5; ...
    1.0 0.0 0.5 0.5; ...
    1.0 0.0 1.0 0.5; ...
    0.0 0.0 1.0 0.0; ...
    0.4 0.0 0.6 0.0; ...
    0.6 0.0 0.6 0.2; ...
    0.6 0.2 0.4 0.2; ...
    0.4 0.2 0.4 0.0; ...
    ]';
letterGrid([2 4],1:18) = letterGrid([2 4],1:18)*1.6;
letterGrid(:,1:18) = letterGrid(:,1:18) * .75;
letterAspectRatio = 1.6*.75;

% letterGrid = [letterGrid [ ...
%     0 1 1 1; ...
%     1 1 1 0; ...
%     1 0 0 0; ...
%     0 0 0 1; ...
%     ]]';

letterFont{double('A')} = [1 2 6 7 8 9 13];
letterFont{double('B')} = [1 4 6 8 11 13 14];
letterFont{double('C')} = [1 2 9 14];
letterFont{double('D')} = [1 4 6 11 13 14];
letterFont{double('E')} = [1 2 7 8 9 14];
letterFont{double('F')} = [1 2 7 8 9];
letterFont{double('G')} = [1 2 8 9 13 14];
letterFont{double('H')} = [2 6 7 8 9 13];
letterFont{double('I')} = [1 4 11 14];
letterFont{double('J')} = [6 9 13 14];
letterFont{double('K')} = [2 5 7 9 12];
letterFont{double('L')} = [2 9 14];
letterFont{double('M')} = [2 3 5 6 9 13];
letterFont{double('N')} = [2 3 6 9 12 13];
letterFont{double('O')} = [1 2 6 9 13 14];
letterFont{double('P')} = [1 2 6 7 8 9];
letterFont{double('Q')} = [1 2 6 9 12 13 14];
letterFont{double('R')} = [1 2 6 7 8 9 12];
letterFont{double('S')} = [1 2 7 8 13 14];
letterFont{double('T')} = [1 4 11];
letterFont{double('U')} = [2 6 9 13 14];
letterFont{double('V')} = [2 5 9 10];
letterFont{double('W')} = [2 6 9 10 12 13];
letterFont{double('X')} = [3 5 10 12];
letterFont{double('Y')} = [3 5 11];
letterFont{double('Z')} = [1 5 10 14];
letterFont{double('0')} = [1 2 5 6 9 10 13 14];
letterFont{double('1')} = [5 6 13];
letterFont{double('2')} = [1 6 7 8 9 14];
letterFont{double('3')} = [1 6 8 13 14];
letterFont{double('4')} = [2 6 7 8 13];
letterFont{double('5')} = [1 2 7 8 13 14];
letterFont{double('6')} = [1 2 7 8 9 13 14];
letterFont{double('7')} = [1 5 11];
letterFont{double('8')} = [1 2 6 7 8 9 13 14];
letterFont{double('9')} = [1 2 6 7 8 13];
letterFont{double(' ')} = [];
letterFont{double('-')} = [7 8];
letterFont{double('+')} = [4 7 8 11];
letterFont{double('*')} = [3 4 5 10 11 12];
letterFont{double('/')} = [5 10];
letterFont{double('=')} = [7 8 14];
letterFont{double('%')} = [2 5 10 13];
letterFont{double('.')} = [15 16 17 18];
letterFont{double(',')} = [10];
letterFont{double('(')} = [5 12];
letterFont{double(')')} = [3 10];
letterFont{double('[')} = [1 2 9 14];
letterFont{double(']')} = [1 6 13 14];