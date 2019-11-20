function AssembleFiles(filetype,sequence,x,y)

%filetype must be a string: '*.filetype'

if nargin == 1 sequence = [1 2 3 4 5 6 7 8]; end
Directory = pwd;

% Session1 = strcat('Begin',num2str(sequence(1)),'\',imagedirectory);
% Session2 = strcat('Begin',num2str(sequence(2)),'\',imagedirectory);
% Session3 = strcat('Begin',num2str(sequence(3)),'\',imagedirectory);
% Session4 = strcat('Begin',num2str(sequence(4)),'\',imagedirectory);
% Session5 = strcat('Begin',num2str(sequence(5)),'\',imagedirectory);
% Session6 = strcat('Begin',num2str(sequence(6)),'\',imagedirectory);
% Session7 = strcat('Begin',num2str(sequence(7)),'\',imagedirectory);
% Session8 = strcat('Begin',num2str(sequence(8)),'\',imagedirectory);

% TPath = strcat(Directory,'\',TTList)
% cluster = ReadFileList(TPath);
files = dir(filetype);
filetype = filetype(end-2:end);

numcols = 7;
numrows = 7;
J = cell(numrows, numcols);
col = 1;
row = 1;

for i = 1:size(files,1)

    %    bmpfile = strrep(char(cluster(i)),'.t','.bmp')
    typefile = files(i).name;
    
    ImgPath = strcat(Directory,'\',typefile);
    [I,map] = imread(ImgPath);
    I = ind2rgb(I,map);

    %    ImgPathB = strcat(Directory,'\',Session2,'\',bmpfile);
    %    ImgPathC = strcat(Directory,'\',Session3,'\',bmpfile);
    %    ImgPathD = strcat(Directory,'\',Session4,'\',bmpfile);
    %    ImgPathE = strcat(Directory,'\',Session5,'\',bmpfile);
    %    ImgPathF = strcat(Directory,'\',Session6,'\',bmpfile);
    %    ImgPathG = strcat(Directory,'\',Session7,'\',bmpfile);
    %    ImgPathH = strcat(Directory,'\',Session8,'\',bmpfile);
    
    
    %    if (fopen(ImgPathB) >= 0) B = imread(ImgPathB);empty = ones(size(A)); else B = empty.*255; end;
    %    if (fopen(ImgPathC) >= 0) C = imread(ImgPathC);empty = ones(size(A)); else C = empty.*255; end;
    %    if (fopen(ImgPathD) >= 0) D = imread(ImgPathD);empty = ones(size(A)); else D = empty.*255; end;
    %    if (fopen(ImgPathE) >= 0) E = imread(ImgPathE);empty = ones(size(A)); else E = empty.*255; end;
    %    if (fopen(ImgPathF) >= 0) F = imread(ImgPathF);empty = ones(size(A)); else F = empty.*255; end;
    %    if (fopen(ImgPathG) >= 0) G = imread(ImgPathG);empty = ones(size(A)); else G = empty.*255; end;
    %    if (fopen(ImgPathH) >= 0) H = imread(ImgPathH);empty = ones(size(A)); else H = empty.*255; end;
    
    %    I = [A B C D E F G H];
    
    if row*col == numrows*numcols, 
        J{row,col} = I;
        row = 1;col=1;
        J = cell2mat(J);
        imwrite(J,strcat(Directory,'\',filetype,num2str(i/(numrows*numcols)),'.jpg'));
        J = cell(numrows,numcols);
        continue
    end
    J{row,col} = I;
    if col < numcols
        col = col + 1;
    else
        col = 1;
        row = row + 1;
    end
end
if row*col < numrows*numcols && row*col > 1
    while row*col <= numrows*numcols
       empty = ones(size(J{1,1}));
       J{row, col} = empty*255;
       if col < numcols
        col = col + 1;
       else
           col = 1;
           row = row + 1;
       end
    end
end
if row*col > 1
    if nargin == 1 
        J(end,:) = [];
        J = cell2mat(J); 
        imwrite(J,strcat(Directory,'\',filetype,'last','.jpg'));
    else
        imwrite(J,strcat(Directory,'\',filetype,'_ordered','last','.jpg'))
    end
end
    

%image(J)
%axis off
%axis image
