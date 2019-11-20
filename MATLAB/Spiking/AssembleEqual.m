function AssembleEqual(TTList,imagedirectory,sequence,x,y)


if nargin == 2 sequence = [1 2 3 4 5 6 7 8]; end
Directory = pwd

Session1 = strcat('begin',num2str(sequence(1)),'\',imagedirectory);
Session2 = strcat('begin',num2str(sequence(2)),'\',imagedirectory);
Session3 = strcat('begin',num2str(sequence(3)),'\',imagedirectory);
Session4 = strcat('begin',num2str(sequence(4)),'\',imagedirectory);
Session5 = strcat('begin',num2str(sequence(5)),'\',imagedirectory);
Session6 = strcat('begin',num2str(sequence(6)),'\',imagedirectory);
Session7 = strcat('begin',num2str(sequence(7)),'\',imagedirectory);
Session8 = strcat('begin',num2str(sequence(8)),'\',imagedirectory);

TPath = strcat(Directory,'\',TTList)
cluster = ReadFileList(TPath);

cellspersheet = 7;
%empty = ones(402,402,3);
empty = ones(420,560,3);

for (i = 1:length(cluster))

   bmpfile = strrep(char(cluster(i)),'.t','.bmp')
%    bmpfile = strcat(char(cluster{i}(1:end)),'.bmp')

   ImgPathA = strcat(Directory,'\',Session1,'\',bmpfile);
   ImgPathB = strcat(Directory,'\',Session2,'\',bmpfile);
   ImgPathC = strcat(Directory,'\',Session3,'\',bmpfile);
   ImgPathD = strcat(Directory,'\',Session4,'\',bmpfile);
   ImgPathE = strcat(Directory,'\',Session5,'\',bmpfile);
   ImgPathF = strcat(Directory,'\',Session6,'\',bmpfile);
   ImgPathG = strcat(Directory,'\',Session7,'\',bmpfile);
   ImgPathH = strcat(Directory,'\',Session8,'\',bmpfile);

   if (fopen(ImgPathA) >= 0) A = imread(ImgPathA);empty = ones(size(A)); else A = empty.*255; end;
   if (fopen(ImgPathB) >= 0) B = imread(ImgPathB);empty = ones(size(A)); else B = empty.*255; end;
   if (fopen(ImgPathC) >= 0) C = imread(ImgPathC);empty = ones(size(A)); else C = empty.*255; end;
   if (fopen(ImgPathD) >= 0) D = imread(ImgPathD);empty = ones(size(A)); else D = empty.*255; end;
   if (fopen(ImgPathE) >= 0) E = imread(ImgPathE);empty = ones(size(A)); else E = empty.*255; end;
   if (fopen(ImgPathF) >= 0) F = imread(ImgPathF);empty = ones(size(A)); else F = empty.*255; end;
   if (fopen(ImgPathG) >= 0) G = imread(ImgPathG);empty = ones(size(A)); else G = empty.*255; end;
   if (fopen(ImgPathH) >= 0) H = imread(ImgPathH);empty = ones(size(A)); else H = empty.*255; end;
   % keyboard
   I = [A B C D E F G H];

   if (mod(i,cellspersheet) == 1) J = I; else J = [J; I]; end;
   if nargin == 2
       if (mod(i,cellspersheet) == 0) imwrite(J,strcat(Directory,'\',strrep(TTList,'.',''),imagedirectory,num2str(i/cellspersheet),'.jpg')); end;
   else
       if (mod(i,cellspersheet) == 0) imwrite(J,strcat(Directory,'\',strrep(TTList,'.',''),imagedirectory,'_ordered',num2str(i/cellspersheet),'.jpg')); end;
   end
end
if nargin == 2
   imwrite(J,strcat(Directory,'\',strrep(TTList,'.',''),imagedirectory,'last','.jpg'));
else
   imwrite(J,strcat(Directory,'\',strrep(TTList,'.',''),imagedirectory,'_ordered','last','.jpg'))
end
    

%image(J)
%axis off
%axis image
