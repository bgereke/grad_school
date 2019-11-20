% COMPU_INSTA_INSTF: 
% compute instantaneous ampitudes instA and frequencies instF of an IMF
%
% Note: 
% the file ¡°hhspectrum.m¡±is downloaded from http://perso.ens-lyon.fr/patrick.flandrin/software2.html.
% 
% By W.Chenxing, 2015.11.7 at Singapore

function [instA, instF] = Compu_instA_instF(imf)

% calculate instA & instF
[instA,instF,tt] = hhspectrum(imf);

% correct result on borders
instA = [instA(:,1), instA, instA(:,end)];
instF = [instF(:,1), instF, instF(:,end)];

% correct false large values excpet imf1
instFtemp = instF(2:end,:);
instFtemp( find( instFtemp >= 0.4 ) ) = 0.001;
instF(2:end,:) = instFtemp;