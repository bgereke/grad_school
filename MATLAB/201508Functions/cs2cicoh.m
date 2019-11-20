function cicoh=cs2cicoh(cs)

%function [cicoh]=cs2cicoh(cs)
%
% Calculates the "corrected" imaginary part of coherency for all channel
% pairs from a given complex cross-spectrum
%
% see A. Ewald, L. Marzetti, F. Zappasodi, F. C. Meinecke, and G. Nolte.
% Estimating true brain connectivity from EEG/MEG data invariant to linear
% and static transformations in sensor space. NeuroImage, 60:476 – 488, 2012.
% http://dx.doi.org/10.1016/j.neuroimage.2011.11.084
%
%
% IN   cs        - sensor level complex cross-spectrum at a single frequency bin;
%                  might have been overaged over frequencies
%                 (chans x chans)
%
% OUT  cicoh     - "corrected" imaginary part of coherency
%
% History
% Arne Ewald (mail@aewald.net), 15.03.2015

% License
%   Copyright (C) 2015  Arne Ewald
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see http://www.gnu.org/licenses/.

[nch,~]=size(cs);
coh=cs./sqrt(diag(cs)*diag(cs)');

imcoh=imag(coh);
rcoh=real(coh);

cicoh=zeros(nch,nch);
for p=1:nch
    for q=1:nch
        if p~=q
            cicoh(p,q)=imcoh(p,q)/sqrt(1-(rcoh(p,q))^2);
        end
    end
end
