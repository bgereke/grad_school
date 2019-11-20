function [GIMFs, j,errone] = RDBR(y,ins_phases,ins_amps,ins_freq,err,J)

%Implements recursive diffeomorphism-based regression (RBDR) 
%Inputs:
%y - 1D signal
%ins_phases - intantaneous phases of each mode [modes x observations]
%ins_amps - instantaneous amplitudes of each mode [modes x observations]
%ins_freq - instantaneous frequencies of each mode [modes x observations]
%err - accuracy parameter
%J - maximum iteration number
%Outputs:
%GIMFs - general intrinsic mode functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = (0:length(y)-1)/length(y);
K = size(ins_phases,1);
GIMFs = zeros(size(ins_phases));
res = y;
err = err*norm(y);
errzero = 2*norm(y);
errone = 1*norm(y);
errtwo = 1*norm(y);
j = 0;

%compute wrapped phases
wrphase = ins_phases;
for k = 1:K
   wrphase(k,:) = mod(ins_phases(k,:),1); 
end

%set options for BSFK
options = struct('periodic',1,'knotremoval_factor',1.01);

%mean correction vars
dt = 0.001;
phgrid = 0:dt:1;

%main loop
while j<J && errone>err && errtwo>err && abs(errone-errzero)>err
    s_tmp = zeros(size(ins_phases));
    tmp_norms = zeros(1,K);
    for k = 1:K
        %line 5
        nonuniform_sample = ins_phases(k,:)/mean(ins_freq(k,:));
        h_k = spline(nonuniform_sample,res./ins_amps(k,:),t);
        %line 7
        pp = BSFK(wrphase(k,:),h_k,2,20,[],options);
%         pp = BSFK(wrphase(k,:),res./ins_amps(k,:),3,30,[],options);
        %line 8 
        mu_tmp = mean(ppval(pp,phgrid));
        s_tmp(k,:) = ppval(pp,wrphase(k,:))-mu_tmp; 
%         s_tmp(k,:) = s_tmp(k,:)-mean(s_tmp(k,:));
        tmp_norms(k) = norm(s_tmp(k,:));
    end
    %line 9
    GIMFs = GIMFs + ins_amps.*(s_tmp);
    %line 10
    res = res - sum(ins_amps.*s_tmp);
    %line 11
    errzero = errone;
    errone = norm(res);
    errtwo = max(tmp_norms);
    %line 12
    j = j+1;
end
