div = 5.0;
for dt=0.050:0.050:0.100       

    phase = 0;
    cont  = 0;
    fname = makename(dt,div,phase,cont); 
    fname
    [xv,tv,pv,ev,sv,nall,occ] =recon(dt,div,phase,cont);
    save(fname,'xv','tv','pv','ev','sv','nall','occ'); 

    phase = 1;
    cont  = 0;
    fname = makename(dt,div,phase,cont); 
    fname
    [xv,tv,pv,ev,sv,nall,occ] =recon(dt,div,phase,cont);
    save(fname,'xv','tv','pv','ev','sv','nall','occ'); 



%   phase = 0;
%   cont  = 1;
%   fname = makename(dt,div,phase,cont); 
%   fname
%
%   [xv,tv,pv,ev,sv,nall,occ] =recon(dt,div,phase,cont);
%   save(fname,'xv','tv','pv','ev','sv','nall','occ'); 
%
%   phase = 1;
%   cont  = 1;
%   fname = makename(dt,div,phase,cont); 
%   fname
%
%   [xv,tv,pv,ev,sv,nall,occ] =recon(dt,div,phase,cont);
%   save(fname,'xv','tv','pv','ev','sv','nall','occ'); 
end




