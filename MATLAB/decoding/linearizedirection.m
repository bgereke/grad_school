function xl = linearizedirection(xnonL,tnonL,x,t,v)

xspan = max(x) - min(x);
xl = xnonL;
% for i = 1:length(xnonL)
%     
%    tdiff = (t - tnonL(i)).^2;
%    [~,ind] = min(tdiff);
%   
%    %for separating out leftward runs
% %    if v(ind) < 0
% %     xl(i) = (xnonL(i) *-1 ) + xspan*2;
% %    end
%        
% end