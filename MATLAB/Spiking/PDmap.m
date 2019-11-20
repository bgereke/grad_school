function [PDxt,Mx,Mt] = PDmap(apath,atshift,invh,kappa,xbins,tbins)

PDxt = zeros(length(xbins),length(tbins));
TMxt = zeros(length(xbins),length(tbins));
Mt = zeros(length(tbins),1);
Mx = zeros(length(xbins),1);

for t = 1:length(tbins)
    for x = 1:length(xbins)     
        [PDxt(x,t),TMxt(x,t)] = pd_estimator(apath,atshift,xbins(x),tbins(t),invh,kappa,tbins); 
    end
    Mt(t) = (TMxt(:,t)/sum(TMxt(:,t)).*PDxt(:,t)/(PDxt(:,t)'*TMxt(:,t)/sum(TMxt(:,t))))'*log2(PDxt(:,t)/(PDxt(:,t)'*TMxt(:,t)/sum(TMxt(:,t)))); %spatial information (bit/spike)
end

% PDxt = PDxt/sum(sum(PDxt));
% Mt = max(PDxt,[],1);
Mx = max(PDxt,[],2);

% Calculate the rate for one position value
function [pd,tm] = pd_estimator(apath,atshift,x,t,invh,kappa,tbins)
% edge-corrected kernel density estimator
deltax = angle(exp(1j*apath)*conj(exp(1j*x)));
deltat = atshift-t;
deltatb = tbins-t;
% conv_sum = sum(gaussian_kernel(delta,invh));
pd = sum(vmg_kernel(deltax,deltat,invh,kappa));
tdiff = diff(atshift);tdiff(end+1) = median(tdiff);
tdiff(tdiff>1) = median(tdiff);
tm = tdiff'*vm_kernel(deltax,kappa);
pd = pd/tm;
% delta = abs(atpath-i);
% tdiff(delta>cut) = [];
% delta(delta>cut) = [];
% tmap =  tdiff*gaussian_kernel(delta,invh)';
% r = conv_sum/tmap;

% Gaussian kernel for the rate calculation
function pd = vmg_kernel(x,t,invh,kappa)
pd =  exp(-0.5*t.*t*invh^2).*exp(kappa*cos(x))/exp(kappa);

function tm = vm_kernel(x,kappa)
tm = exp(kappa*cos(x))/exp(kappa);

function nm = gaussian_kernel(t,invh)
nm = exp(-0.5*t.*t*invh^2);
