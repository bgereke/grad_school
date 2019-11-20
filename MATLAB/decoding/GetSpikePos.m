function [spkx,spky,spkv,spka] = GetSpikePos(ts,posx,posy,post,posv,posa)

N = length(ts);
spkx = zeros(N,1);
spky = zeros(N,1);
spkv = zeros(N,1);
spka = zeros(N,1);
for ii = 1:N
    tdiff = (post-ts(ii)).^2;
    [~,ind] = min(tdiff);
%     if v(ind(1))<0 %if want to dpownsample by direction
    spkx(ii) = posx(ind(1));
    spky(ii) = posy(ind(1));
    spkv(ii) = posv(ind(1));
    spka(ii) = posa(ind(1));
%     end
end