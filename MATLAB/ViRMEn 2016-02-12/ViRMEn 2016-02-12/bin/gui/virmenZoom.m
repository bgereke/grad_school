function rect = virmenZoom(guifig,ax)

startPt = get(ax,'currentpoint');

backupCallback = get(guifig,'windowbuttonupfcn');
set(guifig,'windowbuttonupfcn',@endZoom);

zoomEnded = false;
p = plot(ax,NaN,NaN,':b');
tic;
col = hsv(100);
while ~zoomEnded
    endPt = get(ax,'currentpoint');
    set(p,'xdata',[startPt(1,1) startPt(1,1) endPt(1,1) endPt(1,1) startPt(1,1)], ...
        'ydata',[startPt(1,2) endPt(1,2) endPt(1,2) startPt(1,2) startPt(1,2)]);
    t = mod(fix(toc*100),100)+1;
    set(p,'color',col(t,:));
    drawnow;
end
delete(p);
set(guifig,'windowbuttonupfcn',backupCallback);

rect = [min(startPt(1,1),endPt(1,1)) min(startPt(1,2),endPt(1,2)) abs(endPt(1,1)-startPt(1,1)) abs(endPt(1,2)-startPt(1,2))];

    function endZoom(varargin)
        endPt = get(ax,'currentpoint');
        zoomEnded = true;
    end
end