function figs = figureLayout(figs,layout)

margin = 5;
global guifig;
set(guifig,'units','pixels');
pixelpos = get(guifig,'position');
set(guifig,'units','normalized');
marginX = margin/pixelpos(3);
marginY = margin/pixelpos(4);

allFigs = fieldnames(figs);
showFigs = fieldnames(layout);
noshowFigs = setdiff(allFigs,showFigs);
for ndx = 1:length(showFigs)
    if strcmp(get(figs.(showFigs{ndx}),'visible'),'off')
        set(figs.(showFigs{ndx}),'visible','on');
    end
    pos = layout.(showFigs{ndx});
    pos(1) = pos(1)+marginX/2;
    pos(2) = pos(2)+marginY/2;
    pos(3) = pos(3)-marginX;
    pos(4) = pos(4)-marginY;
    set(figs.(showFigs{ndx}),'position',pos);
end
for ndx = 1:length(noshowFigs)
    if strcmp(get(figs.(noshowFigs{ndx}),'visible'),'on')
        set(figs.(noshowFigs{ndx}),'visible','off');
    end
end

updateFigures('LAYOUT');
rotate3d(guifig,'off');