function [w, trans] = virmenLoadWindows(exper)

w = zeros(5,length(exper.windows));
trans = zeros(1,length(exper.windows));
mnt = virmenOpenGLMonitors();

for ndx = length(exper.windows):-1:1
    if exper.windows{ndx}.primaryMonitor
        indx = 1;
    else
        indx = exper.windows{ndx}.monitor+1;
    end
    if indx > size(mnt,2)
        indx = 1;
    end
    if exper.windows{ndx}.fullScreen
        w(1:4,ndx) = mnt(:,indx);
    else
        w(1:4,ndx) = [exper.windows{ndx}.left+mnt(1,indx) mnt(4,indx)-(exper.windows{ndx}.bottom+mnt(2,indx))-exper.windows{ndx}.height ...
            exper.windows{ndx}.width exper.windows{ndx}.height]';
    end
    w(5,ndx) = exper.windows{ndx}.antialiasing;
    
    if exper.windows{ndx}.rendering3D
        trans(ndx) = exper.windows{ndx}.transformation;
    else
        trans(ndx) = NaN;
    end
end