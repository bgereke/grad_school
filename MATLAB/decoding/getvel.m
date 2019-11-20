function [data] = getvel(sessions,data)

%this version finds the error and average TFR for each time window moving
%successively across the session
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameters
scale = .26; % cm per pixel
binsize = 3; % in cm
x = [];
y = [];
xL = [];
v = [];
a = [];
vt = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load position, compute velocity, acceleration for all sessions
for s = 1:length(sessions)
    vfile = [sessions{s},'VT1.nvt'];
    [x_temp,y_temp,vt_temp] = readVideoData(vfile,scale);
    %[x_temp,y] = axesRotater(x_temp,y);  %modify for 2D
    %x_temp = x_temp - min(x_temp);
    x_temp = x_temp - mean(x_temp);
    y_temp = y_temp - mean(y_temp);
    %rotate data to take linear projection
    [p,~] = polyfit(x_temp,y_temp,1); p = p(1);
    xl = [cos(atan(p)) sin(atan(p))]*[x_temp;y_temp];
    yl = [-sin(atan(p)) cos(atan(p))]*[x_temp;y_temp];
    x_temp = xl;y_temp = yl;
    %estimate and smooth velocity (set max to 80 cm/s)
    %v_temp = findVelLinear(x_temp,vt_temp);   
    v_temp = zeros(1,length(x_temp));    
    for i = 2:1:(length(x_temp)-1)
        v_temp(i) = sqrt((x_temp(i+1)-x_temp(i-1))^2 + (y_temp(i+1)-y_temp(i-1))^2)/(vt_temp(i+1)-vt_temp(i-1));
    end
    v_temp(end) = sqrt((x_temp(end)-x_temp(end-1))^2 + (y_temp(end)-y_temp(end-1))^2)/(vt_temp(length(x_temp))-vt_temp(length(x_temp)-1));
%     for i = 2:1:(length(x_temp)-1)
%         v_temp(i) = (x_temp(i+1)-x_temp(i-1))/(vt_temp(i+1)-vt_temp(i-1));
%     end
%     v_temp(end) = (x_temp(end)-x_temp(end-1))/(vt_temp(length(x_temp))-vt_temp(length(x_temp)-1));
    v_temp(v_temp>=80) = 0.5*(v_temp(circshift((v_temp>=80),-3))+v_temp(circshift((v_temp>=80),3)));
    v_temp = smooth(v_temp,15);

    %estimate acceleration
    a_temp = zeros(size(v_temp));
    for i = 2:1:(length(v_temp)-1)
        a_temp(i) = (abs(v_temp(i+1))-abs(v_temp(i-1)))/(vt_temp(i+1)-vt_temp(i-1)); 
    end
    a_temp(length(v_temp)) = (abs(v_temp(length(v_temp)))-abs(v_temp(length(v_temp)-1)))...
        /(vt_temp(length(vt_temp))-vt_temp(length(vt_temp)-1));
    a_temp(a_temp>=80) = 0.5*(a_temp(circshift((a_temp>=80),-3))+a_temp(circshift((a_temp>=80),3)));
    a_temp = smooth(a_temp,15);

    %[xL_temp] = linearizedirection(x_temp,vt_temp,x_temp,vt_temp,v_temp); %consider reflipping leftward runs
    xL_temp = x_temp;
    if size(x_temp,1)==1
        x = [x x_temp];xL = [xL xL_temp];
    else
        x = [x;x_temp];xL = [xL;xL_temp];
    end
    if size(v_temp,1) == 1
        v = [v v_temp];a = [a a_temp];
    else
        v = [v; v_temp];a = [a; a_temp];
    end
    if size(vt_temp,1) == 1
        vt = [vt vt_temp];
    else
        vt = [vt; vt_temp];
    end
        
end
'stillworking3'
%get "measured position" in terms of bin for the time window
data.x = nan(size(data.w,1),1); %mean position for window
data.vel = nan(size(data.w,1),1); %mean velocity for window
data.acc = nan(size(data.w,1),1); %mean accerlation for window
data.wt = nan(size(data.w,1),1); %mean time for window
tLength = max(xL) - min(xL);
numBins = ceil(tLength / binsize);
data.xbins(:,1) = binsize*(1:numBins)+min(xL)-0.5*binsize; %location bins

for ww = 1:size(data.w,1)  
    xspan = xL(vt>data.w(ww,1) & vt<=data.w(ww,2));
    vspan = v(vt>data.w(ww,1) & vt<=data.w(ww,2));
    aspan = a(vt>data.w(ww,1) & vt<=data.w(ww,2));
    tspan = vt(vt>data.w(ww,1) & vt<=data.w(ww,2));
    x1 = mean(xspan);
    xdiff = (data.xbins - x1).^2;
    [~,bin_idx] = min(xdiff);
    data.x(ww) = data.xbins(bin_idx);
    data.vel(ww) = mean(vspan);
    data.acc(ww) = mean(aspan);
    data.wt(ww) = mean(tspan);
end
'stillworking4'        

    
