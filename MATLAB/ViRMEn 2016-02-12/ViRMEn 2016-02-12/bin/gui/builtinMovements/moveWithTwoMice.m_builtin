function velocity = moveWithTwoMice(vr)

velocity = [0 0 0 0];

% Read data from NIDAQ
data = peekdata(vr.ai,50);

% Remove NaN's from the data (these occur after NIDAQ has stopped)
f = isnan(mean(data,2));
data(f,:) = [];
data = mean(data,1)';
data(isnan(data)) = 0;

% Update velocity
offsetTheta = 135;
data = [cosd(offsetTheta) -sind(offsetTheta); sind(offsetTheta) cosd(offsetTheta)]*data([2 1]);
if ~isfield(vr,'scaling')
    vr.scaling = [13 13];
end

velocity(1) = -vr.scaling(1)*data(1);
velocity(2) = -vr.scaling(2)*data(2);

% velocity(3) = atan2(data(2),data(1));

velocity(1:2) = [cos(vr.position(4)) -sin(vr.position(4)); sin(vr.position(4)) cos(vr.position(4))]*velocity(1:2)';