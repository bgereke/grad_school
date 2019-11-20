function mymap = mycolormap(N,color)

mh=[0,0;0.3,1;1,1];
mm=[0,0;0.3,0;0.7,1;1,1];
ml=[0,0;0.7,0;1,1];
x = linspace(0,1,N);
if strcmp(color,'red')
    rv = interp1( mh(:,1), mh(:,2), x);
    gv = interp1( mm(:,1), mm(:,2), x);
    mv = interp1( ml(:,1), ml(:,2), x);
elseif strcmp(color,'blue')
    rv = interp1( ml(:,1), ml(:,2), x);
    gv = interp1( mm(:,1), mm(:,2), x);
    mv = interp1( mh(:,1), mh(:,2), x);
else
    rv = interp1( mm(:,1), mm(:,2), x);
    gv = interp1( mh(:,1), mh(:,2), x);
    mv = interp1( ml(:,1), ml(:,2), x);
end
mymap = [ rv', gv', mv'];
%exclude invalid values that could appear
mymap( isnan(mymap) ) = 0;
mymap( (mymap>1) ) = 1;
mymap( (mymap<0) ) = 0;