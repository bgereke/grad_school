function s = makename(dt,div,phase,bay)

if bay == 1
    s = 'bayp_';
else
    s = 'tmpp_';
end

s = strcat(s,'d');
s = strcat(s,num2str(div));

s = strcat(s,'dt');
s = strcat(s,num2str(1000*dt));

if phase == 0
    s = strcat(s,'r');
else
    s = strcat(s,'p');
end

