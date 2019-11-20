function s = makename4(div,spk,phase)

s = 'bt_';
s = strcat(s,'d');                     
s = strcat(s,num2str(div));
s = strcat(s,'s');                     
s = strcat(s,num2str(spk));

if phase  == 0
    s = strcat(s,'r');                     
else
    s = strcat(s,'p');                     
end

