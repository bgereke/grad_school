function s = makename3(method,div,spk,phase)

if method == 1
    s = 'bay3_';
else
    s = 'tem3_';
end
s = strcat(s,'d');                     
s = strcat(s,num2str(div));
s = strcat(s,'s');                     
s = strcat(s,num2str(spk));

if phase  == 0
    s = strcat(s,'r');                     
else
    s = strcat(s,'p');                     
end

