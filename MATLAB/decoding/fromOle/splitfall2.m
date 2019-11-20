load pos
load occu
load ttheta

fall = cellfiring(1.0e4,1.08e4);
save fall1a fall;

fall = cellfiring(1.08e4,1.10e4);
save fall2a fall;


for div=1:4:5   

    load fall1a 
    fall = splitfield(fall,ttheta,div); 
    fname = strcat('fall2a1d',num2str(div))  
    save(fname,'fall'); 

    load fall2a 
    fall = splitfield(fall,ttheta,div); 
    fname = strcat('fall2a1d',num2str(div))  
    save(fname,'fall'); 
end
