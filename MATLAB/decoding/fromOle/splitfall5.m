load pos
load occu
load ttheta

for div=1:4:5   

    load falla1 
    fall = splitfield(fall,ttheta,div); 
    fname = strcat('falla1d',num2str(div))  
    save(fname,'fall'); 

    load falla2 
    fall = splitfield(fall,ttheta,div); 
    fname = strcat('falla2d',num2str(div))  
    save(fname,'fall'); 

    load falla3 
    fall = splitfield(fall,ttheta,div); 
    fname = strcat('falla3d',num2str(div))  
    save(fname,'fall'); 

    load falla4 
    fall = splitfield(fall,ttheta,div); 
    fname = strcat('falla4d',num2str(div))  
    save(fname,'fall'); 

    load falla5 
    fall = splitfield(fall,ttheta,div); 
    fname = strcat('falla5d',num2str(div))  
    save(fname,'fall'); 


end
