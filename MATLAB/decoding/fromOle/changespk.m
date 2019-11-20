function w = changespk        



for spk=1:10
     iname = makename3(0,1,spk,1)
     load(iname)
     ot =  100*length(find(ev > 10.0))/length(ev);
     et =  mean(abs(ev)); 
     iname = makename3(1,1,spk,1)

     load(iname)
     ob =  100*length(find(ev > 10.0))/length(ev);
     eb =  mean(abs(ev)); 

     w = [w [spk ot ob et eb]'] 
end

