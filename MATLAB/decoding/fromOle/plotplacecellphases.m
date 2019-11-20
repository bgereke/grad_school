function plotplacecellphases(w)


tetrode = [1 1 1 1 2 2 3 3 3 3 3 3 3 4 4 5 5 6 6 6 6 6 7 7 7 7 7 8 8 8 8 8 8  11 11 11 11 11]; 
cn      = [1 3 6 7 2 3 1 2 4 5 6 7 8 1 5 1 2 3 4 5 6 7 1 2 6 7 8 1 2 3 4 5 10 1 2  4  6  9]; 
close
for i=1:6
    for j=1:7
        k = j+(i-1)*7; 
        if k<=38
            subplot(7,6,k);
            bar((1:25)/25,w(k,:));
            set(gca,'Xlim',[0 1]);
            s1 = strcat('Tetrode ',num2str(tetrode(k)));
            s2 = strcat(' Cell ',num2str(cn(k)));
            set(gca,'FontSize',8);
            title(strcat(s1,s2)); 
            axis off
        end
    end
end
orient landscape
