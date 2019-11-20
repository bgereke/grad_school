[peakRate, Vb6, pathCoord] = equalPlot_max_original('infile.txt',0);
mean_Vb6 = zeros(3,1);
std_Vb6 = zeros(3,1);
stder_Vb6 = zeros(3,1);
med_Vb6 = zeros(3,1);
sub_Vb6 = cell(3,1);
sub = 150;

for i=1:3
    mean_Vb6(i) = mean(Vb6{i,1});
    std_Vb6(i) = std(Vb6{i,1});
    stder_Vb6(i) = std(Vb6{i,1})/sqrt(length(Vb6{i,1}));
    med_Vb6(i)= median(Vb6{i,1});
    
    for ii=1:floor(length(Vb6{i,1})/sub)
        sub_Vb6{i,1}(ii) = mean(Vb6{i,1}((sub*(ii-1)+1):sub*ii));
    end

end





boxplot(data,{'B6','3XTG'},'notch','on')


[peakRate, Vad, pathCoord] = equalPlot_max_original('infile.txt',0);
mean_Vad = zeros(3,1);
std_Vad = zeros(3,1);
stder_Vad = zeros(3,1);
med_Vad = zeros(3,1);
sub_Vad = cell(3,1);
sub = 150;

for i=1:3
    mean_Vad(i) = mean(Vad{i,1});
    std_Vad(i) = std(Vad{i,1});
    stder_Vad(i) = std(Vad{i,1})/sqrt(length(Vad{i,1}));
    med_Vad(i)= median(Vad{i,1});

    for ii=1:floor(length(Vad{i,1})/sub)
        sub_Vad{i,1}(ii) = mean(Vad{i,1}((sub*(ii-1)+1):sub*ii));
    end

end


data=[Vb6{1,1},Vad{1,1}]; 
boxplot(data,{'B6','3XTG'},'notch','on')

plot(plot:Boxplot(data1, data2, Notched == TRUE)):


