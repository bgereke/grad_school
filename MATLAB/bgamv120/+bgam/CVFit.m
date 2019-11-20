classdef CVFit < bgam.Fit
    %bgam.CVFit - A class that stores a cross-validated Fit
    %
    % bgam.CVFit Properties:
    %    cvdeviance - 
    %    cvdeviances - 
    %    cvd2s - 
    %    subfits - 
    %
    %See Also bgam.Fit

    properties
        %The model CV deviance at the "best" iteration
        cvdeviance;
        
        %The sequence of cross-validated deviance values
        cvdeviances;
        
        %The cross-validated D^2 values
        cvd2s;
        
        %The individual fits performed for each fold
        subfits;
    end
end