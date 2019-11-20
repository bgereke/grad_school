classdef Dlink < handle
    %bgam.Dlink - An abstract class for distribution/link combos
    %
    %See also bgam.dlink.BinomialLogistic, bgam.dlink.NormalIdentity, bgam.dlink.PoissonExponential
    %
    %The following is of interest if you want to implement a new
    %distribution/link combo:
    %These functions must be implemented:
    %
    % ll = computeLikelihoodAfterNl(this,y,x) 
    % Compute the likelihood an observation x given an observation y
    %
    % ll = computeLikelihood(this,y,eta) 
    % Compute the likelihood of eta
    %
    % [g,H] = computeDlikelihoodX(this,y,eta,X) 
    % Compute the derivative of the likelihood DL/Dw for a model with 
    % internal variable eta+X*w, evaled at w = 0
    %
    % [g] = computeDlikelihoodEta(this,y,eta)
    % Compute the derivative of the likelihood DL/Deta evaled at eta

    properties
        isContinuous = true; %Specifies whether a link/distro combo plays nice with Gauss-Newton
    end
    
    methods(Abstract)
        ll = computeLikelihoodAfterNl(this,y,x) %Compute the likelihood given a reference observation and another
        ll = computeLikelihood(this,y,eta) %Compute the likelihood of eta
        [g,H,L] = computeDlikelihoodX(this,y,eta,X) %Compute the derivative of the likelihood DL/Dw for a model with internal variable eta+X*w, evaled at w = 0
        [g] = computeDlikelihoodEta(this,y,eta) %Compute the derivative of the likelihood DL/Deta
    end
end