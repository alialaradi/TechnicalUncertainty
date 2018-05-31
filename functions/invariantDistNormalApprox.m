function [ p ] = invariantDistNormalApprox( V, mu, sigma )
%==========================================================================
% Match the invariant distribution of the true reserve amount without 
% learning to a discrete approximation of a normal random variable with 
% mean and variance equal to the initial reserve estiamte and its variance 
% - see Equation (29).
%
% Inputs:
% =======
%   V     = vector of Markov chain states (m x 1)
%   mu    = initial reserve level estimate
%   sigma = variance of initial estimate
%
% Outputs:
% ========
%   p = vector of invariant distribution probabilities   
%==========================================================================

    nStates = length(V);
    
    % initialize invariant distribution
    p = zeros(nStates, 1);
    
    % solve for state 1 probability - Equation (29a)
    p(1) = normcdf((V(1)+V(2))*0.5, mu, sigma) - normcdf((3*V(1)-V(2))*0.5, mu, sigma);
    
    % solve for intermediate states' probabilities - Equation (29b)
    for k = 2: nStates - 1
        p(k) =  normcdf((V(k+1)+V(k))*0.5, mu, sigma) ...
                - normcdf((V(k)+V(k-1))*0.5, mu, sigma);
    end

    % solve for state m probability - Equation (29c)
    p(end) = normcdf((3*V(end)-V(end-1))*0.5, mu, sigma)- normcdf((V(end)+V(end-1))*0.5, mu, sigma);
    
    p = p/sum(p);

end

