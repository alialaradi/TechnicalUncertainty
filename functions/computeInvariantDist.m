function [ p ] = computeInvariantDist( A )
%==========================================================================
% Simulate invariant distribution of Markov chain with generator matrix A
%
% Inputs:
% =======
%   A = Markov chain generator matrix
%
% Outputs:
% ========
%   p = vector of invariant distribution probabilities
%==========================================================================
    
    % transition probabilities for t = 1 and t = 10
    P0 = expm(A);
    P  = expm(10*A);
    
    % initialize invariant distribution 
    p = P(1,:);
    iterate = true;
    
    % compute transition probabilities for t = 11,12,... check if invariant
    % distribution has been reached, iterate until convergence
    while iterate        
        P = P*P0;        
        iterate = max(abs(p - P(1,:))) > 1e-6;       
        p = P(1,:);        
    end

end

