function [ A ] = buildGeneratorMatrix( gamma )
%==========================================================================
% Compute generator matrix given transition rates according to Equation (28) 
%
% Inputs:
% =======
%   gamma = vector of transition rates of each state
%
% Outputs:
% ========
%   A = generator matrix for Markov chain
%==========================================================================
           
    % compute A according to Equation (28)
    % gamma(1) are the edge transitions
    % gamme(2:m-1) are interior transitions and reflected around the mid
    % gamma(m) is mid transition
    
    L = length(gamma)-1;
    m = 2*L+1;
    
    A = zeros(m,m);

    A(1,2) = gamma(1);
    A(m,m-1) = gamma(1);
    
    
    for k = 2 : L
        
        A(k,k+1) = gamma(k);
        A(k,k-1) = gamma(k);
        
        A(m-k+1,m-k+2) = gamma(k);
        A(m-k+1,m-k) = gamma(k);
        
    end
    
    A(L+1,L+2) = gamma(L+1);
    A(L+1,L) = gamma(L+1);

    
    A = A - diag( sum(A,2) );

end

