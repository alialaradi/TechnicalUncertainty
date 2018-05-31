function  p_ij = computeTransitionProb(a,b,A,i,t1,t2)
%==========================================================================
% Given a markov chain Z_t with parameters (a,b,A) computes the probability 
% P(Z_(t1) = j | Z_(t2) = i) = [expm( \int_{t1}^{t2} a*exp(-bs) ds )]_{i,j}
% where the generator matrix is given by Eq. 4: G_t = a e^{-bt} A
%
% Inputs:
% =======
%   a,b,A = parameters of generator matrix
%   i     = initial state
%   t1,t2 = time points for transition
%
% Outputs:
% ========
%   p_ij = vector of transition probabilities   
%==========================================================================
        
    if b > 0         
        probMatrix = expm(a/b * (exp(-b*t1) - exp(-b*t2)) * A);            
    else
        probMatrix =  computeInvariantDist( A );
    end
    
    p_ij = probMatrix(i,:);

end

