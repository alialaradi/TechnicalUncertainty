function [ lambda ] = computeTransitionRates( V, mu, sigma )
%==========================================================================
% Compute transition rates for Markov chain generator matrix by matching
% the implied (simulated) invariant distribution to a target invariant 
% distribution approximated using a normal distribution with mean/variance 
% equal to the initial reserve estimate/estimate variance - see Section 4.
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
    L = ceil(nStates-1)/2;
    
    % compute target invariant distribution
    p = invariantDistNormalApprox( V, mu, sigma );        
    
    % define helper function to use in solver
    function f = myfunc(ll)        
        
        % compute lambda vector for given input
        gamma = 1./(1+exp(ll));
        
        % compute corresponding generator matrix
        A = buildGeneratorMatrix(gamma);
        
        % simulate invariant distribution
        inv_dist = computeInvariantDist( A );
        
        % compute distance between simulated and target invariant dist.
        f = inv_dist' - p;
        
    end

    % solve for optimal lambda to match target invariant distribution
    myoptimset = optimset('MaxFunEvals', 5000*L, 'TolFun',1e-6);        
    l = fsolve(@myfunc, cumsum(0.05*ones(L+1,1)), myoptimset);        
    lambda = 1./(1+exp(l));
    
end


