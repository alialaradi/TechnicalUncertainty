function MM = computeTrueReserveMoments(a, b, A, V, t)
%==========================================================================
% Compute first four moments of \vartheta based on initial estimate
%
% Inputs:
% =======
%   a,b,A = parameters of generator matrix
%   V     = vector of Markov chain states
%   t     = current time
%
% Outputs:
% ========
%   MM = vector of first four moments of \vartheta
%==========================================================================
    
    % compute midstate (initial estimate)
    midstate = ceil(length(V)/2);
    
    % compute probability of converging to each state in Markov chain
    p_ij = computeTransitionProb(a,b,A,midstate,t,Inf);
    
    % compute moments using distribution given above
    MM(1,1) = sum(p_ij .* V );
    MM(2,1) = sum(p_ij .* (MM(1)-V).^2 );
    MM(3,1) = sum(p_ij .* (MM(1)-V).^3 );
    MM(4,1) = sum(p_ij .* (MM(1)-V).^4 );       
    
end