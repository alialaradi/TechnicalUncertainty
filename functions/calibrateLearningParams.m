function [ params ] = calibrateLearningParams(V, A, sigma1, sigma2, t)
%==========================================================================
% Compute learning parameters in h_t = a e^{-bt} as described in Section 4.2
%
% Inputs:
% =======
%   V      = vector of Markov chain states (m x 1)
%   A      = Markov chain generator matrix
%   sigma1 = variance of initial estimate
%   sigma2 = variance of estimate after learning 
%   t      = time frame for variance reduction
%
% Outputs:
% ========
%   params = vector of learning parameters a,b
%==========================================================================

    fprintf('calibrating a, b for the learning agent...\n');

    var1 = sigma1^2;
    var2 = sigma2^2;
    
    % define helper function to use in solver
    function f = myfunc_with_learning(b)
        
        b = exp(b);
        
        % compute moments based on guess at time 0 and time t 
        M_0 = computeTrueReserveMoments(b(1), b(2), A, V, 0); 
        M_T = computeTrueReserveMoments(b(1), b(2), A, V, t);
        
        % compute distance between guess-implied second moments and
        % variance parameters given by learning process
        f = [M_0(2); M_T(2);] - [var1; var2;];
        
    end

    % solve for learning parameters to match second moments
    bb = fsolve( @myfunc_with_learning ,[0 -1]);
    
    % save parameters for output
    params(1:2) = exp(bb);

end

