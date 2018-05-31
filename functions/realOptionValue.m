function [spot, price, bd, Investment, payoff, ProjectValue_FixedVol, p_ij, vpA, priceThruTime] = ...
    realOptionValue(a,b,A,MC,T,kappa,theta,volX,r,alpha,beta,gamma,epsilon,cost,notional,fixed_cost,variable_cost)
%==========================================================================
% Compute real option value for given parameter set. Note: to compare the
% case of (a,b,A)-learning with no learning, set a = 1, b = 0, and 
% A = a/b * A - this ensures that the starting distributions of \vartheta
% in the two cases match (but the no-learning case does not evolve) - see
% Section 3 and 4 for details.
%
% Inputs:
% =======
%   a,b,A,MC = Markov chain parameters (a,b,A) and states (MC)
%   T        = real option maturity
%   kappa,theta,volX 
%            = spot price dynamics, see Equation (5)
%   r = discount rate
%   alpha,beta,epsilon,notional,cost
%            = cash flow parameters, see Equation (7)
%   fixed_cost,variable_cost
%            = investment cost parameters, see Equation (12)
%
% Outputs:
% ========
%   spot  = vector of spot prices
%   price = matrix of real option prices (nSpotPrices x nStates) at t = 0
%   bd    = exercise boundary (nPeriods x nStates)
%   Investment = investment cost in each state 
%   payoff = payoff in each spot/state combination at option maturity
%   ProjectValue_FixedVol 
%          = project value at each state - Equation (11)
%   p_ij   = transition prob to/from states from given time to option mat.
%   vpA    = invariant distribution of true reserve amount based on MC
%   priceThruTime 
%          = option price for each spot/state at each time point
%==========================================================================
    
    %% set up 

    M  = 51*T; % number of exercise dates (weekly steps)
    dt = T/M;  % time increment

    % grid size
    N = 2^13;

    % Real space - set up grid from -10 to 10
    xm      = 10;
    x_1_min = -xm; 
    x_1_max = xm;
    dx_1    = (x_1_max-x_1_min)/(N-1);
    x_1     = x_1_min:dx_1:x_1_max;
    
    spot = exp(theta + x_1);

    % Fourier space
    w_1_max = pi/dx_1;
    dw_1    = 2*w_1_max/N;
    w_1     = [0:dw_1:w_1_max, -w_1_max+dw_1:dw_1:-dw_1];
    
    % vector of ones
    onesN = ones(1,N);     
    
    vol = MC;   % reserve estimate (Markov chain states)
    Nstates = size(vol,2);
    
    % invariant distritubion of exp(A)    
    vpA = computeInvariantDist( A );
  
    %% generate intrinsic option value   
    
    % Compute project value and investment cost in each state - Equations
    % (11) and (12)
    ProjectValue_FixedVol = zeros(N, Nstates);
    Investment = zeros(N,Nstates);    
    
    for i = 1 : Nstates
        
        % compute project value for each state - integral in Equation (11)
        % - uses helper function "ProjectVal" defined below
        Delta = -1/beta * log( 1 - beta * gamma * vol(i) / alpha );
        ProjectValue_FixedVol(:,i) = quadvgk( @(u)(ProjectVal(u,x_1,onesN,alpha)), ...
                                                            [epsilon; epsilon + Delta], N);

        % Investment costs
        Investment(:,i) = fixed_cost + vol(i) * variable_cost;
        
    end
    
    % compute intrinsic value of option   
    if b ~= 0  % in the case with learning  
        
        payoff = max( ProjectValue_FixedVol * expm( (a/b*exp(-b*T)) * A ) - Investment, 0);
    
    else % no-learning case
        
        % use mid state investment cost (since we don't transition)
        Investment = repmat( Investment(:,ceil(Nstates/2)), 1, Nstates);
        % in no-learning case, use invariant distribution to compute exp.        
        payoff = max( repmat( ProjectValue_FixedVol * vpA',1, Nstates) - Investment, 0);
    end
        
    price = payoff;

    % helper function: project value integrand if exercised at level x - integrand in Equation (11)
    function f = ProjectVal(u,x,v, alpha)
        F = exp( theta + x' * exp(-kappa*u)  + volX^2/(4*kappa) * v' * ( 1 - exp( - 2*kappa*u )) );
        f = notional * alpha * exp( v' * (- r * u - beta * (u - epsilon)) )  .* (F - cost);
    end
     
    % characteristic function  
    expa = (exp(2*kappa*dt) - 1)/(2*kappa);

    char_exp_factor = zeros(N,Nstates);
    for j = 1 : Nstates
        char_exp_factor(:,j) = exp( - 1/2 * volX^2 * expa * w_1.^2 );
    end

    % interpolation points
    x_1interp = exp(-kappa*dt) * x_1;

    %% search for trigger function
    bd = NaN * ones(M+1,Nstates);

    exer = (payoff > 0);

    for j = 1: Nstates 
        
        idx = find(exer(:,j) == 1, 1, 'first');
        if ~isempty(idx)
           bd(1,j) = x_1(idx); 
        end
        
    end

    %% mrFST - backward stepping; step 4 of algorithm described in Figure 2
    fftw('planner', 'measure');

    tic
    for l = 1 : M

        % interpolate
        interp_price = interp1(x_1, price, x_1interp);

        % transition rate per step
        if b == 0
            % in no-learning case, there are no transitions
            Lambda = 0;
        else
            Lambda = a * ( exp(-b* (T-(l)*dt) )  - exp(-b* (T-(l-1)*dt) ) )/b;            
        end

        % mrFST
        hold_price = exp(-r*dt) * real(ifft( (fft(interp_price) * expm(Lambda * A ) ) .* char_exp_factor));

        % exercise
        if b ~= 0
            intrinsic = ProjectValue_FixedVol * expm( (a/b*exp(-b*(T-dt*l))) * A ) - Investment;
            p_ij(:,:,l) = expm( (a/b*exp(-b*(T-dt*l))) * A );
        else            
            intrinsic = repmat( ProjectValue_FixedVol * vpA',1, Nstates) - Investment;
            p_ij = 0;
        end

        exer = ( hold_price <= intrinsic - 1e-5 );

        price = ( intrinsic - hold_price) .* exer + hold_price ;

        % search for boundary function
        for j = 1: Nstates 
            
            idx = find(exer(:,j) == 1, 1, 'first');
            if ~isempty(idx)
               bd(1+l,j) = x_1(idx); 
            end
          
        end
        
        priceThruTime(:,:,l) = price;

    end   
end