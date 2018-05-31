addpath('functions')

%% Inputs

% reserve estimate parameters
mu = 10;             % initial estimate
sigma0 = 3;          % initial estimate variance  
sigma1_slow = 2.5;   % reduced variance due to slow learning
sigma1_fast = 1;     % reduced variance due to fast learning
t = 2;               % time to reduced variance

% X_t process model parameters (mean-reverting)
kappa = 0.5;
theta = log(100);
volX  = 0.5;
r     = 0.05;   

% cash-flow parameters
alpha    = 1;
beta 	 = 0.05;
gamma    = 0.9;
epsilon  = 2;
notional = 10^8;
cost     = 20;

% investment cost
fixed_cost    = 10^8;
variable_cost = 3*10^6;

% real option maturity
T = 5; 

% Markov chain parameters
nStates = 31;   % must be odd

%% solve for the Markov chain model
% calibration
nsigma = 2;
V = mu - nsigma * sigma0 + [0 : nStates-1] * 2*nsigma*sigma0 / (nStates-1);

midstate = ceil(nStates/2);

% find the transition rates that replicate the t = 0 distribution
lambda = computeTransitionRates( V, mu, sigma0 );
A      = buildGeneratorMatrix(lambda);


%% price contract using slow and fast learning
% slow learning
params = calibrateLearningParams(V, A, sigma0, sigma1_slow, t);
a_slow = params(1);
b_slow = params(2);

[spot, ~, bd_slow, ~, ~, ~, ~, ~, priceThruTime_slow] = ... 
    realOptionValue(a_slow, b_slow, A, V, T, kappa, theta, volX, r, ...
                    alpha, beta, gamma, epsilon, cost, notional, fixed_cost, variable_cost);

% fast learning
params = calibrateLearningParams(V, A, sigma0, sigma1_fast, t);
a_fast = params(1);
b_fast = params(2);

[~, ~, bd_fast, ~, ~, ~, ~, ~, priceThruTime_fast] = ...
    realOptionValue(a_fast, b_fast, A, V, T, kappa, theta, volX, r, ...
                    alpha, beta, gamma, epsilon, cost, notional, fixed_cost, variable_cost);

% no learning
[~, ~, bd_NL, ~, ~, ~, ~, ~, priceThruTime_NL] = ...
    realOptionValue(1, 0, A, V, T, kappa, theta, volX, r, ...
                    alpha, beta, gamma, epsilon, cost, notional, fixed_cost, variable_cost);

% no running costs
[~, ~, bd_slow_noCost] = realOptionValue(a_slow, b_slow, A, V, T, kappa, theta, volX, r, ...
                            alpha, beta, gamma, epsilon, 0, notional, fixed_cost, variable_cost);
                    
[~, ~, bd_fast_noCost] = realOptionValue(a_fast, b_fast, A, V, T, kappa, theta, volX, r, ...
                            alpha, beta, gamma, epsilon, 0, notional, fixed_cost, variable_cost);

[~, ~, bd_NL_noCost] = realOptionValue(1, 0, A, V, T, kappa, theta, volX, r, ...
                            alpha, beta, gamma, epsilon, 0, notional, fixed_cost, variable_cost);

clc
%% Plots
plotScript
