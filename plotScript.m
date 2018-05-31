close all

% time grid for x-axis
timeGrid = [0 : 1 : 51*T]/51;

% latex text in figures
set(0,'defaultTextInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

% states to plot
statesToPlot = midstate-4 : 2 : midstate+4;

%% conditional distribution of true reserve amount (calibrated vs. approx) - FIGURE 3
f = figure;
set(f,'Position',[100 100 900 400]);

% normal approximation of invariant distributions
normApprox0     = invariantDistNormalApprox(V, mu, sigma0);
normApprox_slow = invariantDistNormalApprox(V, mu, sigma1_slow);
normApprox_fast = invariantDistNormalApprox(V, mu, sigma1_fast);

% transition probabilities at time 0 and t for starting at midstate (slow)
reserveDist_slow_0 = computeTransitionProb(a_slow,b_slow,A,midstate,0,Inf);
reserveDist_slow_t = computeTransitionProb(a_slow,b_slow,A,midstate,t,Inf);

% transition probabilities at time 0 and t for starting at midstate (fast)
reserveDist_fast_0 = computeTransitionProb(a_fast,b_fast,A,midstate,0,Inf);
reserveDist_fast_t = computeTransitionProb(a_fast,b_fast,A,midstate,t,Inf);

% smooth version of normal approximation (for smooth lines)
V_fine  = linspace(V(1),V(end),100);
normApprox0     = pchip(V,normApprox0,V_fine);
normApprox_slow = pchip(V,normApprox_slow,V_fine);
normApprox_fast = pchip(V,normApprox_fast,V_fine);

% slow learning plot
subplot(1,2,1)
hold on
plot(V, reserveDist_slow_0, 'ob', V, reserveDist_slow_t, '*r')
plot(V_fine, normApprox0, '-b', V_fine, normApprox_slow, '-r')
xlabel('Reserve Estimate $v^{(j)}$', 'fontsize', 14);
ylabel('$P \left( \vartheta = v^{(j)} \; | \; V_t = \mu \right)$', 'fontsize', 14);
leg = legend('$t=0$ Model','$t=2$ Model','$t=0$ Normal Approx.','$t=2$ Normal Approx.');
set(leg,'fontsize',10)
grid on
ylim([0 0.18])
title('\textbf{True Reserve Dist. - Slow learning}','fontsize',12)

% fast learning plot
subplot(1,2,2)
hold on
plot(V, reserveDist_slow_0, 'ob', V, reserveDist_fast_t, '*r')
plot(V_fine, normApprox0, '-b', V_fine, normApprox_fast, '-r')
xlabel('Reserve Estimate $v^{(j)}$', 'fontsize', 14);
ylabel('$P \left( \vartheta = v^{(j)} \; | \; V_t = \mu \right)$', 'fontsize', 14);
% leg = legend('$t=0$ Model','$t=2$ Model','$t=0$ Normal Approx.','$t=2$ Normal Approx.');
% set(leg,'fontsize',10,'position',[0.781 0.716 0.188 0.183])
grid on
title('\textbf{True Reserve Dist. - Fast learning}','fontsize',12)

%% plot exercise boundaries - FIGURE 4
% compare 
% statesToPlot = 2:midstate and statesToPlot = midstate:nStates-1
% to compare what's happening in low vs. high states
% in low states:  exercise boundary shifts down with increased estimate with
%                 faster acceleration in the early part
% in high states: exercise boundary ends up lower with increased estimate
%                 but starts out higher with faster deceleration
% again something weird is happening in edge states (1 and nStates)
% with no learning, exercise boundary is not state-dependent (investment 
% cost is not state-dependent and expectation is based on invariant 
% distribution from time t=0)

f = figure;
set(f,'Position',[100 100 900 400]);

rval = linspace(0.3,0.95,length(statesToPlot));
bval = 0.2;
gval = 0.2;

% slow learning
subplot(1,2,1)
hold on;
for i = 1:length(statesToPlot)
    currState = statesToPlot(i);
    plot(timeGrid(1:end-1),exp(bd_slow(end:-1:2,currState) + theta), 'Color', [rval(i),gval,bval], 'LineWidth', 1.5)
    legendName{i} = ['$V_t$ = ' num2str(currState)];
end
plot(timeGrid(1:end-1),exp(bd_NL(end:-1:2,midstate) + theta), 'b--', 'LineWidth', 1.5)
plot(timeGrid(1:end-1),exp(bd_slow(end:-1:2,midstate) + theta), 'w--', 'LineWidth', 2)
grid on
xlabel('Time','FontSize',12)
ylabel('Exercise Price','FontSize',12)
title('\textbf{Slow Learning - with running costs}','FontSize',12)
legend(legendName,'Location','SouthEast','fontsize',10)
ylim([20 60])

% fast learning
subplot(1,2,2)
hold on;
for i = 1:length(statesToPlot)
    currState = statesToPlot(i);
    plot(timeGrid(1:end-1),exp(bd_fast(end:-1:2,currState) + theta), 'Color', [rval(i),gval,bval], 'LineWidth', 1.5)
end
plot(timeGrid(1:end-1),exp(bd_fast(end:-1:2,midstate) + theta), 'w--', 'LineWidth', 2)
plot(timeGrid(1:end-1),exp(bd_NL(end:-1:2,midstate) + theta), 'b--', 'LineWidth', 1.5)
grid on
xlabel('Time','FontSize',12)
ylabel('Exercise Price','FontSize',12)
title('\textbf{Fast Learning - with running costs}','FontSize',12)
ylim([20 60])

%% plot exercise boundaries - no running costs - FIGURE 5
% compare 
% statesToPlot = 2:midstate and statesToPlot = midstate:nStates-1
% to compare what's happening in low vs. high states
% in low states:  exercise boundary shifts down with increased estimate with
%                 faster acceleration in the early part
% in high states: exercise boundary ends up lower with increased estimate
%                 but starts out higher with faster deceleration
% again something weird is happening in edge states (1 and nStates)
% with no learning, exercise boundary is not state-dependent (investment 
% cost is not state-dependent and expectation is based on invariant 
% distribution from time t=0)

%%% states to plot
statesToPlot = midstate-4 : 2 : midstate+4;
%%%

f = figure;
set(f,'Position',[100 100 900 400]);

rval = linspace(0.3,0.95,length(statesToPlot));
bval = 0.2;
gval = 0.2;

% slow learning
subplot(1,2,1)
hold on;
for i = 1:length(statesToPlot)
    currState = statesToPlot(i);
    plot(timeGrid(1:end-1),exp(bd_slow_noCost(end:-1:2,currState) + theta), 'Color', [rval(i),gval,bval], 'LineWidth', 1.5)
    legendName{i} = ['$V_t$ = ' num2str(currState)];
end
plot(timeGrid(1:end-1),exp(bd_NL_noCost(end:-1:2,midstate) + theta), 'b--', 'LineWidth', 1.5)
plot(timeGrid(1:end-1),exp(bd_slow_noCost(end:-1:2,midstate) + theta), 'w--', 'LineWidth', 2)
grid on
xlabel('Time','FontSize',12)
ylabel('Exercise Price','FontSize',12)
title('\textbf{Slow Learning - no running costs}','FontSize',12)
legend(legendName,'Location','Northwest','fontsize',10)
ylim([20 60])

% fast learning
subplot(1,2,2)
hold on;
for i = 1:length(statesToPlot)
    currState = statesToPlot(i);
    plot(timeGrid(1:end-1),exp(bd_fast_noCost(end:-1:2,currState) + theta), 'Color', [rval(i),gval,bval], 'LineWidth', 1.5)
end
plot(timeGrid(1:end-1),exp(bd_fast_noCost(end:-1:2,midstate) + theta), 'w--', 'LineWidth', 2)
plot(timeGrid(1:end-1),exp(bd_NL_noCost(end:-1:2,midstate) + theta), 'b--', 'LineWidth', 1.5)
grid on
xlabel('Time','FontSize',12)
ylabel('Exercise Price','FontSize',12)
title('\textbf{Fast Learning - no running costs}','FontSize',12)
ylim([20 60])

%% option value thru time  - FIGURE 6
startState = 16;
timesToPlot = [1 50 100 150 255];

f = figure;
set(f,'Position',[100 100 900 400]);

rval = linspace(0.1,0.95,length(timesToPlot));
bval = 0.2;
gval = 0.2;

% slow learning
subplot(1,2,1)
hold on;
for i = 1:length(timesToPlot)
    timeIdx = T*51 - timesToPlot(i) + 1;    
    plot(spot,squeeze(priceThruTime_slow(:,startState,timeIdx)), 'Color', [rval(i),gval,bval], 'LineWidth', 1)
    legendName{i} = ['$t$ = ' num2str(round(timesToPlot(i)/51))];
end
plot(spot,squeeze(priceThruTime_NL(:,startState,1)), 'b--', 'LineWidth', 1)
legendName{i+1} = 'No Learning';
xlabel('Spot Price','fontsize',12)
ylabel('Option Value at Mid-State','fontsize',12)
grid on
title(['\textbf{Slow Learning}'],'fontsize',12)
xlim([5 200])
legend(legendName,'Position',[0.36 0.018 0.118 0.281],'fontsize',10);

% fast learning
subplot(1,2,2)
hold on;
for i = 1:length(timesToPlot)
    timeIdx = T*51 - timesToPlot(i) + 1;    
    plot(spot,squeeze(priceThruTime_fast(:,startState,timeIdx)), 'Color', [rval(i),gval,bval], 'LineWidth', 1)
    legendName{i} = ['$t$ = ' num2str(round(timesToPlot(i)/51))];
end
plot(spot,squeeze(priceThruTime_NL(:,startState,1)), 'b--', 'LineWidth', 1)
legendName{i+1} = 'No Learning';
xlabel('Spot Price','fontsize',12)
ylabel('Option Value at Mid-State','fontsize',12)
grid on
title(['\textbf{Fast Learning}'],'fontsize',12)
xlim([5 200])
