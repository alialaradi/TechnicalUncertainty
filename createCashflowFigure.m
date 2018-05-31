%% CREATE ILLUSTRATIVE CASH FLOW PLOT (FIGURE 1)

%% PLOT PARAMETERS
% time parameters
t0 = 1;         % investment time
t1 = 4;         % extraction end
t2 = 5;         % horizon
epsilon = 0.5;  % time to extraction start

% extraction rate parameters
beta = 1;       
alpha = 9;

%% CREATE PLOT
% work with figure object 1
fig = figure(1);
clf(fig);
set(gca,'fontsize',18,'Position',[0.16 0.16 0.75 0.77]);

hold on

% define extraction rate curve - Equation (7)
u = [t0+epsilon : 0.1 : t1];
y = alpha*exp(-beta*(u-(t0+epsilon)));

% plot curve and fill area underneath
fill([u u(end) u(end:-1:1)], [y 0 zeros(1,length(u))],[0.8 0.8 0.8],'edgecolor','none');
plot(u, y,'-k','linewidth',2);

% insert vertical lines at investment time start and end of extraction
plot([1 1]*t0, [0 10], '--k');
plot([1 1]*(t0+epsilon), [0 10], '--k');
plot([1 1]*t1, [0 10], '--k');
xlim([0 t2]);

%% PLOT TEXT
% x and y labels
xlabel('Time','fontsize',24,'interpreter','latex');
ylabel('Extraction Rate','fontsize',24,'interpreter','latex');

% add text
text(t0-0.15,2, 'Investment made','rotation',90,'fontsize',14,'interpreter', 'latex');
text(t0+epsilon-0.15,2, 'Extraction begins','rotation',90,'fontsize',14,'interpreter', 'latex');
text(t1-0.15,2, 'Extraction ends','rotation',90,'fontsize',14,'interpreter', 'latex');
text((t0+epsilon+0.15),1.25, {'Total Volume','Extraced = $\gamma\;\vartheta$'},'fontsize',14,'interpreter','latex');

%% PLOT MODIFICATIONS
% remove x-ticks
set(gca, 'xtick',[t0 t0+epsilon+0.4 t1+0.7], 'xticklabel',{'','',''});

% re-label x-ticks manually 
xTicks = get(gca, 'xtick');
yTicks = get(gca,'ytick');
ax = axis; %Get left most x-position
minY = min(yTicks);
verticalOffset =0.4;

mylbl = {'$t_0$','$t_0+\epsilon$','$t_0+\epsilon+\Delta$'};
for xx = 1:length(xTicks)
%Create text box and set appropriate properties
     text(xTicks(xx), minY - verticalOffset, mylbl{xx},...
         'HorizontalAlignment','Right','interpreter', 'latex','fontsize',18);   
end

% remove y-ticks
set(gca, 'ytick',[]);