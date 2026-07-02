% --------------------- Medium / model parameters ---------------------
clc;
mu  = 25e9;   
nu  = 0.25;
nuu = 0.30;
b   = 0.75;   % Biot B 
d   = 0.2;    % diffusivity D
q0  = 5e-3;    % m^2/s


stage  = 18;      % so total stages = stage+1 = 19 
ITime  = 3;       % hours injection
PTime  = 15;      % hours pause
Period = ITime + PTime;  % hours per stage
Ho  = 3600;
Day = 24*Ho;
Month = 30*Day;
% Total duration (hours): 19 stages * 18 h = 342 h = 14.25 days
Tend_hours = (stage+1)*Period;           % 342 h
Tend_sec   = Tend_hours * Ho;
% Time vector for evaluation
dt_hour = 0.5; 
tSec    = (0 : dt_hour*Ho : Tend_sec).';  % seconds
tDay    = tSec / Day;

% --------------------- Injection migration path ---------------------
x0 = 0;  
y0 = 0;
dDist = 50;                             % m per stage
xi = x0 : dDist : x0 + stage*dDist;   % length stage+1
yi = y0 : -dDist : y0 - stage*dDist;  % length stage+1

t1 = (0:Period:Period*stage);     % start times of each stage (hours)
t2 = t1 + ITime;                  % end of injection within each stage (hours)
% --------------------- CFS / fault parameters ---------------------
deg = 10;
theta = pi/2 - deg2rad(deg);
FRL = -1;
Friction = 0.6;

InjParam=zeros(stage+1,5);
for i=1:stage+1
    InjParam(i,:) = [ ...
        xi(i), yi(i), ...                     % InjX, InjY
        (t1(i)*Ho)/Month, ...            % start [month]
        (t2(i)*Ho)/Month, ...            % end   [month]
        q0];                                  % rate
end


% --------------------- Observation points ---------------------
obs = [100 100; 100 200];   % [x y] in meters
nObs = size(obs,1);
PorePressure = nan(numel(tSec), nObs);   % MPa
ShearStress  = nan(numel(tSec), nObs);   % MPa
NormalStress = nan(numel(tSec), nObs);   % MPa
CFS          = nan(numel(tSec), nObs);   % MPa

SigmaXX = nan(numel(tSec), nObs);        % optional, MPa
SigmaYY = nan(numel(tSec), nObs);        % optional, MPa
SigmaXY = nan(numel(tSec), nObs);        % optional, MPa

for j = 1:nObs
    x = obs(j,1); 
    y = obs(j,2);
    for it = 1:numel(tSec)
        [P, S, cfsVal, tauVal, sigmaVal] = IM_InjectionHistoryNew2d( ...
            x, y, tSec(it), InjParam, mu, nu, nuu, b, d, ...
            theta, FRL, Friction);

        PorePressure(it,j) = P / 1e6;          % Pa -> MPa
        ShearStress(it,j)  = tauVal / 1e6;     % Pa -> MPa
        NormalStress(it,j) = sigmaVal / 1e6;   % Pa -> MPa
        CFS(it,j)          = cfsVal / 1e6;     % Pa -> MPa

        SigmaXX(it,j) = S(1,1) / 1e6;
        SigmaYY(it,j) = S(2,2) / 1e6;
        SigmaXY(it,j) = S(1,2) / 1e6;
    end
end


%%
% --------------------- Plot CFS / shear / normal / pore ---------------------
figure('Color','w');
tiledlayout(nObs,1,'Padding','compact','TileSpacing','compact');
colsCFS    = [0.0314 0.7804 0.4314];
colsTau    = [0.9098 0.0627 0.4039];
colsNormal = [252 151 13] ./ 256;
colsPore   = [0.0745 0.6235 1.0000];

for j = 1:nObs
    ax = nexttile;
    hold(ax,'on');
    g1 = plot(ax, tDay, CFS(:,j), ...
        'Color', colsCFS, 'LineWidth', 3.5, ...
        'DisplayName', '\DeltaCFS');

    g2 = plot(ax, tDay, ShearStress(:,j), ...
        'Color', colsTau, 'LineWidth', 2.5, ...
        'DisplayName', 'Shear');

    g3 = plot(ax, tDay, Friction * NormalStress(:,j), ...
        'Color', colsNormal, 'LineWidth', 2.5, ...
        'DisplayName', 'f*Normal');

    g4 = plot(ax, tDay, Friction * PorePressure(:,j), ...
        'Color', colsPore, 'LineWidth', 2.5, ...
        'DisplayName', 'f*Pressure');
    yl = ylim(ax);
    ylim(ax, yl);
    ylabel(ax, 'Stress (MPa)');
    title(ax, sprintf('$(x_i,y_i)=(%d,%d)$', obs(j,1), obs(j,2)), ...
        'Interpreter','latex');
%     title(sprintf('$(x_i,y_i)=(%d,%d), q_0=%.3f, D=%.1f, Deg=%d$', ...
%         obs(j,1), obs(j,2), q0,D,deg),'Interpreter','latex');
    if j == 1
        legend(ax, [g1 g2 g3 g4], ...
            'Location','best', ...
            'FontSize',16, ...
            'Box','off','Orientation','horizontal');
    end
    if j < nObs
        ax.XTickLabel = [];
    else
        xlabel(ax, 'Time (day)');
    end
    xlim(ax, [0, Tend_sec/Day]);
    Fun_defaultAxes;
end
set(gcf,'Position',[250,80,900,700]);


%%
figure;
tiledlayout(1,1,"TileSpacing","compact","Padding","compact");
nexttile;hold on;
scatter(xi,yi,190,1:19,"filled",'Marker','square','LineWidth',0.6,'MarkerEdgeColor',[.2 .2 .2]); 
for i=1:size(obs,1)
    scatter(obs(i,1),obs(i,2),300,"filled",'Marker','pentagram','LineWidth',0.6, ...
        'MarkerEdgeColor',[.2 .2 .2],'MarkerFaceColor',[0.9098    0.0627    0.4039]);
end
% plot(xi,yi,'-','Marker','square','LineWidth',0.9,'Color',[.2 .2 .2],'LineStyle','none', ...
%     'MarkerSize',8,'MarkerFaceColor',[.8 .8 .8]); 
plot(xi,yi,'-','Marker','x','LineWidth',0.5,'Color',[.2 .2 .2],'LineStyle','none', ...
    'MarkerSize',8,'MarkerFaceColor',[.8 .8 .8]); 
% xticks(0:400:800);
% yticks(-800:400:0);
xlabel('$X_0 (m)$','Interpreter','latex'); 
ylabel('$Y_0 (m)$','Interpreter','latex');
c = colorbar('Location', 'eastoutside');
c.LineWidth=0.8;
c.Label.FontSize=13;
c.Label.FontWeight = 'normal';
clim([1 19])
c.Ticks=([1:4:19]);
c.TickLabels = arrayfun(@(x) sprintf('%.0f', x), c.Ticks, 'UniformOutput', false);
Fun_defaultAxes;
cmp=slanCM('viridis',21);
cmp=cmp(1:19,:);
colormap(cmp);
set(gcf,'Position',[250,80,590,500]);