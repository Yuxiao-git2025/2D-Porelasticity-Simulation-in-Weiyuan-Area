% --------------------- Medium / model parameters ---------------------
clc;
mu  = 25e9;   
nu  = 0.25;
nuu = 0.30;
b   = 0.75;
d   = 0.2;
q0  = 5e-3;

stage  = 18;
ITime  = 3;
PTime  = 15;
Period = ITime + PTime;

Ho = 3600;
Day = 24 * Ho;
Month = 30 * Day;

Tend_hours = (stage + 1) * Period;
Tend_sec = Tend_hours * Ho;

% --------------------- Time vector ---------------------
dt_day = 0.1;
tDay = (0 : dt_day : Tend_sec / Day).';
tSec = tDay * Day;
nTime = numel(tSec);

% --------------------- Injection migration path ---------------------
x0 = 0;
y0 = 0;
dDist = 50;

xi = x0 : dDist : x0 + stage * dDist;
yi = y0 : -dDist : y0 - stage * dDist;

t1 = 0 : Period : Period * stage;
t2 = t1 + ITime;

% --------------------- CFS / fault parameters ---------------------
deg = 10;
theta = pi/2 - deg2rad(deg);
FRL = -1;
Friction = 0.6;

% --------------------- Injection parameters ---------------------
InjParam = zeros(stage + 1, 5);

for i = 1:stage + 1
    InjParam(i,:) = [ ...
        xi(i), yi(i), ...
        (t1(i) * Ho) / Month, ...
        (t2(i) * Ho) / Month, ...
        q0];
end

% --------------------- 2D grid ---------------------
xVec = -205 : 10 : 200; xVec(xVec==0)=1e-1; % Protection
yVec = -205 : 10 : 800; yVec(yVec==0)=1e-1; % Protection

[XGrid, YGrid] = meshgrid(xVec, yVec);

nx = numel(xVec);
ny = numel(yVec);

% --------------------- Allocate fields ---------------------
PorePressure = nan(ny, nx, nTime);
ShearStress  = nan(ny, nx, nTime);
NormalStress = nan(ny, nx, nTime);
CFS          = nan(ny, nx, nTime);

% --------------------- Calculate fields ---------------------
tic;
for it = 1:nTime
    fprintf('Time step %d / %d, t = %.2f day\n', it, nTime, tDay(it));
    for iy = 1:ny
        for ix = 1:nx

            x = XGrid(iy, ix);
            y = YGrid(iy, ix);

            [P, ~, cfsVal, tauVal, sigmaVal] = IM_InjectionHistoryNew2d( ...
                x, y, tSec(it), InjParam, mu, nu, nuu, b, d, ...
                theta, FRL, Friction);

            PorePressure(iy, ix, it) = P / 1e6;
            ShearStress(iy, ix, it)  = tauVal / 1e6;
            NormalStress(iy, ix, it) = sigmaVal / 1e6;
            CFS(iy, ix, it)          = cfsVal / 1e6;
        end
    end
end
toc;

%% Seismicity
% [R, rateOut]=IM_CalSeismicity(tDay, CFS, ...
%     'a', 0.003, ...
%     'SigmaEff', 16.67, ...
%     'TauDot0', 0.001 / 365.25, ...
%     'R0', 1);


%% Save data
save('Injection2d_Weiyuan.mat', ...
    'xVec', 'yVec', 'XGrid', 'YGrid', 'tDay', 'tSec', ...
    'PorePressure', 'ShearStress', 'NormalStress', 'CFS', ...
    'xi', 'yi', 't1', 't2', 'InjParam', ...
    'mu', 'nu', 'nuu', 'b', 'd', 'q0', ...
    'theta', 'FRL', 'Friction', ...
    '-v7.3');

%%
% Continue from MainGridSearch Script
% --------------------- Plot 2D snapshots of four quantities ---------------------
itPlot = 51;
figure('Color','w');
arr=tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
dataList = {
    CFS(:,:,itPlot),          '\DeltaCFS (MPa)';
    ShearStress(:,:,itPlot),  'Shear (MPa)';
    NormalStress(:,:,itPlot), 'Normal (MPa)';
    PorePressure(:,:,itPlot), 'Pore-Pressure (MPa)';
    };
colormap(flip(slanCM('RdBu',200)));
for k = 1:4
    ax = nexttile;
    hold(ax,'on');
    imagesc(ax, xVec, yVec, dataList{k,1});
    set(ax,'YDir','normal');
    axis(ax,'tight');
    cb = colorbar(ax);
    cb.Label.String = dataList{k,2};
    xticks([-200 0 200]);
    yticks([-200 0 800]);
    Fun_Decorat;
end
xlabel(arr, '$XDist$ (m)', 'Interpreter','latex','FontSize',22,'FontName','Times New Roman');
ylabel(arr, '$YDist$ (m)', 'Interpreter','latex','FontSize',22,'FontName','Times New Roman');
title(arr, sprintf('Time Tick=%.2f (day)',tDay(itPlot)),'FontSize',22,'FontName','Times New Roman');
set(gcf,'Position',[200,50,1000,750]);