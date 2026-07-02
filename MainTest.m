% =========================================================================
% Parameters:
D=1;
B=0.75;
nu=0.25;
nuu=0.3;
q0=5e-3;
G=25e9;
Xi=1e2;
Yi=1e2;
X0=0;
Y0=0;
deg=10;
theta=pi/2-deg2rad(deg);
FRL=-1;  % Lateral direction
Friction=0.6;
% Set the time sequence; Unit: Hour
Ho=3600; % one hour
tt1=3;   % Inject time
tt2=15;  % Pause time

% =========================================================================
% Pore-pressure
t1=0*Ho:Ho/60:tt1*Ho; % Every 1 minite
N1=length(t1);
p1=zeros(N1,1);
for i=1:length(t1)
    [p1(i)]=IM_Pressure2d(Xi,Yi,X0,Y0,t1(i),t1(1),q0,G,nu,nuu,B,D);
end
p(1)=0; % Ensure continuity
tend=t1(end);

t2=tend+(0:Ho/60*5:tt2*Ho); % Every 5 minite
N2=length(t2);
p11=zeros(N2,1);
p22=zeros(N2,1);
for i=1:length(t2)
    % Notice the start time are different here
    [p11(i)]=IM_Pressure2d(Xi,Yi,X0,Y0,t2(i),t1(1),q0,G,nu,nuu,B,D); % Hold 
    [p22(i)]=IM_Pressure2d(Xi,Yi,X0,Y0,t2(i),tend,-q0,G,nu,nuu,B,D); % Reverse
end
p22(1)=0;       % Ensure continuity
p2=p11+p22;     % Stack the (+,-) Pressure (i.e. q=0)




% =========================================================================
% CFS Change
sxx0=zeros(N1,1);
syy0=zeros(N1,1);
sxy0=zeros(N1,1);
tau0=zeros(N1,1);
sigma0=zeros(N1,1);
for i=1:length(t1)
    [s]=IM_Sigma2d(Xi,Yi,X0,Y0,t1(i),t1(1),q0,G,nu,nuu,B,D);
    sxx0(i)=s(1,1);
    syy0(i)=s(2,2);
    sxy0(i)=s(1,2);
    tau0(i)=FRL*(((-sxx0(i)+syy0(i))/2)*sin(2*theta)+sxy0(i)*cos(2*theta));
    sigma0(i)=+((sxx0(i)+syy0(i))/2-((sxx0(i)-syy0(i))/2)*cos(2*theta)-sxy0(i)*sin(2*theta));
end
tau0(1)=0;
tend=t1(end);
sxx1=zeros(N2,1);
syy1=zeros(N2,1);
sxy1=zeros(N2,1);
sxx2=zeros(N2,1);
syy2=zeros(N2,1);
sxy2=zeros(N2,1);
tau1=zeros(N2,1);
tau2=zeros(N2,1);
sigma1=zeros(N2,1);
sigma2=zeros(N2,1);
for i=1:length(t2)
    % Tensor to single force-value
    [s01]=IM_Sigma2d(Xi,Yi,X0,Y0,t2(i),t1(1),q0,G,nu,nuu,B,D);
    [s02]=IM_Sigma2d(Xi,Yi,X0,Y0,t2(i),tend,-q0,G,nu,nuu,B,D);
    sxx1(i)=s01(1,1);
    syy1(i)=s01(2,2);
    sxy1(i)=s01(1,2);
    sxx2(i)=s02(1,1);
    syy2(i)=s02(2,2);
    sxy2(i)=s02(1,2);
    % Compression negative, tension positive
    tau1(i)=FRL*(((-sxx1(i)+syy1(i))/2)*sin(2*theta)+sxy1(i)*cos(2*theta));
    tau2(i)=FRL*(((-sxx2(i)+syy2(i))/2)*sin(2*theta)+sxy2(i)*cos(2*theta));
    sigma1(i)=+((sxx1(i)+syy1(i))/2-((sxx1(i)-syy1(i))/2)*cos(2*theta)-sxy1(i)*sin(2*theta));
    sigma2(i)=+((sxx2(i)+syy2(i))/2-((sxx2(i)-syy2(i))/2)*cos(2*theta)-sxy2(i)*sin(2*theta));
end
tau2(1)=0;
% Ontained the conbined value: tau and sigma
tauNew=tau1+tau2;
sigmaNew=sigma1+sigma2;
% Ontained the CFS, using empirical friction constant
CFS0=tau0+Friction*(sigma0+p1);
CFSNew=tauNew+Friction*(sigmaNew+p2);


% Plots
figure;
tiledlayout(1,1,'TileSpacing','loose','Padding','compact');
nexttile;hold on;
ll=xline(tend./Ho,'LineWidth',2,'Color','k','LineStyle','--');
gcfs=plot(t1./Ho,CFS0./(1e6),'Color',[ 0.0314    0.7804    0.4314],'LineWidth',5);
plot(t2./Ho,CFSNew./(1e6),'Color',[ 0.0314    0.7804    0.4314],'LineWidth',5);
gtau=plot(t1./Ho,tau0./(1e6),'Color',[0.9098    0.0627    0.4039],'LineWidth',3);
plot(t2./Ho,tauNew./(1e6),'Color',[0.9098    0.0627    0.4039],'LineWidth',3);
gsigma=plot(t1./Ho,Friction*sigma0./(1e6),'Color',[252 151 13]./256,'LineWidth',3);
plot(t2./Ho,Friction*sigmaNew./(1e6),'Color',[252 151 13]./256,'LineWidth',3);
gpore=plot(t1./Ho,Friction*p1./(1e6),'Color',[0.0745    0.6235    1],'LineWidth',3);
plot(t2./Ho,Friction*p2./(1e6),'Color',[0.0745    0.6235    1],'LineWidth',3);
legend([ll,gcfs,gtau,gsigma,gpore],{'Pause Tick','\Delta{CFS}','Shear','f*Normal','f*Pore'},...
    'fontsize',24,'location','NE','box','on');
xlabel('Time (hours)');
ylabel('\DeltaCFS (MPa)');
Fun_defaultAxes;
set(gcf,'position',[300,50,800,700]);
title(sprintf('$(x_i,y_i)=(%d,%d), q_0=%.3f, D=%.1f, Deg=%d$', ...
    Xi,Yi,q0,D,deg),'Interpreter','latex');
xlim('tight');