function Chaos;

% 3-variable Rossler model - chaos
% Didier Gonze
% 8/7/2008

clc;

%%%% Number of variable and initial conditions:

nbvar=17; 
for i=1:17
    xini(i)=1e4;
end

xini(1)=1e5;
xini(15)=1e6;

%%%% Time parameters:

trans=200000;
tend=202500;
tstep=1;

%%%% Range (for bifurcation diagram as a function of b):

b=0.2;   % (default value for chaos) 

bmin=0;
bmax=0.8;
bint=0.0025;



brange=[bmin bint bmax];



%%%% Task:

%integration(xini,trans,tend,tstep,b);
bifurcation(xini,trans,tend,tstep,brange);



%====================================================================
% Integration
%====================================================================

function output=integration(x0,trans,tend,tstep,b);

[t,x] = run(x0,trans,tend,tstep,b);


set(figure(1),'Position', [400 400 500 300]);  
clf;

plot(t,x(:,1:3));           
xlabel('Time','fontsize',18);
ylabel('Variables','fontsize',18);
xlim([0 tend]);
legend('X','Y','Z');


set(figure(2),'Position', [400 400 500 300]);  
clf;

plot3(x(:,1),x(:,2),x(:,3));
xlabel('X','fontsize',18);
ylabel('Y','fontsize',18);
zlabel('Z','fontsize',18);
box on;



%====================================================================
% Bifurcation
%====================================================================

function output=bifurcation(x0,trans,tend,tstep,range);

D=[];  % data (bifurcation diagram)

for b=range(1):range(2):range(3)

    fprintf('b=%g...\n',b);

    [t,x] = run(x0,trans,tend,tstep,b);
    
%     for i=2:length(x(:,2))-1
%         if((x(i,2)>x(i-1,2))&&(x(i,2)>x(i+1,2)))
%             D=[D; b x(i,2)]
%         end
%     end

D=[D; b max(x(:,3)+x(:,17)) min(x(:,3)+x(:,17))]
    
end

figure(3)
plot(D(:,1),D(:,2),'ro','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',1.5)
hold on;
plot(D(:,1),D(:,3),'ro','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',1.5)
xlabel('b','fontsize',18);
ylabel('max(X)','fontsize',18);
hold off;


% ===================================================================
% Run
% ===================================================================

function [t,x]=run(x0,trans,tend,tstep,b)

ttrans = [0:tstep:trans];
tspan = [0:tstep:tend];

option = [];
%option = odeset('RelTol', 1e-5);
%option=odeset('OutputS',[1:3],'OutputF','odeplot');

if trans > 0 
    [t x] = ode45(@dxdt,ttrans,x0,option,b);
    x0=x(end,:);
end

[t x] = ode45(@dxdt,tspan,x0,option,b);



% ===================================================================
% dxdt
% ===================================================================

function dy = dxdt(t,y,b)

%%% parameters
piH=26; 
piV=2600;
eta1=1.4;
eta2=200; 
muH=(1/65)/365; 
muV=26/365;
phi1=73/365;
phi2=73/365;
gamma1=52/365;
gamma2=52/365;
deltaH1=0.18/365;
deltaH2=0.18/365;
deltaV=0/365;
omega1=12/365;
omega2=12/365;
%b=0.3;
u2=0.3;
alphaH1=100/365;
alphaH2=25/365;
alphaV1=1000/365;   
alphaV2=800/365;

NH=y(1)+y(2)+y(3)+y(4)+y(5)+y(6)+y(7)+y(8)+y(9)+y(10)+y(11)+y(12)+y(13)+y(14);
NV=y(15)+y(16)+y(17);

lambdaH1 = alphaH1*(b*y(13)+y(3))/NH;
lambdaV1 = alphaV1*(y(16))/NH;

lambdaH2 = alphaH2*(eta2*y(7)+y(9))/NH;
lambdaV2 = alphaV2*(y(17))/NH;

dy(1)=piH-lambdaV1*y(1)-lambdaV2*y(1)-muH*y(1);

dy(2)=lambdaV1*y(1)-(phi1 + muH)*y(2);

dy(3)=phi1*y(2)-(gamma1 + muH+ deltaH1)*y(3);

dy(4)=gamma1*y(3)-(muH+omega1)*y(4);

dy(5)=omega1*y(4)-(b*lambdaV2+(1-b)+muH)*y(5);

dy(6)=b*lambdaV2*y(5)-(phi2+muH)*y(6);

dy(7)=phi2*y(6)-(gamma2 + muH+ deltaH2)*y(7);

dy(8)=lambdaV2*y(1)-(phi2 + muH)*y(8);

dy(9)=phi2*y(8)-(gamma2 + muH+ deltaH2)*y(9);

dy(10)=gamma2*y(9)-(muH+omega2)*y(10);

dy(11)=omega2*y(10)-(u2*lambdaV1+(1-u2)+muH)*y(11);

dy(12)=u2*lambdaV1*y(11)-(phi1+muH)*y(12);

dy(13)=phi1*y(12)-(gamma1 + muH+ deltaH1)*y(13);

dy(14)=gamma1*y(13) +gamma2*y(7) +(1-b)*y(5) +(1-u2)*y(11) -muH*y(14);

dy(15)=piV-lambdaH1*y(15)-lambdaH2*y(15)-muH*y(15);

dy(16)=lambdaH1*y(15)-(deltaV+muV)*y(16);

dy(17)=lambdaH2*y(15)-(deltaV+muV)*y(17);

%%% equations

dy=dy';


