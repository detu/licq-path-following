function [JEqmodel] = CstrDistAEqModel(Xk,Uk)
%CSTRDISTAEQMODEL Summary of this function goes here
% 
% Compute equality constraint values for rotated cost function 
%
% [OUTPUTARGS] = CSTRDISTAEQMODEL(INPUTARGS) Explain usage here
%  
%
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2016/07/19 15:34:33 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2016

import casadi.*

u = [Xk;Uk];
NT = 41;
% Location of feed stage (stages are counted from the bottom):
NF = 21;
% Relative volatility
alpha = 1.5;
% Nominal liquid holdups
M0(1)=0.5;      	% Nominal reboiler holdup (kmol)
i=2:NT-1; M0(i)=0.5*ones(1,NT-2);% Nominal stage (tray) holdups (kmol)
M0(NT)=0.5;      	% Nominal condenser holdup (kmol)
% Data for linearized liquid flow dynamics (does not apply to reboiler and condenser):
taul = 0.063;     	% time constant for liquid dynamics (min)
F0   = 1;	 	    % Nominal feed rate (kmol/min) 
qF0  = 1; 		    % Nominal fraction of liquid in feed 
L0   = 2.70629;     % Nominal reflux flow (from steady-state data)
L0b  = L0 + qF0*F0;	% Nominal liquid flow below feed (kmol/min)
%lambda = 0;		    % Effect of vapor flow on liquid flow ("K2-effect")


% Inputs and disturbances
LT   = u(2*NT+3);               % Reflux
VB   = u(2*NT+4);               % Boilup
F    = u(2*NT+5);               % Feedrate
D    = u(2*NT+6);               % Distillate
B    = u(2*NT+7);               % Bottoms flow 
F_0  = 0.3;    
zF   = 0.45;                    % Feed composition                         
qF   = 1.0;                     % Feed liquid fraction


% Vapor-liquid equilibria
%i=1:NT-1;    y(i)=alpha*u(i)./(1+(alpha-1)*u(i));
 
for i=1:NT-1
   y{i}  = MX.sym(['y_' num2str(i)],1);
   V{i}  = MX.sym(['V_' num2str(i)],1);
   L{i}  = MX.sym(['L_' num2str(i)],1);
   dMdt{i}  = MX.sym(['dMdt_' num2str(i)],1);
   dMxdt{i} = MX.sym(['dMxdt_' num2str(i)],1);
   dxdt{i}  = MX.sym(['dxdt_' num2str(i)],1);
end
L{NT}     = MX.sym(['L_' num2str(NT)],1);
dMdt{NT}  = MX.sym(['dMdt_' num2str(NT)],1);
dMxdt{NT} = MX.sym(['dMxdt_' num2str(NT)],1);
dxdt{NT}  = MX.sym(['dxdt_' num2str(NT)],1);
dMdt{NT+1}  = MX.sym(['dMdt_' num2str(NT+1)],1);
dMxdt{NT+1} = MX.sym(['dMxdt_' num2str(NT+1)],1);
dxdt{NT+1}  = MX.sym(['dxdt_' num2str(NT+1)],1);

% % collect equality constraints
% ceq = zeros(2*NT+2,1);
% % prepare necessary variables
% y     = zeros(NT-1,1);
% V     = zeros(NT-1,1);
% L     = zeros(NT,1);
% dMdt  = zeros(NT+1,1);
% dMxdt = zeros(NT+1,1);


y = vertcat(y{:});
for i=1:NT-1
    y(i)=alpha*u(i)/(1+(alpha-1)*u(i));
end

% Vapor Flows assuming constant molar flows
V = vertcat(V{:});
for i=1:NT-1
    V(i) = VB;
    if i >= NF
        V(i) = V(i) + (1-qF)*F;
    end
end

% Liquid flows assuming linearized tray hydraulics with time constant taul
% Also includes coefficient lambda for effect of vapor flow ("K2-effect").
L = vertcat(L{:});
L(NT) = LT;
for i=2:NT-1
    if i <=NF
        L(i) = L0b + (u(NT+i)-M0(1,i))./taul;
    else
        L(i) = L0  + (u(NT+i)-M0(1,i))./taul;
    end
end



% Time derivatives from  material balances for 
% 1) total holdup and 2) component holdup

% Column
dMdt  = vertcat(dMdt{:});
dMxdt = vertcat(dMxdt{:});
for i=2:NT-1
    dMdt(i) = L(i+1)         - L(i)       + V(i-1)         - V(i);
    dMxdt(i)= L(i+1).*u(i+1) - L(i).*u(i) + V(i-1).*y(i-1) - V(i).*y(i);
end

% Correction for feed at the feed stage
% The feed is assumed to be mixed into the feed stage
dMdt(NF) = dMdt(NF)  + F;
dMxdt(NF)= dMxdt(NF) + F*u(NT+1);

% Reboiler (assumed to be an equilibrium stage)
dMdt(1) = L(2)      - V(1)      - B;
dMxdt(1)= L(2)*u(2) - V(1)*y(1) - B*u(1);

% Total condenser (no equilibrium stage)
dMdt(NT) = V(NT-1)         - LT       - D;
dMxdt(NT)= V(NT-1)*y(NT-1) - LT*u(NT) - D*u(NT);

% CSTR model
k1          = 34.1/60.0;
dMdt(NT+1)  = F_0 + D - F;
dMxdt(NT+1) = F_0*zF + D*u(NT) - F*u(NT+1) - k1*u(2*NT+2)*u(NT+1);

dxdt  = vertcat(dxdt{:});
for i=1:NT+1
    dxdt(i) = (dMxdt(i) - u(i).*dMdt(i) ) / u(NT+1+i,1);
end

ceq = [dxdt;dMdt];

% load Lagrange multipliers from steady-state optimization
load LamdaCstrDist.mat; % lamda

JEqmodel = lamda'*ceq;
end
