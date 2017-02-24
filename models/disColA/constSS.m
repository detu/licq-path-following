function [cin, ceq, gradcin, gradceq] = constSS(u,dist)
%CONSTSS Summary of this function goes here
% 
% [OUTPUTARGS] = CONSTSS(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2016/06/23 14:19:15 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2016


NT = dist.NT;
xB = u(1); 
xD = u(NT);

%% inequality constraints
cin = [xB-0.008;
       0.95-xD];

%% equality constraint
% WRITE DOWN MASS BALANCE EQUATIONS in COLAMOD.M HERE !!!
%Ceq = [];

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

% P-Controllers for control of reboiler and condenser hold up.
% controller gains
KcB = 10;  
KcD = 10;         
% Nominal holdups - these are rather small
MDs = 0.5; 
MBs = 0.5;          
% Nominal flows
Ds = 0.5; 
Bs = 0.5;          
% Actual reboiler and condenser holdup
MB = u(NT+1);  
MD = u(2*NT); 
D  = Ds +(MD-MDs)*KcD;       % Distillate flow
B  = Bs +(MB-MBs)*KcB;       % Bottoms flow    


% Inputs and disturbances
%LT = U(1);                            % Reflux
LT = u(2*NT+1);
%VB = U(2);                            % Boilup
VB = u(2*NT+2);
%D  = U(3);                            % Distillate
%B  = U(4);                            % Bottoms

%F  = U(5);                            % Feedrate
F  = dist.F;
%zF = U(6);                            % Feed composition
zF = dist.zF;
%qF = U(7);                            % Feed liquid fraction
qF = dist.qF;

% THE MODEL

% Vapor-liquid equilibria
i=1:NT-1;    y(i)=alpha*u(i)./(1+(alpha-1)*u(i));

% Vapor Flows assuming constant molar flows
i=1:NT-1;    V(i)=VB*ones(1,NT-1);
i=NF:NT-1;   V(i)=V(i) + (1-qF)*F;

% Liquid flows assuming linearized tray hydraulics with time constant taul
% Also includes coefficient lambda for effect of vapor flow ("K2-effect").
%i=2:NF;      L(i) = L0b + (M(i)-M0(i))./taul;
%i=2:NF;      L(i) = L0b + (u(NT+i,1)-M0(1,i))./taul;
%i=NF+1:NT-1; L(i) = L0  + (M(i)-M0(i))./taul;
%i=NF+1:NT-1; L(i) = L0  + (u(NT+i,1)-M0(1,i))./taul;

L(NT)=LT;
for i=2:NT-1
    if i <=NF
        L(i) = L0b + (u(NT+i,1)-M0(1,i))./taul;
    else
        L(i) = L0  + (u(NT+i,1)-M0(1,i))./taul;
    end
end



% Time derivatives from  material balances for 
% 1) total holdup and 2) component holdup

% Column
%i=2:NT-1;
%dMdt(i) = L(i+1)         - L(i)       + V(i-1)         - V(i);
%dMxdt(i)= L(i+1).*u(i+1) - L(i).*u(i) + V(i-1).*y(i-1) - V(i).*y(i);
for i=2:NT-1
    dMdt(i) = L(i+1)         - L(i)       + V(i-1)         - V(i);
    dMxdt(i)= L(i+1).*u(i+1,1) - L(i).*u(i,1) + V(i-1).*y(i-1) - V(i).*y(i);
end


% Correction for feed at the feed stage
% The feed is assumed to be mixed into the feed stage
dMdt(NF) = dMdt(NF)  + F;
dMxdt(NF)= dMxdt(NF) + F*zF;

% Reboiler (assumed to be an equilibrium stage)
dMdt(1) = L(2)      - V(1)      - B;
dMxdt(1)= L(2)*u(2) - V(1)*y(1) - B*u(1);

% Total condenser (no equilibrium stage)
dMdt(NT) = V(NT-1)         - LT       - D;
dMxdt(NT)= V(NT-1)*y(NT-1) - LT*u(NT) - D*u(NT);

% Compute the derivative for the mole fractions from d(Mx) = x dM + M dx
%i=1:NT;
%ceq(i) = (dMxdt(i) - u(i).*dMdt(i) )./u(NT+i,1);
%ceq(i) = (dMxdt(i) - u(i).*dMdt(i) )./M(i);
%dxdt(i) = (dMxdt(i) - x(i).*dMdt(i) )./M(i);
for i=1:NT
    ceq(i) = (dMxdt(i) - u(i).*dMdt(i) )./u(NT+i,1);
end

i=1:NT;
ceq(NT+i) = dMdt(i);

if nargout > 2
    % inequality constraints
    gradcin       = zeros(2*NT+2,2);
    gradcin(1,1)  = 1;
    gradcin(NT,2) = -1;
    % equality constraints 
    gradceq     = computeJacEq(u,dist);
end

end

function gradceq = computeJacEq(u,dist)
import casadi.*
NT = dist.NT;
% Symbolic primitives
x = {};
for i=1:2*NT
   x{i} = SX.sym(['x_' num2str(i)],1);
end
u1  = SX.sym('u1');   % LT
u2  = SX.sym('u2');   % VB
% concatenate states and controls 
x   = vertcat(x{:});
x   = [x;u1;u2];

% the dynamics
NT = dist.NT;
NF = 21;            % Location of feed stage (stages are counted from the bottom):
alpha = 1.5;        % Relative volatility
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
%lambda = 0;		% Effect of vapor flow on liquid flow ("K2-effect")

% P-Controllers for control of reboiler and condenser hold up.
% controller gains
KcB = 10;  
KcD = 10;         
% Nominal holdups - these are rather small
MDs = 0.5; 
MBs = 0.5;          
% Nominal flows
Ds = 0.5; 
Bs = 0.5;          
% Actual reboiler and condenser holdup
MB = u(NT+1);  
MD = u(2*NT); 
D  = Ds +(MD-MDs)*KcD;       % Distillate flow
B  = Bs +(MB-MBs)*KcB;       % Bottoms flow    

% Inputs and disturbances
LT = u(2*NT+1);
VB = u(2*NT+2);
F  = dist.F;
zF = dist.zF;
qF = dist.qF;

% THE MODEL

% Vapor-liquid equilibria
for i=1:NT-1
   y{i}  = SX.sym(['y_' num2str(i)],1);
   V{i}  = SX.sym(['V_' num2str(i)],1);
   L{i}  = SX.sym(['L_' num2str(i)],1);
   dMdt{i}  = SX.sym(['dMdt_' num2str(i)],1);
   dMxdt{i} = SX.sym(['dMxdt_' num2str(i)],1);
end
L{NT}     = SX.sym(['L_' num2str(NT)],1);
dMdt{NT}  = SX.sym(['dMdt_' num2str(NT)],1);
dMxdt{NT} = SX.sym(['dMxdt_' num2str(NT)],1);

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
L(NT)=LT;
for i=2:NT-1
    if i <=NF
        L(i) = L0b + (u(NT+i,1)-M0(1,i))./taul;
    else
        L(i) = L0  + (u(NT+i,1)-M0(1,i))./taul;
    end
end

% Time derivatives from  material balances for 
% 1) total holdup and 2) component holdup

% Column
dMdt  = vertcat(dMdt{:});
dMxdt = vertcat(dMxdt{:});
for i=2:NT-1
    dMdt(i) = L(i+1)         - L(i)       + V(i-1)         - V(i);
    dMxdt(i)= L(i+1).*u(i+1,1) - L(i).*u(i,1) + V(i-1).*y(i-1) - V(i).*y(i);
end


% Correction for feed at the feed stage
% The feed is assumed to be mixed into the feed stage
dMdt(NF) = dMdt(NF)  + F;
dMxdt(NF)= dMxdt(NF) + F*zF;

% Reboiler (assumed to be an equilibrium stage)
dMdt(1) = L(2)      - V(1)      - B;
dMxdt(1)= L(2)*u(2) - V(1)*y(1) - B*u(1);

% Total condenser (no equilibrium stage)
dMdt(NT) = V(NT-1)         - LT       - D;
dMxdt(NT)= V(NT-1)*y(NT-1) - LT*u(NT) - D*u(NT);

% Compute the derivative for the mole fractions from d(Mx) = x dM + M dx
for i=1:2*NT
    ceq{i}  = SX.sym(['ceq_' num2str(i)],1);
end
ceq  = vertcat(ceq{:});
for i=1:NT
    ceq(i) = (dMxdt(i) - u(i).*dMdt(i) ) / u(NT+i,1);
end

i=1:NT;
ceq(NT+i) = dMdt(i);

cons    = Function('Const', {x}, {ceq}, char('x'), char('cons'));
Jcon    = Function(cons.jacobian('x','cons'));
Jac     = Jcon(u);
gradceq = full(Jac);
gradceq = gradceq';
end
