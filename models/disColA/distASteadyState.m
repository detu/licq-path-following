function distASteadyState
%DISTASTEADYSTATE Summary of this function goes here
% 
% [OUTPUTARGS] = DISTASTEADYSTATE(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2016/06/22 21:05:37 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2016

% TO DO:
% - let use CasADi for computing Hessian !
% - continue with Greshgorin convexification 
format long;

%% parameter values
NT  = 41;                             % number of trays
LT  = 2.827;                          % Reflux
VB  = 3.454;                          % Boilup
F   = 1.2;                            % Feedrate
%zF  = 0.65;                           % Feed composition at CSTR 
zF  = 0.45;
qF  = 1.0;                            % Feed liquid fraction
dist.F  = F;
dist.NT = NT;
dist.zF = zF;
dist.qF = qF;

%% optimization preparation
% decision variables are states and controls
X = 0.5*ones(2*NT,1);
x = X(1:NT,1);                          % Liquid composition from btm to top
M = X(NT+1:2*NT,1);                     % Liquid hold up from btm to top
Uinit = [x;M;LT;VB];

% bound constraints
lb_u = [0.1; 0.1];
ub_u = [10; 4.008];
lb_x = zeros(2*NT,1);
ub_x = ones(2*NT,1);
lb   = [lb_x;lb_u];
ub   = [ub_x;ub_u];

% price setting
price.pf = 1; 
price.pV = 0.01; 
price.pB = 1; 
price.pD = 2;

%% call optimization
%options = optimset('Display','iter','TolFun',1e-6,'TolCon',1e-5,'Algorithm','active-set');
options = optimset('Display','iter','TolFun',1e-6,'TolCon',1e-5);

% options = optimoptions(@fmincon,'Algorithm','interior-point',...
%                        'Display','off','SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,...
%                        'HessianFcn',@(u,lambda)hessinterior(u,lambda,dist,price));

[xsol,objVal,exitFlag,output,lambda,~,Hess] = fmincon(@(u)objSS(u,dist,price), Uinit, [], [], [], [], lb, ub, @(u)constSS(u,dist), options);

%% Hessian regularization
import casadi.* 
% Symbolic primitives
x = {};
l = {};
for i=1:2*NT
   x{i} = SX.sym(['x_' num2str(i)],1);
   l{i} = SX.sym(['l_' num2str(i)],1);
end
u1  = SX.sym('u1');   % LT
u2  = SX.sym('u2');   % VB
% concatenate states and controls 
x   = vertcat(x{:});
x   = [x;u1;u2];
% extract Lagrange multiplier
l   = vertcat(l{:});
 
% define the dynamics as equality constraints
% eq1   = 1 - a1*exp(-1/x3)*x1^alpha - a2*exp(-delta/x3)*x1 - x1;
% eq2   = a1*exp(-1/x3)*x1^alpha - x2;
% eq3   = u - x3;
[obj,eq] = buildModelEq(x,dist,price);

% define Lagrangian
%L    = -x2 + l1*eq1 + l2*eq2 + l3*eq3;
L  = obj + l'*eq;

Lagr = Function('Lagr', {x,l}, {L}, char('x','l'), char('Lagr'));
Jac  = Function(Lagr.jacobian('x','Lagr'));
H    = Function(Lagr.hessian('x','Lagr'));
cons = Function('Const', {x}, {eq}, char('x'), char('cons'));
Jcon = Function(cons.jacobian('x','cons'));

Hx   = H(xsol,lambda.eqnonlin);
Hx   = full(Hx);
Jac  = Jcon(xsol);
Jac  = full(Jac);
rH   = null(Jac)'*Hx*null(Jac)

% Hitung Qmax untuk xsol ! 
% lalu baru compare dengan 
Qmax = zeros(84,1);
[Hxxl,Q,Qmax]   = Greshgorin(Hx,Qmax);
save Qmax.mat Qmax;

end

function [J,ceq] = buildModelEq(u,dist,price)
import casadi.* 
NT = dist.NT;
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
% objective function
F  = dist.F;    
J = price.pf*F + price.pV*u(2*NT+2) - price.pB*B - price.pD*D;

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
%i=1:NT-1;    y(i)=alpha*u(i)./(1+(alpha-1)*u(i));
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
%i=1:NT-1;    V(i)=VB*ones(1,NT-1);
%i=NF:NT-1;   V(i)=V(i) + (1-qF)*F;

V = vertcat(V{:});
for i=1:NT-1
    V(i) = VB;
    if i >= NF
        V(i) = V(i) + (1-qF)*F;
    end
end

% Liquid flows assuming linearized tray hydraulics with time constant taul
% Also includes coefficient lambda for effect of vapor flow ("K2-effect").
%i=2:NF;      L(i) = L0b + (M(i)-M0(i))./taul;
%i=2:NF;      L(i) = L0b + (u(NT+i,1)-M0(1,i))./taul;
%i=NF+1:NT-1; L(i) = L0  + (M(i)-M0(i))./taul;
%i=NF+1:NT-1; L(i) = L0  + (u(NT+i,1)-M0(1,i))./taul;

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
%i=2:NT-1;
%dMdt(i) = L(i+1)         - L(i)       + V(i-1)         - V(i);
%dMxdt(i)= L(i+1).*u(i+1) - L(i).*u(i) + V(i-1).*y(i-1) - V(i).*y(i);
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
%i=1:NT;
%ceq(i) = (dMxdt(i) - u(i).*dMdt(i) )./u(NT+i,1);
%ceq(i) = (dMxdt(i) - u(i).*dMdt(i) )./M(i);
%dxdt(i) = (dMxdt(i) - x(i).*dMdt(i) )./M(i);
for i=1:2*NT
    ceq{i}  = SX.sym(['ceq_' num2str(i)],1);
end
ceq  = vertcat(ceq{:});
for i=1:NT
    ceq(i) = (dMxdt(i) - u(i).*dMdt(i) ) / u(NT+i,1);
end

i=1:NT;
ceq(NT+i) = dMdt(i);
end

function [H,Q,Qm] = Greshgorin(H,Qm)
numH    = size(H,1);
Q       = zeros(numH,numH);
%Q       = eye(numH);
%delta   = 1e-2;
%delta   = 1;
delta   = 1e-6;
for i=1:numH  % iterate all row of Hessian
    sumRow = 0;
    for j=1:numH
        if j ~= i
            sumRow = sumRow + abs(H(i,j));
        end
    end
    %sumRow = sumRow - abs(H(i,i));
    
    if H(i,i) <= sumRow   % include equality 
        Q(i,i) = sumRow - H(i,i) + delta;
    end
end

% loop over Qm to obtain maximum number
for i=1:numH
    if Q(i,i) > Qm(i,1)
        Qm(i,1) = Q(i,i);
    end
end
end

function h = hessinterior(u,lambda,dist,price)
import casadi.* 

NT = dist.NT;
%NT = 41; % HARD-CODE !!!
x = {};
l = {};
for i=1:2*NT
   x{i} = SX.sym(['x_' num2str(i)],1);
   l{i} = SX.sym(['l_' num2str(i)],1);
end
u1  = SX.sym('u1');   % LT
u2  = SX.sym('u2');   % VB
% concatenate states and controls 
x   = vertcat(x{:});
x   = [x;u1;u2];
% extract Lagrange multiplier
l   = vertcat(l{:});

% inequality constraints
mu1  = SX.sym('mu1');   
mu2  = SX.sym('mu2');   
mu   = [mu1;mu2];
xB   = x(1); 
xD   = x(NT);
cin  = [xB-0.008;
        0.95-xD];

[obj,eq] = buildModelEq(x,dist,price);

% define Lagrangian
%L  = obj + l'*eq;
L   = obj + l'*eq + mu'*cin;

%Lagr = Function('Lagr', {x,l}, {L}, char('x','l'), char('Lagr'));
%H    = Function(Lagr.hessian('x','Lagr'));
Lagr = Function('Lagr', {x,l,mu}, {L}, char('x','l','mu'), char('Lagr'));
H    = Function(Lagr.hessian('x','Lagr'));

Hx   = H(u,lambda.eqnonlin,lambda.ineqnonlin);
h    = full(Hx);

end