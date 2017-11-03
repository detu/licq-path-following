function distACstrSS
%DISTACSTRSS Summary of this function goes here
%
% Steady-state CSTR + Distillation Column A
% 
% [OUTPUTARGS] = DISTACSTRSS(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2016/06/24 14:30:29 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2016


format long;
import casadi.*

%% parameter values
NT  = 41;                            % number of trays
LT  = 2.827;                         % Reflux
VB  = 3.454;                         % Boilup
F   = 1.0;                           % Feedrate
zF  = 1.0;                           % Feed composition at CSTR 
D   = 0.5;                           % Distillate flow
B   = 0.5;                           % Bottoms flow 
qF  = 1.0;                           % Feed liquid fraction

% create a struct to make them nice
dist.F_0  = 0.3;
dist.NT   = NT;
dist.zF   = zF;
dist.qF   = qF;

% price setting
price.pf = 1; 
price.pV = 0.02;
price.pB = 2; 
price.pD = 0;

% Symbolic primitives
x = {};
l = {};
for i=1:2*NT+2       % +2 due to CSTR
   x{i} = SX.sym(['x_' num2str(i)],1);
   l{i} = SX.sym(['l_' num2str(i)],1);
end
u1  = SX.sym('u1');   % LT
u2  = SX.sym('u2');   % VB
u3  = SX.sym('u3');   % F
u4  = SX.sym('u4');   % D
u5  = SX.sym('u5');   % B
% concatenate states and controls 
x   = vertcat(x{:});
x   = [x;u1;u2;u3;u4;u5];

% decision variables are states and controls
Xinit = 0.5*ones(84,1);
Uinit = [Xinit;LT;VB;F;D;B];

% define the dynamics as equality constraints and additional inequality
% constraints (lbx, ubx, lbg, and ubg)
[obj,eq, lbx, ubx, lbg, ubg] = buildModelEq(x,dist,price);

prob   = struct('f', obj, 'x', x, 'g', eq);
options = struct;
%options.ipopt.tol       = 1e-12;
%options.acceptable_compl_inf_tol    = 1e-6;
solver = nlpsol('solver', 'ipopt', prob, options);

% Solve the NLP
startnlp = tic;
sol   = solver('x0', Uinit, 'lbx', lbx, 'ubx', ubx, 'lbg', lbg, 'ubg', ubg);
elapsednlp = toc(startnlp);
fprintf('IPOPT solver runtime = %f\n',elapsednlp);

u      = full(sol.x);
lamda  = full(sol.lam_g);
lamda(43:end) = -1*lamda(43:end);

Xinit  = u;
save CstrDistXinit.mat Xinit;
save LamdaCstrDist.mat lamda;


%% Compute Hessian and perform Greshgorin convexification
xsol = u;
lambda.eqnonlin = lamda;
% extract Lagrange multiplier
l   = vertcat(l{:});
L  = obj + l'*eq;

Lagr = Function('Lagr', {x,l}, {L}, char('x','l'), char('Lagr'));
H    = Function(Lagr.hessian('x','Lagr'));
cons = Function('Const', {x}, {eq}, char('x'), char('cons'));
Jcon = Function(cons.jacobian('x','cons'));

eqVal = cons(xsol);
Hx   = H(xsol,lambda.eqnonlin);
Hx   = full(Hx);
Jac  = Jcon(xsol);
Jac  = full(Jac);
% display nullspace of the constraint and its eigenvalue
rH   = null(Jac)'*Hx*null(Jac)
erH  = eig(rH)                   


[Hxxl,Qmax]   = Greshgorin(Hx);
save Qmax.mat Qmax;

% check at initial point for optimization
xstat  = Xinit(1:84);
u0     = [2.5; 3.5; 0.6; 0.5; 0.5];
xeval  = [xstat;u0];
Jeval  = Jcon(xeval);
Jeval  = full(Jeval);
Hxxl   = H(xeval,lambda.eqnonlin);
Hxxl   = full(Hxxl);
Hconv  = Hxxl + diag(Qmax);
rHe    = null(Jeval)'*Hconv*null(Jeval);

keyboard;

end

function [J, ceq, lbx, ubx, lbg, ubg] = buildModelEq(u,dist,price)
import casadi.* 
NT = dist.NT;
% Location of feed stage (stages are counted from the bottom):
NF = 21;
% Relative volatility
alpha = 1.5;
% Nominal liquid holdups
Muw   = 0.5;
% Data for linearized liquid flow dynamics (does not apply to reboiler and condenser):
taul = 0.063;     	% time constant for liquid dynamics (min)
F0   = 1;	 	    % Nominal feed rate (kmol/min) 
qF0  = 1; 		    % Nominal fraction of liquid in feed 
L0   = 2.70629;     % Nominal reflux flow (from steady-state data)
L0b  = L0 + qF0*F0;	% Nominal liquid flow below feed (kmol/min)
%lambda = 0;		    % Effect of vapor flow on liquid flow ("K2-effect")

% Inputs and disturbances
LT  = u(2*NT+3);                       % Reflux
VB  = u(2*NT+4);                       % Boilup
D   = u(2*NT+6);                       % Distillate
B   = u(2*NT+7);                       % Bottoms
F   = u(2*NT+5);                       % Feedrate
F_0 = dist.F_0;                     
zF  = dist.zF;                         % Feed composition         
qF  = dist.qF;                         % Feed liquid fraction



% THE MODEL
% objective function  
J = price.pf*F_0 + price.pV*VB - price.pB*B - price.pD*D;

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
dMdt{NT+1}  = SX.sym(['dMdt_' num2str(NT+1)],1);
dMxdt{NT+1} = SX.sym(['dMxdt_' num2str(NT+1)],1);

y = vertcat(y{:});
for i=1:NT-1
    y(i)=alpha*u(i)/(1+(alpha-1)*u(i));
end

% Vapor Flows assuming constant molar flows
V = vertcat(V{:});
for i=1:NT-1
    if i >= NF
        V(i) = VB + (1-qF)*F;
    else
        V(i) = VB;
    end
end

L = vertcat(L{:});
L(NT)=LT;
for i=2:NT-1
    if i <=NF
        L(i) = L0b + (u(NT+1+i)-Muw)./taul;
    else
        L(i) = L0  + (u(NT+1+i)-Muw)./taul;
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
dMxdt(NF)= dMxdt(NF) + F*u(NT+1);

% Reboiler (assumed to be an equilibrium stage)
dMdt(1) = L(2)      - V(1)      - B;
dMxdt(1)= L(2)*u(2) - V(1)*y(1) - B*u(1);

% Total condenser (no equilibrium stage)
dMdt(NT) = V(NT-1)         - LT       - D;
dMxdt(NT)= V(NT-1)*y(NT-1) - LT*u(NT) - D*u(NT);

% Compute the derivative for the mole fractions from d(Mx) = x dM + M dx
for i=1:(2*NT+2)
    ceq{i}  = SX.sym(['ceq_' num2str(i)],1);
end


% CSTR model
k1          = 34.1/60.0;
dMdt(NT+1)  = F_0 + D - F;
dMxdt(NT+1) = F_0*zF + D*u(NT) - F*u(NT+1) - k1*u(2*NT+2)*u(NT+1);

ceq  = vertcat(ceq{:});
for i=1:NT+1
    ceq(i) = dMxdt(i);
end

for i=1:NT+1
    ceq(NT+1+i) = dMdt(i);
end

% bound constraints
lb_u = [0.1; 0.1; 0.1; 0.1; 0.1];
ub_u = [10; 4.008; 10; 1.0; 1.0];

% State bounds and initial guess
x_min     = zeros(84,1);
x_max     = ones(84,1);
xB_max    = 0.1;
x_max(1)  = xB_max;
x_min(84) = 0.3;
x_max(84) = 0.7;
lbx  = [x_min;lb_u];
ubx  = [x_max;ub_u];
lbg  = zeros(2*NT+2,1);
ubg  = zeros(2*NT+2,1);

end

%function [H,Q,Qm] = Greshgorin(H,Qm)
function [H,Q] = Greshgorin(H)
numH    = size(H,1);
Q       = zeros(numH,numH);
%delta   = 2.5e-1;  % normal case
delta   = 2.5;      % with measurement noise 1 percent
for i=1:numH  % iterate all row of Hessian
    sumRow = 0;
    for j=1:numH
        if j ~= i
            sumRow = sumRow + abs(H(i,j));
        end
    end
    
    if H(i,i) <= sumRow   % include equality 
        Q(i,i) = sumRow - H(i,i) + delta;
    end
end
Q = diag(Q);
end
