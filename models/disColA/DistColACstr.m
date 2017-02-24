function [t,states,xdot,inputs] = DistColACstr(U)
%DISTCOLACSTR Summary of this function goes here
% 
% Distillation column A with CSTR model included.
%
% [OUTPUTARGS] = DISTCOLACSTR(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2016/03/14 06:29:36 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2016

import casadi.* 

% model parameters
NC    = 2; 
NF    = 21;
NT    = 41; 
q     = 1.0; 
alpha = 1.5;
z_F0  = 1.0;
Muw   = 0.5;
F_0   = U;

% states and control inputs
x      = SX.sym('x',NT+1,NC-1);   % include CSTR model
M      = SX.sym('M',NT+1,1);      % include CSTR model
states = [x;M];
L_T    = SX.sym('L_T');
V_B    = SX.sym('V_B');
F      = SX.sym('F');
D      = SX.sym('D');
B      = SX.sym('B');
inputs = [L_T;V_B;F;D;B];


t  = SX.sym('t');
y  = SX.sym('y', NT-1, NC-1);             % Vapor composition
Li = SX.sym('Li', NT, 1);                 % Liquid flow
Vi = SX.sym('Vi', NT, 1);                 % Vapor flow

dMdt  = SX.sym('dMdt', NT+1, 1);            % Total molar holdup
dMxdt = SX.sym('dMxdt', NT+1, NC-1);        % Componentwise molar holdup
dxdt  = SX.sym('dxdt', NT+1, NC-1);         % Rate of change of composition

% Vapor flows (Constant, no dynamics) V[NT] skipped since D given
for i=1:NT-1
    Vi(i) = V_B;
    if i >= NF
        Vi(i) = Vi(i) + (1-q)*F;
    end
end


Vi(NT) = inf;

% Liquid flows (Wier formula) L[1] not required since B given
taul  = 0.063;        % time constant for liquid dynamics (min)
F0    = 1;            % Nominal feed rate (kmol/min) (NEED TO CHECK THE VALUE !!!)
qF0   = 1;            % Nominal fraction of liquid in feed  (NEED TO CHECK THE VALUE !!!)
L0    = 2.70629;      % Nominal reflux flow (from steady-state data) (NEED TO CHECK THE VALUE !!!)
L0b   = L0 + qF0*F0;  % Nominal liquid flow below feed (kmol/min)
Li(1) = inf;

for i=2:NT
    if i<=NF
        Li(i) = L0b + (M(i) - Muw)/taul;
    else
        Li(i) = L0 + (M(i) - Muw)/taul;
    end
end
% Top tray liquid
Li(NT) = L_T;

% Vapor liquid equilibrium
for i=1:NT-1
    for j=1:NC-1
        y(i,j) = (x(i,j)*alpha) / ( 1 + ((alpha-1)*x(i,j)) );
    end
end


% Partial Reboiler (for simplicity one can use L(0) = B)
dMdt(1) = Li(2) - Vi(1) - B;
for j=1:NC-1
    dMxdt(1,j)  = Li(2)*x(2,j) - Vi(1)*y(1,j) - B*x(1,j);
end

% Stripping and  Enrichment sections
for i=2:NT-1
    dMdt(i)  = Li(i+1) - Li(i)+ Vi(i-1) - Vi(i);
    for j=1:NC-1
        dMxdt(i,j) = Li(i+1)*x(i+1,j) - Li(i)*x(i,j) + Vi(i-1)*y(i-1,j) - Vi(i)*y(i,j);
    end
end

% Correction for the feed stage:
dMdt(NF) = dMdt(NF) + F;
for j=1:NC-1
    dMxdt(NF,j) =  dMxdt(NF,j) + F*x(NT+1);
end

% Total Condenser (for simplicity one can use V(NT-1) = D )
dMdt(NT) = Vi(NT-1) - Li(NT) - D;
for j=1:NC-1
    dMxdt(NT,j) = Vi(NT-1)*y(NT-1,j) - Li(NT)*x(NT,j)  - D*x(NT,j);
end

% CSTR model
k1          = 34.1/60.0;
dMdt(NT+1)  = F_0 + D - F;
for j=1:NC-1
    dMxdt(NT+1,j) = F_0*z_F0(j,1) + D*x(NT,j) - F*x(NT+1,j) - k1*M(NT+1)*x(NT+1,j);
end

for i=1:NT+1
    for j=1:NC-1
        dxdt(i,j) = (dMxdt(i,j)-x(i,j)*dMdt(i))/M(i);
    end
end

xdot = [dxdt;dMdt];
end
