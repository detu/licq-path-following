function Jobj = computeObjectiveFunctionValues(uOpt,xActual)
%COMPUTEOBJECTIVEFUNCTIONVALUES Summary of this function goes here
% 
% [OUTPUTARGS] = COMPUTEOBJECTIVEFUNCTIONVALUES(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2017/01/04 19:26:57 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2017

% prices
pf  = 1;
pV  = 0.02;
pB  = 2;
pD  = 0;
F_0 = 0.3;

% steady-state values
load CstrDistXinit.mat;
xs    = Xinit(1:84);
us    = Xinit(85:89);
nx    = size(xs,1);
nu    = size(us,1);
% weights
%load Qmax.mat;
load Q.mat;
Qmax = Q;
c1  = -0.05; % noisy case
lss = -0.256905910000000 + c1; %steady-state objective function value

%J = (pf*F_0 + pV*uOpt(2) - pB*uOpt(5) - pD*uOpt(4)) + (Qmax(nx+1:nx+nu,1).*(uOpt - us))' * (uOpt - us) + (Qmax(1:nx,1).*(xActual - xs))' * (xActual - xs);
Jecon    = (pf*F_0 + pV*uOpt(2) - pB*uOpt(5) - pD*uOpt(4));
Jcontrol = (Qmax(nx+1:nx+nu,1).*(uOpt - us))' * (uOpt - us);
Jstate   = (Qmax(1:nx,1).*(xActual - xs))' * (xActual - xs);
J        = Jecon + Jcontrol + Jstate - lss;
fprintf('-------------------------------------------------------------------- \n');
fprintf('Jecon: %f \t, Jcontrol: %f, \t Jstate: %f \n', Jecon,Jcontrol,Jstate);
Jobj.reg    = J;
Jobj.econ   = Jecon;
end
