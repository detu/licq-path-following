function J = computeNetCostDistCstr(Uopt,Uss,mpciterations)
%COMPUTENETCOSTDISTCSTR Summary of this function goes here
% 
% Net cost relative to steady state
%
% [OUTPUTARGS] = COMPUTENETCOSTDISTCSTR(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2016/09/09 17:49:28 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2016

pf   = 1; 
pV   = 0.01; 
pB   = 2; 
pD   = 0;
numU = 5; %[LT;VB;F;D;B]
Uopt = reshape(Uopt,mpciterations,numU);
Uss  = reshape(Uss,mpciterations,numU);
F_0  = 0.3;
J    = 0;
for i=1:mpciterations
    % compute cost function from optimization result
    VB = Uopt(i,2);
    D  = Uopt(i,4);
    B  = Uopt(i,5);
    Jo = pf*F_0 + pV*VB - pB*B - pD*D;
    
    % compute cost from steady-state solution
    VBs = Uss(i,2);
    Ds  = Uss(i,4);
    Bs  = Uss(i,5);
    Js  = pf*F_0 + pV*VBs - pB*Bs - pD*Ds;
    
    % compute net relative cost
    J  = J + Jo - Js;
end

end
