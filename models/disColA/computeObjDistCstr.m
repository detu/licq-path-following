function J = computeObjDistCstr(Uopt, mpciterations)
%COMPUTEOBJDISTCSTR Summary of this function goes here
% 
% compute objective function for NMPC simulation
%
% [OUTPUTARGS] = COMPUTEOBJDISTCSTR(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2016/08/09 15:14:06 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2016

pf   = 1; 
pV   = 0.01; 
pB   = 2; 
pD   = 0;
numU = 5; %[LT;VB;F;D;B]
Uopt = reshape(Uopt,mpciterations,numU);
F_0  = 0.3;
J    = 0;
for i=1:mpciterations
    VB = Uopt(i,2);
    D  = Uopt(i,4);
    B  = Uopt(i,5);
    J  = J + pf*F_0 + pV*VB - pB*B - pD*D;
end

end
