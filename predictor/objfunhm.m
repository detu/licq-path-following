function f = objfunhm(x,p1)
%OBJFUNHM Summary of this function goes here
% 
% [OUTPUTARGS] = OBJFUNHM(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2015/10/22 18:19:06 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2015

%p1 = 8;
f = p1*x(1)^3 + x(2)^2;

end
