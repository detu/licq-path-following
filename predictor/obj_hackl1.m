function f = obj_hackl1(x,p)
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

f = x(1)^2 + (4*p-2)*x(2) + 0.5*x(2)^2;
end