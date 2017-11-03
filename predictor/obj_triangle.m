function f = obj_triangle(x,p)
%OBJ_TRIANGLE Summary of this function goes here
% 
% [OUTPUTARGS] = OBJ_TRIANGLE(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2016/01/05 13:19:54 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2016

f = -p*x(1)*x(2)*x(4) + (1-p)*( ( x(1) - 1.5 )^2 + x(2)^2 + ( x(3) - 1.5 )^2 + ( x(4) - 1 )^2 + x(5)^2 ); 
end
