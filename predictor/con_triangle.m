function [c, ceq] = con_triangle(x,p)
%CON_TRIANGLE Summary of this function goes here
% 
% [OUTPUTARGS] = CON_TRIANGLE(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2016/01/05 13:19:26 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2016

% Nonlinear inequality constraints
c = [];
% Nonlinear equality constraints
ceq = [ x(1) + x(2) + x(3) - 3; ...
        x(1)^2 + x(2)^2 - x(3)^2 - 2*x(1)*x(2)*x(5);...
        x(4)^2 + x(5)^2 - 1];

end
