function [c, ceq] = con_hackl2(x,p)
%CONFUNHM Summary of this function goes here
% 
% [OUTPUTARGS] = CONFUNHM(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2015/10/22 18:19:48 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2015

%p2 = 1;
% Nonlinear inequality constraints
%c = [-x(1); -x(2); -x(1) + 2*x(2) - 2 + 4*p; 2*x(1) + x(2) - 1 - 5*p];
c = [-x(1) + 2*x(2) - 2 + 4*p; 2*x(1) + x(2) - 1 - 5*p];
% Nonlinear equality constraints
ceq = [];

end
