function [c, ceq] = confunhm(x,p2)
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
c = [exp(-x(1)) - x(2);     
     p2 - x(1)];
% Nonlinear equality constraints
ceq = [];

% test for equality constraint
% Nonlinear inequality constraints
% ceq = [x(2) - exp(-x(1));     
%        x(1) - p2 ];
% ceq = [exp(-x(1)) - x(2);     
%        p2 - x(1) ];

% Nonlinear equality constraints
% c = [];

end
