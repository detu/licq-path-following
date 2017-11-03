%CHECK_PROB4 Summary of this function goes here
% 
% [OUTPUTARGS] = CHECK_PROB4(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2015/10/22 18:08:41 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2015

%x0 = [0.5,0.6];     % Make a starting guess at the solution
x0 = [0.496811686457911, 0.608467572142176];
%x0 = [1.000000000000000, 0.367879441171442];
%x0 = [-4.000000000000000; 54.598150033144236];
% p1 = 1;
% p2 = -4;
% p1 = 8;
% p2 = 1;
p1 = 1.7;
p2 = -3.5;
options = optimset('LargeScale','off','Display', 'iter');
tic;
[x, fval, exit,~, lamda] = fmincon(@(x)objfunhm(x,p1),x0,[],[],[],[],[],[],@(x)confunhm(x,p2),options);
toc;