%CHECK_TRIANGLE Summary of this function goes here
% 
% [OUTPUTARGS] = CHECK_TRIANGLE(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2016/01/05 13:18:56 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2016

format long;
x0 = [1.5, 0.0, 1.5, 1.0, 0.0];     % Make a starting guess at the solution
p = 0;
options = optimset('LargeScale','off','Display', 'iter');
tic;
[x, fval, ~,~, lamda] = fmincon(@(x)obj_triangle(x,p),x0,[],[],[],[],[],[],@(x)con_triangle(x,p),options)
%[x,fval,exitflag,output,lambda,grad,hessian] = knitromatlab(@(x)obj_triangle(x,p),x0,[],[],[],[],[],[],@(x)con_triangle(x,p))
toc;
