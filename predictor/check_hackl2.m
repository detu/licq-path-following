%CHECK_HACKL2 Summary of this function goes here
% 
% [OUTPUTARGS] = CHECK_HACKL2(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: detu $	$Date: 2015/11/04 23:43:52 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2015

x0 = [0,0];     % Make a starting guess at the solution
lb = [0, 0];
ub = [inf, inf];
p = 1;
A = [];
b = [];
Aeq = [];
beq = [];
options = optimset('LargeScale','off','Display', 'iter');
tic;
%[x, fval, ~,~, lamda] = fmincon(@(x)obj_hackl2(x),x0,[],[],[],[],lb,[],@(x)con_hackl2(x,p),options)
% Call Knitro to solve the optimization model.
% Specify some extra Knitro-specific options in "nlpoptions.opt".  Any Knitro
% options can be specified by passing in a Knitro options file. 
[x,fval,exitflag,output,lambda,grad,hess] = ...
   knitromatlab(@(x)obj_hackl2(x),x0,A,b,Aeq,beq,lb,ub,@(x)con_hackl2(x,p),[],options,'nlpoptions.opt');
toc;