%CHECK_HACK1 Summary of this function goes here
% 
% [OUTPUTARGS] = CHECK_HACK1(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: detu $	$Date: 2015/11/04 22:47:06 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2015

x0 = [0,1];     % Make a starting guess at the solution
p = 1;
options = optimset('LargeScale','off','Display', 'iter');
tic;
[x, fval, ~,~, lamda] = fmincon(@(x)obj_hackl1(x,p),x0,[],[],[],[],[],[],@(x)con_hackl1(x,p),options)
toc;