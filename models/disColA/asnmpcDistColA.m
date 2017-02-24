function asnmpcDistColA
%ASNMPCDISTCOLA Summary of this function goes here
% 
% Implementation of Advanced-step NMPC for Distillation Column A
% 
% [OUTPUTARGS] = ASNMPCDISTCOLA(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2016/03/30 15:01:17 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2016

import casadi.* 
format long;

% advanced-step NMPC algorithm, as described in Huang's PhD thesis (Chapter 3 page 33),
% consists of:
% 1. Background calculation: Having xk and uk at time step k, predict the future state
% of the system zk+1 using the dynamic model. Assuming the computations can be
% completed within one sampling time, solve the NLP (3.9) based on (3.14) with p0 =
% zk+1.
% 2. On-line update: At time step k+1, obtain the true state xk+1. Set p =
% xk+1 and use (3.12) to get the fast updated solution. Extract the control action uk+1 from the
% approximate solution vector and inject to the plant.
% 3. Iterate: Set k k+1 and go to background.

% TO DO
% 1. Parameters should be initial state values ! (CHANGE !)


% QUESTION: how to get mu (last barrier multiplier) from IPOPT ?
% answer: - by reading IPOPT iteration output ?

%===============================================================================
% Initiate guess both for primal and dual variables
%===============================================================================
load nlp0.mat;
xstart = x_opt;
ystart = y_opt;
bstart = b_opt;

% need to iterate in bstart to get VL and VU
VL = zeros(numel(bstart),1);
VU = zeros(numel(bstart),1);
for i=1:numel(bstart)
    if bstart(i) < 0
        VL(i) = bstart(i);
    else
        VU(i) = bstart(i);
    end
end
VL = -VL;  % CHECK AGAIN IF THIS IS NECESSARY
  

% P is initial states 
load cola_init.mat;
p_init = Xinit;
load u_opt_ss.mat;
p_final = xf(1:82);

lb_init = lb;
ub_init = ub;

[primal, dual] = sIPOPT(@(p)distColA_casadi_sens(p), p_init, p_final, xstart, ystart, VL, VU, lb_init, ub_init);


end
