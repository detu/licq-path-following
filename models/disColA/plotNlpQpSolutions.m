%PLOTNLPQPSOLUTIONS Summary of this function goes here
% 
% [OUTPUTARGS] = PLOTNLPQPSOLUTIONS(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2016/10/07 20:40:03 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2016

% load QP solution
%load qp_sol_10s.mat;   % pure-predictor
%load qp_sol_pg.mat;    % pure-predictor + linearized constraint
load  qp_sol_pc.mat;    % predictor-corrector

% load NLP solution
load nlp_sol.mat;

% load parameter setting
%load nlpQpDebug.mat;
% delta_t = 0.2;
% 
% numSteps   = 1/delta_t;

t_iter     = [0.01; 0.03; 0.09; 0.1; 0.15; 0.2; 0.4; 0.6; 0.8; 1.0]; 
numSteps   = size(t_iter,1);

err_primal = [];
err_dual   = [];
t_data     = [];
%t          = delta_t;
for i=1:numSteps
    err_p = norm((nlp_primal(:,i) - qp_primal(:,i)), 2);
    err_d = norm((nlp_dual(:,i) - qp_dual(:,i)), 2);
    err_primal = [err_primal err_p];
    err_dual   = [err_dual err_d];
    t          = t_iter(i,1);
    t_data     = [t_data t];
    %t          = t + delta_t;
end

plot(t_data, err_primal, 'r--o');