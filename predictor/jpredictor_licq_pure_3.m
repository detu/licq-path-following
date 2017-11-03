function [x_init, y_init, elapsedqp] = jpredictor_licq_pure_3(problem, p_init, p_final, x_init, y_init, delta_t, lb_init, ub_init, verbose_level, N)
%JPREDICTOR_LICQ_PURE Summary of this function goes here
% 
% [OUTPUTARGS] = JPREDICTOR_LICQ_PURE(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2016/08/13 01:22:15 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2016

clear prob;
sym pp;
p       = p_init;
theprob = @(pp)problem(pp);
prob    = theprob(p);
t       = 0;
alpha_2 = 0.5;
% number of iteration
iter = 0;
elapsedqp = 0;
numX    = size(x_init,1);
x0      = zeros(numX,1);   
if (verbose_level)
    fprintf('Solving problem %s \n',prob.name);
    fprintf('iteration  delta_t        t        Success\n');
end

p_0 = p_init;

while (t <= 1)
    
    % calculate step s
    tk   = t + delta_t;
    p_t  = (1 - tk)*p_0 + tk*p_final;
    step = p_t - p_init;
    
    % update bound constraint
    if(~isempty(lb_init))
        lb = lb_init - x_init; % (0 - 0.0999)
        ub = ub_init - x_init; %(0.1 - 0.0999)
    else
        lb = [];
        ub = [];
    end

    
    % solve QP problem
    % output from QP is y (directional derivative)
    [y, ~, qp_exit, lamda, qp_run] = solve_qp(prob, p_init, x_init, y_init, step, lb, ub, N, x0, lb_init, ub_init);  % supply initial guess
    elapsedqp = elapsedqp + qp_run;
    if (qp_exit < 0) % QP is infeasible
        
        % shorten step
        delta_t = alpha_2 * t;
        
        % print out iteration and FAIL
        iter    = iter + 1;
        success = 0;
        if (verbose_level)
            fprintf('%6g     %3.1e     %3.1e   %4g \n',iter, delta_t , t, success);
        end
        
    else
        % QP is feasible
        
        % update states, multipliers, parameter, and time step
        x_init       = x_init + y; 
        y_init.lam_x = y_init.lam_x - lamda.lam_x; % - (minus) because TOMLAB has different sign compared to IPOPT

        t       = t + delta_t;
        p_init  = p_t;
        
        
        % FIX STEPLENGTH
        fprintf('delta_t: %f\n', delta_t);
        fprintf('t: %f\n', t);
        
        
        % print out iteration and SUCCESS
        fprintf('--------------------------------------------------- \n');
        iter    = iter + 1;
        fprintf('iteration number: %d\n', iter);
        success = 1;
        if (verbose_level)
            fprintf('%6g     %3.1e     %3.1e   %4g \n',iter, delta_t , t, success);
        end
        
    end
    if ( (1-t) <= 1e-5 )
        break;
    end
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOLVE QP program
% solution of QP program: [y] (directional derivative)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y, qp_val, qp_exit, lamda, elapsedqp] = solve_qp(prob, p, x_init, y_init, step, lb, ub, N, x0, lb_init, ub_init)
    

    % obtain derivatives information
    [~,g,H,Lxp,cin,~,~,Jeq,dpe,~] = prob.obj(x_init,y_init,p, N);    

    % QP setup
    A   = [];
    b   = [];
    f   = Lxp * step + g;
    
    % Only Equality Constraint
    ceq  = cin;
    Aeq  = Jeq;
    beq  = dpe*step + ceq;   %OK

       
    % CHECK LAGRANGE MULTIPLIERS FROM BOUND CONSTRAINTS 
    lmC = abs(y_init.lam_x);
    bAc = find(lmC >= 1e-3);
    
    % build equality constraint from active bound constraints
    numBaC = size(bAc,1);
    for i=1:numBaC        
        % put strongly active constraint on boundary
        indB         = bAc(i);
        
        ub(indB)     = 0;        % keep upper bound on boundary
        lb(indB)     = 0; 
        
    end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Finally solve QP problem
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    Prob   = qpAssign(H, f, Aeq, beq, beq, lb, ub);
    Prob.optParam.eps_x = 1e-7;
    Prob.optParam.fTol  = 1e-7;
    Prob.optParam.xTol  = 1e-7;
    %Prob.PriLevOpt = 5;
    %Prob.PriLev = 1;
    startqp  = tic;
    Result = tomRun('qp-minos', Prob, 1);

    elapsedqp = toc(startqp);
    fprintf('QP solver runtime: %f\n',elapsedqp);
    qp_exit = Result.ExitFlag;
    if qp_exit == 0
        y       = Result.x_k;
        qp_val  = Result.f_k;
    else
        keyboard;
    end

    numX        = size(x_init,1);
    lamda.lam_x = Result.v_k(1:numX);
    lamda.lam_g = -Result.v_k(numX+1:end);
    fprintf('QP return: %d\n', qp_exit);
    
end




