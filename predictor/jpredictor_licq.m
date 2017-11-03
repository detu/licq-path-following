function [x_init, y_init, elapsedqp] = jpredictor_licq(problem, p_init, p_final, x_init, y_init, delta_t, lb_init, ub_init, verbose_level, N)
%JPREDICTOR_LICQ Summary of this function goes here
% 
% jpredictor with LICQ, i.e., without LP solver.
%
% [OUTPUTARGS] = JPREDICTOR_LICQ(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2016/04/23 17:45:51 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2016
      

% assign problem to be solved 
clear prob;
sym pp;
p       = p_init;
theprob = @(pp)problem(pp);
prob    = theprob(p);
%t       = delta_t;
t       = 0;
alpha_1 = 0.25;
alpha_2 = 0.5;
% number of iteration
iter = 0;
%gamma   = 2.5;
gamma    = 1.5;
qpv_init = 0;
elapsedqp = 0;



% display problem
%fprintf('Solving problem %s with p=%3.1e using estimate at p=%3.1e\n',prob.name,p(1),p(2));
if (verbose_level)
    fprintf('Solving problem %s \n',prob.name);
    fprintf('iteration  delta_t        t        Success\n');
end

% calculate step s
%step = delta_t * (p_final - p_init);
nth = 1/delta_t;

while (t < 1)
    %while ( (t <= t_final) && (delta_t > 1e-6) )   % add stopping criteria delta_t not equal zero
    
%     % calculate step s
    step = delta_t * (p_final - p_init);
    
    % update bound constraint
    if(~isempty(lb_init))
        lb = lb_init - x_init;
        ub = ub_init - x_init;
    else
        lb = [];
        ub = [];
    end
    
    % solve QP problem
    % output from QP is y (directional derivative)
    [y, ~, qp_exit, lamda, qp_run] = solve_qp(prob, p_init, x_init, y_init, step, lb, ub, N);
    %[y, ~, qp_exit, lamda] = solve_qp(prob, p_init, x_init, y_init, step, [], [], N);
    elapsedqp = elapsedqp + qp_run;
    if (qp_exit < 0); % QP is infeasible
        
        % shorten step
        delta_t = alpha_2 * t;
        %t       = t - alpha_2 * delta_t;
        %t       = t - delta_t;
        
        % print out iteration and FAIL
        iter    = iter + 1;
        success = 0;
        if (verbose_level)
            fprintf('%6g     %3.1e     %3.1e   %4g \n',iter, delta_t , t, success);
        end
        
    else
        % QP is feasible
        
        % update states, multipliers, parameter, and time step
        x_init  = x_init + y;
        %y_init  = lamda.eqlin;
        y_init  = lamda;  % TOMLAB
        
        t       = t + delta_t;
        %p_init  = p_init + step + t*(p_final - p_init);
        p_init  = (1 - t)*p_init + t*p_final;
        
        % FIX STEPLENGTH
        fprintf('delta_t: %f\n', delta_t);
        fprintf('t: %f\n', t);
        
        % update delta_t
        nth = nth - 1;
        delta_t = 1/nth;
        
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
function [y, qp_val, qp_exit, lamda, elapsedqp] = solve_qp(prob, p, x_init, y_init, step, lb, ub, N)
    
    % setup objective function for QP problem
    %[~,g,H,Lxp,cin,J,cp,Jeq,dpe] = prob.obj(x_init,y_init,p);
    %[~,g,H,Lxp,cin,~,~,Jeq,dpe,HObj] = prob.obj(x_init,y_init,p, N);
    [~,g,H,Lxp,cin,~,~,Jeq,dpe,~] = prob.obj(x_init,y_init,p, N);
    
    % Gershgorin Convexification
    %rho = 0.1;
    %Q   = gershgorinConvex(H,rho);
    %Q   = gershgorinConvex(HObj,rho);
    %Q   = computeWeightsGreshgorin(x_init,N);
    A   = [];
    b   = [];
    f   = Lxp * step + g;
    
    % Only Equality Constraint
    ceq = cin;
    Aeq  = Jeq;
    beq  = -dpe*step - ceq;  % just for case 1 parameter! need to think further!!!
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Finally solve QP problem
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     %options = optimset('Display','off','Algorithm','active-set');
%     options = optimset('Display','off','Algorithm','interior-point-convex');
%     %options = optimset('Display','off','Algorithm','interior-point-convex', 'TolFun', 1e-6, 'TolX', 1e-6);
% %     options = optimset('Display','off','Algorithm','interior-point-convex', 'MaxIter', 4);
% %     %options = optimset('Display','on', 'MaxIter', 2);
% %     %[y, qp_val, qp_exit, ~, lamda] = quadprog(H,f,A,b,Aeq,beq,lb,ub,[],options);
%     startqp = tic;
%     [y, qp_val, qp_exit, ~, lamda] = quadprog(H,f,A,b,Aeq,beq,lb,ub,[],options);
%     %[y, qp_val, qp_exit, ~, lamda] = quadprog(H,f,A,b,Aeq,beq,lb,ub,x_init,options);  % supply initial solution
%     %[y, qp_val, qp_exit, ~, lamda] = quadprog(H,f,A,b,Aeq,beq,lb,ub,x_init);
%     %[y, qp_val, qp_exit, ~, lamda] = quadprog(H+Q,f,A,b,Aeq,beq,lb,ub,x_init,options);
%     %[y, qp_val, qp_exit, ~, lamda] = quadprog(Q,f,A,b,Aeq,beq,lb,ub,x_init,options);
%     elapsedqp = toc(startqp);
%     fprintf('QP solver runtime: %f\n',elapsedqp);
%     fprintf('QP return: %d\n', qp_exit);
   
%    % test using qpOASES
%    %nC = size(Aeq,2);   
%    aC  = ones(3280,1);
%    auxInput = qpOASES_auxInput('hessianType', 5, 'x0', x_init, 'guessedWorkingSetC', aC);
%    %options  = qpOASES_options( 'printLevel',2, 'maxIter', 100 , 'initialStatusBounds', 0);
%    %auxInput = qpOASES_auxInput('hessianType', 5, 'x0', x_init);
%    options  = qpOASES_options( 'printLevel',2, 'initialStatusBounds', 0);
%    startqp  = tic;
%    %[y,qp_val,qp_exit,iter,lamda,auxOutput] = qpOASES(H,f,Aeq,lb,ub,beq,beq,options,auxInput);
%    %[y,qp_val,qp_exit,iter,lamda,auxOutput] = qpOASES( H,f,Aeq,lb,ub,beq,beq,options);
%    [y,qp_val,qp_exit,iter,lamda,auxOutput] = qpOASES(H,f,Aeq,[],[],beq,beq,options,auxInput);
%    elapsedqp = toc(startqp);
%    fprintf('QP solver runtime: %f\n',elapsedqp);
%    fprintf('QP return: %d\n', qp_exit);

    % test using TOMLAB
    %Prob   = qpAssign(H, f, Aeq, beq, beq, lb, ub, x_init);
    %x0_guess = 0*x_init;
    %Prob   = qpAssign(H, f, Aeq, beq, beq, lb, ub, x0_guess);
    Prob   = qpAssign(H, f, Aeq, beq, beq, lb, ub);
    %Prob   = qpAssign(H+Q, f, Aeq, beq, beq, lb, ub, x_init);
%     Prob.optParam.eps_f = 1e-4;
%     Prob.optParam.eps_x = 1e-4;
%     Prob.optParam.fTol  = 1e-4;
%     Prob.optParam.xTol  = 1e-4;
    startqp  = tic;
    Result = tomRun('qp-minos', Prob, 1);
    elapsedqp = toc(startqp);
    fprintf('QP solver runtime: %f\n',elapsedqp);
    y       = Result.x_k;
    qp_val  = Result.f_k;
    qp_exit = Result.ExitFlag;

    %lamda   = Result.v_k(size(y_init));
    numX    = size(x_init,1);
    %numY    = size(y_init,1);
    lamda   = Result.v_k(numX+1:end);
    fprintf('QP return: %d\n', qp_exit);
end

function xinput = checkBoundAdjust(xinput,lb,ub)
    % check xmeasure with bound constraints
    numX = numel(xinput);
    for i=1:numX
        if xinput(i) < lb(i)
            xinput(i) = lb(i);
        end
        
        if xinput(i) > ub(i)
            xinput(i) = ub(i);
        end
    end
end
