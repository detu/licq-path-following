function [x_init, y_init, elapsedqp] = jpredictor_licq_pure(problem, p_init, p_final, x_init, y_init, delta_t, lb_init, ub_init, verbose_level, N)
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
numX    = size(x_init,1);
%x0      = 1e-5*ones(numX,1);
x0      = zeros(numX,1);   % TRY THIS

% % store result for debugging, do compasion with IPOPT solutions
% qp_primal = [];
% qp_dual   = [];

% display problem
%fprintf('Solving problem %s with p=%3.1e using estimate at p=%3.1e\n',prob.name,p(1),p(2));
if (verbose_level)
    fprintf('Solving problem %s \n',prob.name);
    fprintf('iteration  delta_t        t        Success\n');
end

% t_iter     = [0.01; 0.03; 0.09; 0.1; 0.15; 0.2; 0.4; 0.6; 0.8; 1.0]; 
% numStep    = size(t_iter,1);

% % update bound constraint
% if(~isempty(lb_init))
%     lb = lb_init - x_init;
%     ub = ub_init - x_init;
% else
%     lb = [];
%     ub = [];
% end

p_0 = p_init;

while (t <= 1)
%for i=1:numStep
    %while ( (t <= t_final) && (delta_t > 1e-6) )   % add stopping criteria delta_t not equal zero
    
    % calculate step s
    %step = delta_t * (p_final - p_init);
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
    
    %step = containStepWithinBounds(lb,ub,step);
    
    % solve QP problem
    % output from QP is y (directional derivative)
    %[y, ~, qp_exit, lamda, qp_run] = solve_qp(prob, p_init, x_init, y_init, step, lb, ub, N);
    [y, ~, qp_exit, lamda, qp_run] = solve_qp(prob, p_init, x_init, y_init, step, lb, ub, N, x0, lb_init, ub_init);  % supply initial guess
    elapsedqp = elapsedqp + qp_run;
    if (qp_exit < 0) % QP is infeasible
        
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
        
        % check x_init if outside bound constraint ? 
        
        %y_init  = lamda.eqlin;
        %y_init  = y_init + lamda;
        %y_init  = lamda;  % TOMLAB
        y_init.lam_x = lamda;
        
        % for debugging
        %qp_primal = [qp_primal x_init];
        %qp_dual   = [qp_dual y_init];
        
        t       = t + delta_t;
        p_init  = p_t;
        
        % update initial guess
        %x0   = x_init + 1e-5*ones(numX,1);
        %x0   = x_init;
        %x0    = y;  % try not to update initial guess
        
        
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

%save qp_sol_10s.mat qp_primal qp_dual;
%save qp_sol_pg.mat qp_primal qp_dual;
%save qp_sol_pc.mat qp_primal qp_dual;

%save qp_sol_100s.mat qp_primal qp_dual;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOLVE QP program
% solution of QP program: [y] (directional derivative)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y, qp_val, qp_exit, lamda, elapsedqp] = solve_qp(prob, p, x_init, y_init, step, lb, ub, N, x0, lb_init, ub_init)
    
    % setup objective function for QP problem
    %[~,g,H,Lxp,cin,J,cp,Jeq,dpe] = prob.obj(x_init,y_init,p);
    %[~,g,H,Lxp,cin,~,~,Jeq,dpe,HObj] = prob.obj(x_init,y_init,p, N);
    
%     % update x_init 
%     nS = size(step,1);
%     x_init(1:nS) = x_init(1:nS) + step;
    
    [~,g,H,Lxp,cin,~,~,Jeq,dpe,~] = prob.obj(x_init,y_init,p, N);
    
    % FOR DEBUGGING KEEP THE DERIVATIVES FIX !
    % save derivatives.mat g H Lxp cin Jeq dpe
    %load derivatives.mat;
    
    % Gershgorin Convexification
    %rho = 0.1;
    %Q   = gershgorinConvex(H,rho);
    %Q   = gershgorinConvex(HObj,rho);
    %Q   = computeWeightsGreshgorin(x_init,N);
    A   = [];
    b   = [];
    %f   = Lxp * step + g;
    f   = Lxp * step;
    %f   = Lxp * step + Jeq'*y_init;
    
    % Only Equality Constraint
    ceq = cin;
    Aeq  = Jeq;
    %beq  = -dpe*step - ceq;  % type 2
    beq  = -dpe*step;       % for soft-constraint 
    
    % CHECK LAGRANGE MULTIPLIERS FROM BOUND CONSTRAINTS 
    %lmC = y_init.lam_x;
    lmC = abs(y_init.lam_x);
%     bAc = find(lmC >= 1e-6);
%     bWc = find(lmC < 1e-6);
    %bAc = find(lmC >= 1e-6);
    %bAc = find(lmC >= 1e-7);
    bAc = [];
    bWc = 1:length(lmC);
    
    % build equality constraint from active bound constraints
    numBaC = size(bAc,1);
    numPri = size(x_init,1);
    AeqB   = zeros(numBaC,numPri);
    beqB   = zeros(numBaC,1);
    

    
    for i=1:numBaC
        % populate matrix AeqB and beqB
        indB         = bAc(i);
        AeqB(i,indB) = 1;
        
        % need to check if the Lagrange multiplier is positive or negative?
        lamBC        = y_init.lam_x(indB);
        if lamBC > 0
            beqB(i,1) = ub(indB);    % active upper-bound
        else
            beqB(i,1) = lb(indB);    % active lower-bound
        end
        
% %         % remove bound constraints from the active constraint list
%         lb(indB)     = -inf;
%         ub(indB)     = inf;
        %lb(indB)     = [];
        %ub(indB)     = [];
        %beq(indB,1)  = 0;
        %Aeq(indB,:)  = 0;
        %Aeq(indB,indB) = 1;
        %ub(indB)    = [];
        %lb(indB)    = [];
        %beq(indB)   = [];
        %Aeq(indB,:) = [];
    end
    
    % weakly active constraint
    numBwC = size(bWc,1);
    AeqBwU  = zeros(numBwC,numPri);
    AeqBwL  = zeros(numBwC,numPri);
    beqBwU  = zeros(numBwC,1);
    beqBwL  = zeros(numBwC,1);
    for i=1:numBwC
        indW         = bWc(i);
        AeqBwU(i,indW) = 1;
        %AeqBwL(i,indW) = -1;
        
        % need to check if the Lagrange multiplier is positive or negative?
        lamBC        = y_init.lam_x(indW);
        if lamBC > 0
            beqBwU(i,1) = ub_init(indW) - x_init(indW);
            %beqBwL(i,1) = ub_init(indW) - x_init(indW);
        else
            beqBwL(i,1) = lb_init(indW) - x_init(indW);
            %beqBwL(i,1) = ub_init(indW) - x_init(indW);
        end
    end
    
    % Modify matrix Aeq and beq
    %Aeq = [Aeq; AeqB];
    %beq = [beq; beqB];
    
    Aeq = [Aeq; AeqB; AeqBwU];
    beqU = [beq; beqB; beqBwU];
    beqL = [beq; beqB; beqBwL];
    
    %beqU = [beq; inf(size(beqB)); inf*beqBwU];
    %beqL = [beq; -inf(size(beqB)); -inf*beqBwL];
    
    % try to remove bound constraint here
    lb = [];
    ub = [];
%     lb = -100*ones(size(x_init));
%     ub = 100*ones(size(x_init));
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Finally solve QP problem
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %options = optimset('Display','off','Algorithm','active-set');
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
   
%    % test USING QPOPT
%     global otxProb
%     otxProb = ProbDef;
%     otxProb.SolverQP = 'qpopt';
%     startqp = tic;
%     [y, qp_val, qp_exit, ~, lamda] = quadprog(H,f,A,b,Aeq,beq,lb,ub);
%     elapsedqp = toc(startqp);
%     fprintf('QP solver runtime: %f\n',elapsedqp);
%     fprintf('QP return: %d\n', qp_exit);

%    % test QPAS
%     dsp=1;
%     startqp = tic;
%     %[y, qp_val, qp_exit, ~, lamda] = qpas(H,f,[],[],Aeq,beq,lb,ub,dsp);
%     %[y, qp_exit, lamda] = qpas(full(H),f,[],[],full(Aeq),full(beq),lb,ub,dsp);
%     [y, qp_exit, lamda] = qpip(full(H),f,[],[],full(Aeq),full(beq),lb,ub,dsp);
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
    %Prob   = qpAssign(H, f, Aeq, beq, beq, lb, ub, y);
    
    %Prob   = qpAssign(H, f, Aeq, beq, beq, lb, ub);
    %Prob   = qpAssign(H, f, Aeq, beq, beq);
    
    %Prob   = qpAssign(H, f, Aeq, beq, beq, lb, ub, x0);   %supply initial guess
    %Prob   = qpAssign(H, f, Aeq, beq, beq, lb, ub, x0,'WarmStart','true');
    %Prob   = qpAssign(H, f, Aeq, beqL, beqU, lb, ub, x0,'WarmStart','true');
    Aeq  = sparse(Aeq);
    beqL = sparse(beqL);
    beqU = sparse(beqU);
    f = sparse(f);
    x0 = sparse(x0);
    Prob   = qpAssign(H, f, Aeq, beqL, beqU, lb, ub, x0,'PriLevOpt', 1);
    Prob.optParam.eps_x = 1e-7;
    Prob.optParam.fTol  = 1e-7;
    Prob.optParam.xTol  = 1e-7;
    Prob.PriLevOpt = 5;
    Prob.PriLev = 1;
    startqp  = tic;
    Result = tomRun('qp-minos', Prob, 1);
    %Result = tomRun('qpopt', Prob, 1);
    elapsedqp = toc(startqp);
    fprintf('QP solver runtime: %f\n',elapsedqp);
    qp_exit = Result.ExitFlag;
    if qp_exit == 0
        y       = Result.x_k;
        qp_val  = Result.f_k;
    else
        keyboard;
        numPri = size(x_init,1);
        y = zeros(numPri,1);
        qp_val = Result.f_k;
    end
    

    %lamda   = Result.v_k(size(y_init));
    numX    = size(x_init,1);
    %numY    = size(y_init,1);
    lamda   = Result.v_k(numX+1:end);
    fprintf('QP return: %d\n', qp_exit);
    
%     % set solution to ZERO for active bound constraints
%     y(bAc(:)) = 0;
end

function step = containStepWithinBounds(lb,ub,step)

stepM    = -step;  % because step becomes beq in equality constraint
numStep  = size(step,1);

for i=1:numStep
    if stepM(i,1) > ub(i,1)        % if step greater than bound upper 
        step(i,1) = ub(i,1); 
    elseif stepM(i,1) < lb(i,1)  % if step less than lower bound
        step(i,1) = lb(i,1);
    else
        % do nothing
    end
end
end


