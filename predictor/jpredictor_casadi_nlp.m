function [x_init, y_init, info, nlpResult] = jpredictor_casadi_nlp(problem, p_init, p_final, x_init, y_init, delta_t, lb_init, ub_init, verbose_level, nlpRun)
%JPREDICTOR_CASADI_NLP Summary of this function goes here
% 
% [OUTPUTARGS] = JPREDICTOR_CASADI_NLP(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2016/02/02 19:58:29 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2016

% initial value
in_info   = 1;      % index for info
info(in_info).t = 0;
info(in_info).x = x_init;
info(in_info).y = y_init;       

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

% initial dual variable
delta_l  = y_init;

% Comparison NLP and PF
t_data = [];
x_data = [];
y_data = []; 
nlpResult = [];


% display problem
%fprintf('Solving problem %s with p=%3.1e using estimate at p=%3.1e\n',prob.name,p(1),p(2));
if (verbose_level)
    fprintf('Solving problem %s \n',prob.name);
    fprintf('iteration  delta_t        t        Success\n');
end

while (t < 1)
%while ( (t <= t_final) && (delta_t > 1e-6) )   % add stopping criteria delta_t not equal zero
    
    % calculate step s
    step = delta_t * (p_final - p_init);
    %step = (p_final - p_init);
    
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
    [y, qp_val, qp_exit, oqp, k_zero_tilde, k_plus_tilde, g, lamda] = solve_qp(prob, p_init, x_init, y_init, step, lb, ub);
    
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
        
    else % QP is feasible
        % solve LP
        [delta_l, lp_exit] = solve_lp(prob, oqp, y_init, step, y, k_zero_tilde, k_plus_tilde, g, lamda);
        
        % check exit flag
        if (lp_exit < 0) % LP is infeasible
            
            % shorten step
            %delta_t = alpha_1 * t;
            %t       = t - alpha_1 * delta_t;
            %t       = t - delta_t;
            delta_t = alpha_1 * delta_t;
            
            % print out iteration and FAIL
            iter    = iter + 1;
            success = 0;
            if (verbose_level)
                fprintf('%6g     %3.1e     %3.1e   %4g \n',iter, delta_t , t, success);
            end
            
        else
            % trust-region update
%             [Jo,~] = prob.obj(x_init,p_init);
            x_old  = x_init;
            % update states, multipliers, parameter, and time step
            x_init  = x_init + y;
            %y_init  = y_init + delta_l;  % TRY FOR EQUALITY CONSTRAINT - LARGER ERROR!
            y_init  = delta_l;
            
            
            t       = t + delta_t;
            %p_init  = p_init + step + t*(p_final - p_init);
            p_init  = (1 - t)*p_init + t*p_final;
            
            % FIX STEPLENGTH
            %delta_t = min(delta_t, (1-t));
            fprintf('delta_t: %f\n', delta_t);
            fprintf('t: %f\n', t);

%             [Jk,~] = prob.obj(x_init,p_init);
%             rho    = (Jo - Jk) / (qpv_init-qp_val);
%             if (rho > 0.75)
%                 %p_init  = p_init + step + t*(p_final - p_init);
%                 delta_t = min(gamma*delta_t, (1-t));
%             end
%             
%             if ( (rho <= 0.75) && ( 0.05 <= rho))
%                 %p_init  = p_init + step + t*(p_final - p_init);
%                 delta_t = min(delta_t, (1-t));
% 
%             end
%             
%             if (rho < 0.05)
%                 delta_t = min(0.5*delta_t, (1-t));
%             end

            %fprintf('Current parameter value: %f\n', p_init);
%             %add another stopping criteria based on x steplength
%             x_diff = x_init - x_old;
%             if ( (norm(x_diff,inf) < 1e-4) && (iter ~= 0) )
%                 break;
%             end
            
            %p_init   = p_init + step;
            qpv_init = qp_val;
            % update info
            in_info         = in_info + 1;
            info(in_info).t = t;
            info(in_info).x = x_init;
            info(in_info).y = y_init;
            
            % print out iteration and SUCCESS
            fprintf('--------------------------------------------------- \n');
            iter    = iter + 1;
            fprintf('iteration number: %d\n', iter);
            success = 1;
            if (verbose_level)
                 fprintf('%6g     %3.1e     %3.1e   %4g \n',iter, delta_t , t, success);
            end
            
            % compare with NLP
            if(nlpRun)
                %[nlpPrimal,nlpDual] = optDistillationParam(p_init);
%                 %xstart = [1.000000000000000; 0.367879441171442];
%                 xstart = [-4.000000000000000; 54.598150033144236];
%                 %p1 = p_nlp(1);
%                 %p2 = p_nlp(2);
                xstart = [0.496811686457911; 0.608467572142176];
                %xstart = x_old;
                p1 = p_init(1);
                p2 = p_init(2);
                
                % later change this to IPOPT !
%                 options = optimset('LargeScale','off','Display', 'off');
%                 [nlpPrimal, ~, ~,~, nlpDual] = fmincon(@(x)objfunhm(x,p1),xstart,[],[],[],[],[],[],@(x)confunhm(x,p2),options);
%                 diff_primal = nlpPrimal - x_init;
%                 %diff_dual   = nlpDual.ineqnonlin - y_init;
%                 diff_dual   = nlpDual.eqnonlin - y_init;
%                 %%diff_dual   = nlpDual - y_init;
                
%                 % IPOPT
%                 % objective function
%                 import casadi.*
%                 x    = SX.sym('x',2);
%                 obj  = p1*x(1)^3 + x(2)^2;
%                 cons = [exp(-x(1)) - x(2); p2 - x(1)];
%                 %lbg  = [0;0];
%                 %ubg  = [inf;inf];
%                 %ubg  = [0;0];  % test equality constraint
%                 ubg  = [0;0];
%                 lbg  = [-inf;-inf];
%                 prob1 = struct('f', obj, 'x', x, 'g', cons);
%                 solver = nlpsol('solver', 'ipopt', prob1);
%                 sol    = solver('x0', xstart, 'lbg', lbg, 'ubg', ubg);
%                 nlpPrimal = full(sol.x);
%                 %nlpDual   = abs(full(sol.lam_g));
%                 nlpDual   = full(sol.lam_g);

                % FMINCON
                options = optimset('LargeScale','off','Display', 'iter');
                [nlpPrimal, obval, exitFlag,~, nlpDual] = fmincon(@(x)objfunhm(x,p1),xstart,[],[],[],[],[],[],@(x)confunhm(x,p2),options);
                
                diff_primal = nlpPrimal - x_init;
                diff_dual   = nlpDual.ineqnonlin - y_init;
                
                err_primal  = norm(diff_primal,inf);
                err_dual    = norm(diff_dual,inf);
                fprintf('Primal error: %f\n', err_primal);
                fprintf('Dual error: %f\n', err_dual);
                %fprintf('--------------------------------------------------- \n');
                % collect data for plotting
                %t_data = [t_data; t_nlp];
                t_data = [t_data; t];
                x_data = [x_data; err_primal];
                y_data = [y_data; err_dual];
            end
            
        end
        if ( (1-t) <= 1e-5 )
            break;
        end
        
    end
end

info = qp_val;

% plot graph NLP vs. PF
if(nlpRun)
%     figure(1);
%     plot(t_data,x_data, 'LineWidth',2.5);
%     legend('primal error');
%     xlabel('parameter [t]');
%     ylabel('errors');
%     title('QP solver: Quadprog');
%     
%     figure(2);
%     plot(t_data,y_data, 'LineWidth',2.5);
%     legend('dual error');
%     xlabel('parameter [t]');
%     ylabel('errors');
%     title('QP solver: Quadprog');
    nlpResult.primal = nlpPrimal;
    nlpResult.dual   = nlpDual;
    nlpResult.t_data = t_data;
    nlpResult.x_data = x_data;
    nlpResult.y_data = y_data;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOLVE QP program
% solution of QP program: [y] (directional derivative)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y, qp_val, qp_exit, oqp, k_zero_tilde, k_plus_tilde, g, lamda] = solve_qp(prob, p, x_init, y_init, step, lb, ub)
    
    % setup objective function for QP problem
    [~,g,H,Lxp,cin,J,cp,Jeq,dpe] = prob.obj(x_init,y_init,p);
    
    % Gershgorin Convexification
    %rho = 1e3;
    %Q   = gershgorinConvex(H,rho);
    
    %f   = Lxp * step + g;
    f   = Lxp * step;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ARRANGE THE CONSTRAINT
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1. Inequality Constraint
    %    construct A and b used in QUADPROG
    A   = [];
    b   = [];
    k_zero_tilde = [];
    k_plus_tilde = [];
    if (prob.niq)
        % evaluate the constraint function and its Jacobian
        %[cin, J] = prob.cin(x_init, p);
        %cin
        
        % identify which constraint is active (strongly active) and
        % inactive (weakly active); use information for inital dual
        % variables y_init
        k_plus_tilde = find(y_init > 1e-5);                % strongly active constraint
        k_zero_tilde = find(y_init <= 1e-5);
        nk_pt = size(k_plus_tilde,1);
        nk_zt = size(k_zero_tilde,1);
        
        % the weakly active constrait assigns as inequality constrait
        % check the number of weakly active constraint
        if (nk_zt)
            if (nk_zt == prob.niq)
                % assign all the weakly active constraint as inequality 
                A = J;
                b = -(cp*step)-cin;
            else
                % assign partly of the weakly active constraint as inequality 
                A  = zeros(nk_zt, size(x_init,1));
                b  = zeros(nk_zt, 1);
                for i=1:nk_zt
                    index   = k_zero_tilde(i,1);
                    A(i,:)  = J(index,:);
                    b(i,:)  = -(cp(index,:)*step)-cin(index);
                end
            end
        else
            % no inequality constraint
            A = [];
            b = [];
        end
        
    end
    
    % 2. Equality Constraint
    %    construct Aeq and beq used in QUADPROG
    ceq = cin;  % should be set in different way
    if (prob.neq)
        Aeq  = Jeq;
        beq  = -dpe*step - ceq;  % just for case 1 parameter! need to think further!!!
        %beq  = -dpe*step;
    else
        % there might equality constraint from strongly active inequality
        % constraint; use information from k_plus_tilde
        if (k_plus_tilde)
            Aeq  = zeros(nk_pt, size(x_init,1)); 
            beq  = zeros(nk_pt, 1);
            for i=1:nk_pt
                index    = k_plus_tilde(i,1);
                Aeq(i,:) = J(index,:);
                beq(i,:) = -(cp(index,:)*step)-cin(index);
            end
        else
            Aeq = [];
            beq = [];
        end;
        
    end
    
    % prepare variables to be used for LP problem
    if (~isempty(cp))
        oqp.cp  = cp;
    else
        oqp.cp  = [];
    end
    oqp.dpe     = dpe;
    oqp.Lxp     = Lxp;
    oqp.H       = H;
    %oqp.H       = H+Q;
    if (~isempty(cin))
        oqp.cin = cin;
    else
        oqp.cin = [];
    end
    if(~isempty(J)) 
        oqp.J   = J; 
    end;
    if(~isempty(Jeq)) 
        oqp.Jeq = Jeq; 
    end;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Finally solve QP problem
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %options = optimset('Display','off','Algorithm','active-set');
    %options = optimset('Display','off','Algorithm','interior-point-convex');
    options = optimset('Display','off','Algorithm','interior-point-convex', 'MaxIter', 10);
    %options = optimset('Display','on', 'MaxIter', 2);
    %[y, qp_val, qp_exit, ~, lamda] = quadprog(H,f,A,b,Aeq,beq,lb,ub,[],options);
    startqp = tic;
    [y, qp_val, qp_exit, ~, lamda] = quadprog(H,f,A,b,Aeq,beq,lb,ub,x_init,options);  % supply initial solution
    %[y, qp_val, qp_exit, ~, lamda] = quadprog(H,f,A,b,Aeq,beq,lb,ub,x_init);
    %[y, qp_val, qp_exit, ~, lamda] = quadprog(H+Q,f,A,b,Aeq,beq,lb,ub,x_init,options);
    elapsedqp = toc(startqp);
    fprintf('QP solver runtime: %f\n',elapsedqp);
    fprintf('QP return: %d\n', qp_exit);

   
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOLVE LP program
% solution of LP program = [delta_lamda delta_eta]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [delta_l, lp_exit] = solve_lp(prob, oqp, y_init, step, y,  k_zero_tilde, k_plus_tilde, g, lamda)

    ny      = size(y_init,1);
    delta_l = zeros(ny, 1);
    nx      = size(y,1);
    lb_delta_lamda = [];
    lb_delta_tau   = [];
    a21 = [];
    a22 = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CONSTRUCT CONSTRAINT
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(~isempty(k_zero_tilde))
        % compute k_zero_hat BELONGING TO k_zero_tilde index member !!!
        nk_zt      = size(k_zero_tilde,1);
        k_zero_hat = zeros(nk_zt,1);
        for i=1:nk_zt
            index = k_zero_tilde(i);
            k_zero_hat(i) = oqp.cin(index) + oqp.J(index,:)*y + oqp.cp(index,:)*step;
        end
        
        %index_zh   = find(k_zero_hat < 0);  % 1e-5 numerical error instead of lower than zero
        index_zh   = find(k_zero_hat < 1e-5);
        nk_zh      = size(index_zh,1);
        
        remove_index = zeros(nk_zh,1);
        if (nk_zt == nk_zh)
            for i=1:nk_zt
                remove_index(i)          = k_zero_tilde(i);
                delta_l(k_zero_tilde(i)) = 0;
            end
        end
    end
    
    
    % equality constraint
    % construct matrix Aeq = [a11 a12; a21 a22];
    if (prob.neq)
        a11 = oqp.Jeq';
        
        % QUADPROG - MATLAB
        beq       = -oqp.Lxp*step - oqp.H*y - g - lamda.upper + lamda.lower;
        if (prob.niq)
            a12 = zeros(size(a11));
        else
            a12 = [];
        end
    else
        a11         = [];
        a12         = [];
        delta_lamda = [];
        lb_delta_lamda = [];
    end
    
    if (prob.niq)
        a22          = oqp.J';  %% TRANSPOSE IMPORTANT !!! REMEMBER !
        beq       = -oqp.Lxp*step - oqp.H*y - g;
        
        % add additional constraint delta_eta = 0
        if (ny > nx)
            if (nk_zh > 1)
                adc = zeros(nk_zh,size(a22,2));
                for i=1:nk_zh
                    adc(i,k_zero_tilde(i)) = 1;
                    beq = [beq; 0];
                end
                a22 = [a22;adc];
            else
                adc = zeros(1,size(a22,1));
                adc(k_zero_tilde) = 1;
                a22 = [a22;adc];
                beq = [beq; 0]; 
            end
        end
        
        
        
        % matrix arrangement
        if (prob.neq)
           a21 = zeros(size(a22));
        else
           a21 = [];
        end
    end
    
    % setup Aeq and beq
    Aeq = [a11 a12; a21 a22];
    
    lb = [];
    ub = [];
    
    % Inequality constraint
    A = [];
    b = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CONSTRUCT OBJECTIVE FUNCTION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % consist of inequality and equality parts
    f = 0;
    if (prob.niq)
        in_p = step'*oqp.cp';
        f = f + in_p;
    end
    
    if (prob.neq)
        eq_p = step'*oqp.dpe';
        f = f + eq_p;
    end
    
    % cek residual error
    x0  = lamda.eqlin;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solve LP problem; Maximization problem ! %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %optionslin = optimset('Display','off','Algorithm', 'interior-point-legacy', 'TolCon', 1e-3, 'TolFun', 1e-1, 'TolX', 1e-3, 'LargeScale','on'); %relax constraint tolerance
%     optionslin = optimset('Display','off','Algorithm', 'interior-point-legacy', 'TolCon', 1, 'TolFun', 1, 'TolX', 1, 'LargeScale','on');
%     
%     %optionslin = optimset('Display','off','Algorithm', 'interior-point', 'TolCon', 1, 'TolFun', 1, 'TolX', 1, 'LargeScale','on');
%     optionslin = optimset('Display','off', 'TolCon', 1, 'TolFun', 1, 'TolX', 1, 'LargeScale','on');
%     startlp = tic;
%     [delta_lp, ~, lp_exit] = linprog(-f,A,b,Aeq,beq,lb,ub,x0,optionslin); % supply initial guess
%     %[delta_lp, ~, lp_exit] = linprog(-f,A,b,Aeq,beq,lb,ub,x0);
%     elapsedlp = toc(startlp);
%     fprintf('LP solver runtime = %f\n',elapsedlp);
%     fprintf('LP return: %d\n', lp_exit);

    %optionslin = optimoptions(@linprog,'Algorithm','interior-point', 'ConstraintTolerance',1, 'OptimalityTolerance', 1, 'Display', 'off');
    %optionslin = optimoptions(@linprog,'Algorithm','interior-point', 'Display', 'off');
    %optionslin = optimoptions(@linprog,'Algorithm','active-set','Display', 'off', 'ConstraintTolerance',1e-6);
    %optionslin = optimoptions(@linprog,'Algorithm','active-set','Display', 'off', 'ConstraintTolerance',1e-2, 'OptimalityTolerance', 1e-2);
    optionslin = optimoptions(@linprog,'Algorithm','active-set','Display', 'off', 'ConstraintTolerance',1, 'OptimalityTolerance', 1);
    startlp = tic;
    [delta_lp, ~, lp_exit] = linprog(-f,A,b,Aeq,beq,lb,ub,x0,optionslin);
    elapsedlp = toc(startlp);
    fprintf('LP solver runtime = %f\n',elapsedlp);
    fprintf('LP return: %d\n', lp_exit);


%     % GUROBI
%     clear model;
%     model.A     = sparse(Aeq);
%     model.obj   = f;
%     model.modelsense = 'max';
%     model.rhs   = beq;
%     sign        = '=';
%     signConstraint = repmat(sign,size(Aeq,1),1);
%     model.sense = signConstraint;
% 
%     result = gurobi(model)
% 
%     disp(result.objval);
%     disp(result.x);
%     
    % arrange delta_l based on index 
    if (~isempty(k_zero_tilde))
        if (lp_exit > 0)
            ndel = 0;
            for i=1:ny
                if (any(i == remove_index))
                    %if i == remove_index
                    ndel = ndel + 1;
                    continue;  % already assigned to zero above
                else
                    ni = i - ndel;
                    delta_l(i) = delta_lp(ni);
                end
            end
        end
    else
        delta_l = delta_lp;
    end
    

end
