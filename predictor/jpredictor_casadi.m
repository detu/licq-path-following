function [x_init,y_init,info] = jpredictor_casadi(problem, p_init, p_final, x_init, y_init, delta_t, lb_init, ub_init, verbose_level)
%JPREDICTOR_CASADI Summary of this function goes here
% 
% [OUTPUTARGS] = JPREDICTOR_CASADI(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: detu $	$Date: 2015/11/15 18:31:13 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2015


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


% display problem
%fprintf('Solving problem %s with p=%3.1e using estimate at p=%3.1e\n',prob.name,p(1),p(2));
if (verbose_level)
    fprintf('Solving problem %s \n',prob.name);
    fprintf('iteration  delta_t        t        Success\n');
end

while (t < 1)
%while ( (t <= t_final) && (delta_t > 1e-6) )   % add stopping criteria delta_t not equal zero
    
    %if (iter == 3) break; end
    
    % calculate step s
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
    %[y, qp_exit, oqp, k_zero_tilde] = solve_qp(prob, p, x_init, y_init, step);
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
            %y_init  = y_init + delta_l;
            y_init  = delta_l;
            t       = t + delta_t;
            %p_init  = p_init + step + t*(p_final - p_init);
            p_init  = (1 - t)*p_init + t*p_final;
            
            % FIX STEPLENGTH
            delta_t = min(delta_t, (1-t));
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
            %delta_t
            %p_init
            fprintf('Current parameter value: %f\n', p_init);
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
            %info(in_info).t = t / t_final;
            info(in_info).x = x_init;
            info(in_info).y = y_init;
            
            % print out iteration and SUCCESS
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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOLVE QP program
% solution of QP program: [y] (directional derivative)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y, qp_val, qp_exit, oqp, k_zero_tilde, k_plus_tilde, g, lamda] = solve_qp(prob, p, x_init, y_init, step, lb, ub)
    
    % setup objective function for QP problem
    %H   = prob.hess(x_init,y_init,p);  
    %Lxp = prob.lxp(x_init,y_init,p);
    %f   = Lxp * step;
    % modification of objective function is here !
    %[~,g] = prob.obj(x_init,p);
    %[~,g,H] = prob.obj(x_init,p);
    %f     = g'*step;
    %f     = g;
    [~,g,H,Lxp,cin,J,cp,Jeq,dpe] = prob.obj(x_init,y_init,p);
    f   = Lxp * step + g;
    %f   = Lxp * step;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ARRANGE THE CONSTRAINT
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1. Inequality Constraint
    %    construct A and b used in QUADPROG
    A   = [];
    b   = [];
    %cp  = [];
    %cin = [];
    %J   = [];
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
        %cp = prob.dp_in(x_init,p);
        if (nk_zt)
            if (nk_zt == prob.niq)
                % assign all the weakly active constraint as inequality 
                A = J;
                b = -(cp*step)-cin;
            else
                % assign partly of the weakly active constraint as inequality 
                A  = zeros(nk_zt, size(x_init,1));
                %pd = zeros(nk_zt, size(x_init,1));
                b  = zeros(nk_zt, 1);
                for i=1:nk_zt
                    index   = k_zero_tilde(i,1);
                    A(i,:)  = J(index,:);
                    %pd(i,:) = cp(index,:);
                    %b(i,:)  = -(pd(i,:)*step)-cin(index);
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
    %Jeq = prob.ceq(x_init, p);
    %dpe = prob.dp_eq(x_init,p);
    if (prob.neq)
        Aeq = Jeq;
        %beq = dpe;
        %beq = -dpe'*step;
        beq = -dpe'*step - cin;
    else
        % there might equality constraint from strongly active inequality
        % constraint; use information from k_plus_tilde
        if (k_plus_tilde)
            Aeq  = zeros(nk_pt, size(x_init,1)); 
            %pp   = zeros(nk_pt, size(x_init,1));
            beq  = zeros(nk_pt, 1);
            for i=1:nk_pt
                index    = k_plus_tilde(i,1);
                Aeq(i,:) = J(index,:);
                %pp(i,:)  = cp(index,:);
                %beq = -(cp(index,:)*step);
                %beq(i,:) = -(cp(index,:)*step);
                beq(i,:)  = -(cp(index,:)*step)-cin(index);
                %cin(index)
            end
            %beq = -(pp*step);
        else
            Aeq = [];
            beq = [];
        end;
        
    end
    
    % prepare variables to be used for LP problem
    %oqp.cp      = cp;
    if (~isempty(cp))
        oqp.cp  = cp;
    else
        oqp.cp  = [];
    end
    oqp.dpe     = dpe;
    oqp.Lxp     = Lxp;
    oqp.H       = H;
    %oqp.cin     = cin;
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
    options = optimset('Display','off','Algorithm','interior-point-convex');
    %[y, qp_val, qp_exit, ~, lamda] = quadprog(H,f,A,b,Aeq,beq,lb,ub,[],options);
    [y, qp_val, qp_exit, ~, lamda] = quadprog(H,f,A,b,Aeq,beq,lb,ub,x_init,options);  % supply initial solution
    %qp_exit
    %if ( qp_exit >= 0 )
    %    lamda.eqlin
    %    lamda.ineqlin
    %end
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
            %k_zero_hat(i) = oqp.cin(index) + oqp.J(index,:)'*y + oqp.cp(index,:)'*step;
            k_zero_hat(i) = oqp.cin(index) + oqp.J(index,:)*y + oqp.cp(index,:)*step;
        end
        
        %k_zero_hat = oqp.cin + oqp.J'*y + oqp.cp'*step;
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
        %beq       = -oqp.Lxp*step - oqp.H*y - g;
        %beq       = -oqp.Lxp*step - oqp.H*y;
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
        %a22          = oqp.J;
        a22          = oqp.J';  %% TRANSPOSE IMPORTANT !!! REMEMBER !
        %beq          = -oqp.H*y - g;
        beq       = -oqp.Lxp*step - oqp.H*y - g;
%         lb_delta_tau = -y_init;
        
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
    
    % bound constraint THINK !
    % need to check size of Aeq and beq
    %n_aeq = size(Aeq,1); 
    %ub = inf*ones(n_aeq, 1);
    %ub = inf*ones(ny, 1);
    %lb = [lb_delta_lamda; lb_delta_tau];
    %lb = -inf*ones(ny, 1);
    %lb = zeros(ny, 1);
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
        %in_p = oqp.cp'*step;
        %in_p = oqp.cp*step;
        in_p = step'*oqp.cp';
        f = f + in_p;
        %f = f + in_p + (oqp.cp * y_init)'*step;
        %if (~no_remove)
        %f(remove_index) = [];
        %end
    end
    
    if (prob.neq)
        %eq_p = oqp.dpe'*step;
        %eq_p = step'*oqp.dpe;
        eq_p = step'*oqp.dpe';
        f = f + eq_p;
    end
    
    % cek residual error
    x0  = lamda.eqlin;
%     res = Aeq*x0 - beq;
%     err_res = norm(res,inf);
%     fprintf('residual error: %f\n', err_res);
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solve LP problem; Maximization problem ! %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     optionslin = optimset('Display','off','Algorithm', 'active-set');;
    %optionslin = optimset('Display','off','Algorithm', 'interior-point');
    %optionslin = optimset('Display','off','Algorithm', 'interior-point', 'TolCon', 1e-3, 'TolFun', 1e-1, 'TolX', 1e-3, 'LargeScale','on'); %relax constraint tolerance
    optionslin = optimset('Display','off','Algorithm', 'interior-point', 'TolCon', 1e-3, 'TolFun', 1, 'TolX', 1e-3, 'LargeScale','on');
    %[delta_lp, ~, lp_exit] = linprog(-f,A,b,Aeq,beq,lb,ub,[],optionslin);
    [delta_lp, ~, lp_exit] = linprog(-f,A,b,Aeq,beq,lb,ub,x0,optionslin); % supply initial guess
    %[delta_lp, ~, lp_exit] = linprog(-f,A,b,Aeq,beq,lb,ub,[],optionslin);
    %lp_exit
    %delta_lp
    %delta_lp(k_zero_tilde) = 0;
    fprintf('LP return: %d\n', lp_exit);

   
    fprintf('--------------------------------------------------- \n');
    
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