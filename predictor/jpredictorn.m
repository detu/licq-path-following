function [x_init,y_init,info] = jpredictorn(problem, p_init, p_final, x_init, y_init, delta_t)
%JPREDICTORN Summary of this function goes here
% 
% Modification of jpredictor with tolerance in LP problem
% epsilon1 and epsilon2
% [OUTPUTARGS] = JPREDICTORN(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2015/10/23 15:14:51 $	$Revision: 0.1 $
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
gamma   = 1.25;

% if (delta_t == 1)
%     t = 1;
% else
%     t = 0;
% end

% display problem
%fprintf('Solving problem %s with p=%3.1e using estimate at p=%3.1e\n',prob.name,p(1),p(2));
fprintf('Solving problem %s \n',prob.name);
fprintf('iteration  delta_t        t        Success\n');
while ( (t <= 10) && (delta_t > 1e-6) )   % add stopping criteria delta_t not equal zero 
%while ( (t <= 2) )
    
    % calculate step s
    step = delta_t * (p_final - p_init);
    
    % solve QP problem
    % output from QP is y (directional derivative)
    %[y, qp_exit, oqp, k_zero_tilde] = solve_qp(prob, p, x_init, y_init, step);
    [y, qp_exit, oqp, k_zero_tilde] = solve_qp(prob, p_init, x_init, y_init, step);
    
    if (qp_exit < 0); % QP is infeasible
        
        % shorten step
        delta_t = alpha_2 * t;
        %t       = t - alpha_2 * delta_t;
        t       = t - delta_t;
        
        % print out iteration and FAIL
        iter    = iter + 1;
        success = 0;
        fprintf('%6g     %3.1e     %3.1e   %4g \n',iter, delta_t , t, success);
        
        
    else % QP is feasible
        % solve LP
        [delta_l, lp_exit] = solve_lp(prob, oqp, y_init, step, y, k_zero_tilde);
        
        % check exit flag
        if (lp_exit < 0) % LP is infeasible
            
            % shorten step
            delta_t = alpha_1 * t;
            %t       = t - alpha_1 * delta_t;
            t       = t - delta_t;
            
            % print out iteration and FAIL
            iter    = iter + 1;
            success = 0;
            fprintf('%6g     %3.1e     %3.1e   %4g \n',iter, delta_t , t, success);
            
        else
            
            % update states, multipliers, parameter, and time step
            x_init  = x_init + y;
            y_init  = y_init + delta_l;
            t       = t + delta_t;
            %delta_t = 1 - t; % fix delta_t
            %delta_t = t;
            %delta_t  = min(gamma*delta_t, (1-t));
            
            p_init   = p_init + step;
            
            % update info
            in_info         = in_info + 1;
            info(in_info).t = t;
            info(in_info).x = x_init;
            info(in_info).y = y_init;
            
            % print out iteration and SUCCESS
            iter    = iter + 1;
            success = 1;
            fprintf('%6g     %3.1e     %3.1e   %4g \n',iter, delta_t , t, success);
        end
        
        
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOLVE QP program
% solution of QP program: [y] (directional derivative)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y, qp_exit, oqp, k_zero_tilde] = solve_qp(prob, p, x_init, y_init, step)
    
    % setup objective function for QP problem
    H   = prob.hess(x_init,y_init,p);  
    Lxp = prob.lxp(x_init,y_init,p);
    f   = Lxp * step;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ARRANGE THE CONSTRAINT
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1. Inequality Constraint
    %    construct A and b used in QUADPROG
    A   = [];
    b   = [];
    cp  = [];
    cin = [];
    J   = [];
    k_zero_tilde = [];
    if (prob.niq)
        % evaluate the constraint function and its Jacobian
        [cin, J] = prob.cin(x_init, p);
        
        % identify which constraint is active (strongly active) and
        % inactive (weakly active); use information for inital dual
        % variables y_init
        k_plus_tilde = find(y_init > 1e-5);                % strongly active constraint
        %k_zero_tilde = find(y_init <= 1e-5 & y_init >= 0); % weakly active
        k_zero_tilde = find(y_init <= 1e-5);
        nk_pt = size(k_plus_tilde,1);
        nk_zt = size(k_zero_tilde,1);
        
        % the weakly active constrait assigns as inequality constrait
        % check the number of weakly active constraint
        cp = prob.dp_in(x_init,p);
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
    Jeq = prob.ceq(x_init, p);
    dpe = prob.dp_eq(x_init,p);
    if (prob.neq)
        Aeq = Jeq;
        %beq = dpe;
        beq = -dpe'*step;
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
    %options = optimset('Display','on','Algorithm','active-set');
    %[y, ~, qp_exit] = quadprog(H,f,A,b,Aeq,beq,[],[],x_init,options); % is it correct to use x_init as initial guess ?
    options = optimset('Display','off','Algorithm','interior-point-convex');
    [y, ~, qp_exit, ~, lamda] = quadprog(H,f,A,b,Aeq,beq,[],[],[],options);
    lamda.eqlin
    %[y, ~, qp_exit] = quadprog(H,f,A,b,Aeq,beq,[],[],x_init,options);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOLVE LP program
% solution of LP program = [delta_lamda delta_eta]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [delta_l, lp_exit] = solve_lp(prob, oqp, y_init, step, y,  k_zero_tilde)

    % bound tolerance in LP program for allowing active-set changes
    %eps1 = 0.01;
    eps2 = inf;
    eps3 = 1e-4;

    ny      = size(y_init,1);
    delta_l = zeros(ny, 1);
    lb_delta_lamda = [];
    lb_delta_tau   = [];
    ub_delta_lamda = [];
    ub_delta_tau   = [];
    remove_index   = [];
    
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
        % does not need to eliminate remove_index 
        %na22 = prob.niq - size(remove_index,1);
        na22 = prob.niq;
        a22  = zeros(na22,na22);
        lb_delta_tau = zeros(na22,1);
        ub_delta_tau = zeros(na22,1);
        beq  = zeros(na22,1);
        
        for i=1:prob.niq
            if (any(i == remove_index))
                if (na22 == 0)
                    a22       = [];
                    delta_eta = [];
                    break;
                else
                    invJ            = oqp.J';
                    a22(:,i)        = invJ(:,i);
                    lb_delta_tau(i) = -eps2;
                    ub_delta_tau(i) = eps2;
                    beq(i)          = -oqp.Lxp(i,:)*step - oqp.H(i,:)*y;
                end
                %continue;
            else
                invJ                = oqp.J';
                if na22 > 1
                    a22(:,i)        = invJ(:,i);
                    %lb_delta_tau(i) = -y_init(i);
                    %lb_delta_tau(i) = -y_init(i) - eps3;  % RELAX THE UPPER BOUND
                    lb_delta_tau(i) = -inf;
                    ub_delta_tau(i) = inf;
                    beq(i)          = -oqp.Lxp(i,:)*step - oqp.H(i,:)*y;
                else
                    % find an index from absoluture y vector
                    yabs = abs(y);
                    idy  = find(yabs == max(yabs), 1);
                    a22  = invJ(idy, i);
                    %lb_delta_tau = -y_init(i);  % from the constraint
                    %lb_delta_tau(i) = -y_init(i) - eps3;  % RELAX THE UPPER BOUND
                    lb_delta_tau(i) = -inf;
                    ub_delta_tau(i) = inf;
                    beq             = -oqp.Lxp(idy,:)*step - oqp.H(idy,:)*y;
                end
            end
        end    
        
        if (prob.neq)
           a21 = zeros(size(a22));
        else
           a21 = [];
        end

    else
        a22       = [];
        a21       = [];
        delta_eta = [];
        beq       = -oqp.Lxp*step - oqp.H*y;
    end
    
    % setup Aeq and beq
    Aeq = [a11 a12; a21 a22];
    
    % check number of zero element in delta_eta
    if (~isempty(k_zero_tilde))
        num_zero = size(index_zh,1);
        if (num_zero == prob.niq)
            delta_l = [ delta_lamda; delta_eta];
            %lp_exit = [];
            lp_exit = -1; % just put negative value; error!
            return % do not run optimizer linprog
        end
    end
    
    % bound constraint THINK !
    % need to check size of Aeq and beq
    %n_aeq = size(Aeq,1); 
    %ub = inf*ones(n_aeq, 1);
    ub = [ub_delta_lamda; ub_delta_tau];
    lb = [lb_delta_lamda; lb_delta_tau];
    
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
        %f(remove_index) = [];
    end
    
    if (prob.neq)
        %eq_p = oqp.dpe'*step;
        %eq_p = step'*oqp.dpe;
        eq_p = step'*oqp.dpe';
        f = f + eq_p;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solve LP problem; Maximization problem ! %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     optionslin = optimset('Display','off','Algorithm', 'active-set');
%     [delta_l, ~, lp_exit] = linprog(f,A,b,Aeq,beq,lb,ub,[],optionslin);
    optionslin = optimset('Display','off','Algorithm', 'interior-point');
    [delta_lp, ~, lp_exit] = linprog(-f,A,b,Aeq,beq,lb,ub,[],optionslin);
    delta_lp
    %[delta_lp, ~, lp_exit] = linprog(f,A,b,Aeq,beq,lb,ub,y_init(1),optionslin);
    
    % arrange delta_l based on index
    delta_l = delta_lp;
%     if (~isempty(k_zero_tilde))
%         if (lp_exit > 0)
%             ndel = 0;
%             for i=1:ny
%                 if (any(i == remove_index))
%                     %if i == remove_index
%                     ndel = ndel + 1;
%                     continue;  % already assigned to zero above
%                 else
%                     ni = i - ndel;
%                     delta_l(i) = delta_lp(ni);
%                 end
%             end
%         end
%     else
%         delta_l = delta_lp;
%     end
    

end