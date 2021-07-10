function [prob] = distColACstr_pn(p)
%DISTCOLA_CASADI_PN Summary of this function goes here
% 
% [OUTPUTARGS] = DISTCOLA_CASADI_PN(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2016/04/07 19:01:23 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2016

import casadi.*

prob = struct('neq',{0},'niq',{0},'cin',{0},'ceq',{0},'dp_in',{0},'dp_eq',{0},'hess',{0},'lxp',{0},'x',0,'name',0);
prob.neq  = 2000;         % HARD-CODE !    % number of equality constraint
prob.niq  = 0;            % number of inequality constraint
prob.name = 'Distillation Column A model';
prob.x    = zeros(2,1);

prob.obj  = (@(x,y,p,N)(objective(x,y,p,N)));

end

function [f,g,H,Lxp,cst,J,cp,Jeq,dpe,Hobj] = objective(x,y,p,N)
    import casadi.*
   
    nPrimal = numel(x);
    %nDual   = numel(y);
    nDual   = numel(y.lam_g);
    nParam  = numel(p);
    % model parameters
    NT = 41;
    %Uf = 1.1;              % Feeding rate F
    Uf = 0.3;  
    %P  = SX.sym('P');      % parameter used in path-following algorithm
    % invoke the model
    %[~,state,xdot,inputs] = DistColAP(Uf,p);
    %[~,state,xdot,inputs] = DistColA(Uf);
    [~,state,xdot,inputs] = DistColACstr(Uf);
    sf = Function('sf',{state,inputs}, {xdot});
    
    % dimensions
    global nx nu nk d tf ns;
%     nx = 82;
%     nu = 2;
% 
%     nk = 10;       % control discretization
%     tf = 10.0;     % end time
    h  = tf/nk;
    
    % preparing collocation matrices
    [~,C,D,d] = collocationSetup();
    

    % NLP variable vector 
    V    = {};      % decision variables contain both control and state variables
    obj  = 0;       % objective function
    cons = {};      % nonlinear constraint
    
    delta_time = 1;
    alpha = 1;
    beta  = 1;
    gamma = 1;

    % Initial states and controls 
    %load u_opt_ss.mat;
    %xdot_val_rf_ss = xf;
    load CstrDistXinit.mat;
%     Xinit(85) = sf1*Xinit(85);
%     Xinit(86) = sf2*Xinit(86);
%     global sf1 sf2;
    %xdot_val_rf_ss = Xinit(1:84);
    xf             = Xinit(1:84);
    u_opt          = Xinit(85:89);
    
    % prices
    pf = 1; 
    %pV = 0.01;
    pV = 0.02;
    pB = 2; 
    pD = 0;
    %F  = Uf;
    
    % compact price variable
    params.price.pf = pf;
    params.price.pV = pV;
    params.price.pB = pB;
    params.price.pD = pD;
    params.price.F_0= Uf;
    
    % controller gains
    KcB = 10;  
    KcD = 10;
    % Nominal holdups - these are rather small 
    MDs = 0.5; 
    MBs = 0.5;         
    % Nominal flows
    Ds  = 0.5; 
    Bs  = 0.5;
    
    % compact controller gain variable
    params.gain.KcB = KcB;
    params.gain.KcD = KcD;
    params.gain.MDs = MDs;
    params.gain.MBs = MBs;
    params.gain.Ds  = Ds;
    params.gain.Bs  = Bs;
    
    % compact model variable
    params.model.NT = NT;
    params.model.sf  = sf;
    params.model.xdot_val_rf_ss = xf;
    params.model.x  = x;
    params.model.u_opt = u_opt;
    
    % compact collocation variable
    params.colloc.C = C;
    params.colloc.D = D;
    params.colloc.h = h;
    
    % compact weight variable
    params.weight.delta_time = delta_time;
    params.weight.alpha      = alpha;
    params.weight.beta       = beta;
    params.weight.gamma      = gamma;
    
    % "Lift" initial conditions
    %X0   = SX.sym('X0', nx);
    X0  = MX.sym('X0', nx);
    V    = {V{:}, X0};
    %cons = {cons{:}, X0 - xdot_val_rf_ss};
    cons = {cons{:}, X0 - x(1:nx,1)};
    cons_x0 = X0 - x(1:nx,1);
    
    %cons = {cons{:}, X0 - p};
    %cons_x0 = X0 - p;
    
    % formulate the NLP
    Xk = X0;
    
%     load weights2.mat;
%     params.weights = weights;

    load Qmax.mat;
    params.Qmax = Qmax;
    %load Qmax_n1.mat;  % with measurement noise 1 percent
    %params.Qmax = Qmax_n1;
    %load Qmax_31.mat;
    %params.Qmax = Qmax_31;
    % objective function at initial state value
    %obj   = obj +( Qmax(1:nx,1).*(X0 - xdot_val_rf_ss))' * (X0 - xdot_val_rf_ss);
    ssoftc = 0;
    for i=1:N
        [obj,cons,V,Xk,params,ssoftc] = iterateOnPredictionHorizon(Xk, V, cons, obj, params, i,ssoftc);
    end
    
    V = vertcat(V{:});
    % Concatenate constraints
    cons  = vertcat(cons{:});
    
    % objective function and constraint functions
    f = Function('f', {V}, {obj}, char('V'), char('objective'));
    c = Function('c', {V}, {cons}, char('V'), char('constraint'));
    cx0 = Function('cx0', {X0}, {cons_x0}, char('X0'), char('constraint'));
    
    % construct Lagrangian
    %lag_expr = obj + y'*cons;
    lag_expr = obj + y.lam_g'*cons;
    %lag_expr = obj - y.lam_g'*cons;
    
    g    = f.gradient();
    lagr = Function('lagr', {V}, {lag_expr}, char('V'), char('lag_expr'));
    H    = Function(lagr.hessian('V','lag_expr'));
    Hobj = f.hessian('V','objective');
    J    = c.jacobian('V','constraint');
    Jp   = cx0.jacobian('X0','constraint');
    
    f   = f(x);
    g   = g(x);
    H   = H(x);
    Lxp = H(1:nPrimal,1:nParam);
    J   = J(x);
    %cp  = J(1:nDual,1:nParam);
    Jtemp  = zeros(nDual,nParam);
    cp  = Jp(x(1:nParam));
    Jtemp(1:nParam,1:nParam) = full(cp);
    cp  = sparse(Jtemp);
    cst = c(x);
    
    % Evaluation of objective function used for Greshgorin bound
    Hobj = Hobj(x);
    Hobj = sparse(Hobj);

    f   = full(f);
    g   = sparse(g);    
    H   = sparse(H);
    Lxp = sparse(Lxp);
    J   = sparse(J);
    cp  = sparse(cp);
    cst = full(cst);

    % Equality constraint
    Jeq = J;
    dpe = cp;
    
    
end

function [obj,cons,V,Xk,params,ssoftc] = iterateOnPredictionHorizon(Xk, V, cons, obj, params, iter, ssoftc)

   import casadi.*
   
   % extract compact variables
   NT  = params.model.NT;
   sf  = params.model.sf;
   xdot_val_rf_ss = params.model.xdot_val_rf_ss;
   u_opt = params.model.u_opt;
   
   pf = params.price.pf;
   pV = params.price.pV;
   pB = params.price.pB;
   pD = params.price.pD;
   F_0= params.price.F_0;
   
   KcB = params.gain.KcB;
   KcD = params.gain.KcD;
   MDs = params.gain.MDs;
   MBs = params.gain.MBs;
   Ds  = params.gain.Ds;
   Bs  = params.gain.Bs;
   
   C = params.colloc.C;
   D = params.colloc.D;
   h = params.colloc.h;
   
   delta_time = params.weight.delta_time;
   Qmax = params.Qmax;
%    alpha      = params.weight.alpha;
%    beta       = params.weight.beta;
%    gamma      = params.weight.gamma;
   
   global nx nu nk d ns N;
   %weights = ones(nx+nu,1);
   %params.weights = [ones(nu,1);ones(nx,1)];  
   %load weights.mat;
   %load weights1.mat;

   count = 1;
   for k=0:nk-1
        % New NLP variable for the control
        %Uk  = MX.sym(['U_' num2str(k)], nu);
        Uk  = MX.sym(['U_' num2str((iter-1)*nk+k)], nu);
        V   = {V{:}, Uk};
        
        %J   = J + (Qmax(nx+nu,1).*(Uk - u_opt))' * (Uk - u_opt);
        Jcontrol   = (Qmax(nx+1:nx+nu,1).*(Uk - u_opt))' * (Uk - u_opt);

        % State at collocation points
        Xkj = {};
        %AvgX43 = 0;
        %AvgX83 = 0;
        SumX1 = 0;
        for j=1:d
            %Xkj{j} = MX.sym(['X_' num2str(k) '_' num2str(j)], nx);
            Xkj{j} = MX.sym(['X_' num2str((iter-1)*nk+k) '_' num2str(j)], nx);
            V      = {V{:}, Xkj{j}};
            
            % sum of bottom concentration
            %SumX1 = SumX1 + Xkj{j}(1,1);
            
%             if j == 1
%                 SumX1 = SumX1 + 0.2*Xkj{j}(1,1);
%             end
%             
%             if j == 2
%                 SumX1 = SumX1 + 0.4*Xkj{j}(1,1);
%             end
%             
%             if j == 3
%                 SumX1 = SumX1 + 0.4*Xkj{j}(1,1);
%             end
            

            
%              % add soft-constraint here
%             % introduce slack variable (pikirkan lagi!)
%             Skj{j} = MX.sym(['S_' num2str((iter-1)) '_' num2str(j)], ns);
%             V      = {V{:}, Skj{j}};
%             %cons   = {cons{:}, [Xkj{j}(1,1) - Skj{j}(1,1); Xkj{j}(84,1) - Skj{j}(2,1)] };
%             %cons   = {cons{:}, [Xkj{j}(1,1) - Skj{j}(1,1)]};
%             cons   = {cons{:}, [Xkj{j}(1,1) - Skj{j}(1,1); Xkj{j}(83,1) + Skj{j}(2,1); Xkj{j}(84,1) - Skj{j}(3,1)] };
%             ssoftc = ssoftc + sum(Skj{j});
           
%             % equality constraints for level controllers
%             cons = {cons{:}, [Xkj{j}(43,1) - 0.5; Xkj{j}(83,1) - 0.5]};

%             % compute average holdup 
%             AvgX43 = AvgX43 + Xkj{j}(43,1);
%             AvgX83 = AvgX83 + Xkj{j}(83,1);

        end

        % Loop over collocation points
        Xk_end = D(1)*Xk; 
%         AvgX43 = AvgX43 / d; 
%         AvgX83 = AvgX83 / d;
        for j=1:d
           % Expression for the state derivative at the collocation point
           xp = C(1,j+1)*Xk;
           for r=1:d
               xp = xp + C(r+1,j+1)*Xkj{r};
           end

           % Append collocation equations
           fj   = sf(Xkj{j},Uk);
           cons = {cons{:}, h*fj - xp};

           % Add contribution to the end state
           Xk_end = Xk_end + D(j+1)*Xkj{j};

%            % average holdup constraints
%            cons = {cons{:}, [Xkj{j}(43,1) - AvgX43; Xkj{j}(83,1) - AvgX83] };

        end    

        % New NLP variable for state at end of interval
        %Xk  = MX.sym(['X_' num2str(k+1)], nx);
        Xk  = MX.sym(['X_' num2str((iter-1)*nk+k)], nx);
        V   = {V{:}, Xk};

        % Add equality constraint
        cons= {cons{:}, Xk_end-Xk};
        %cons   = {cons{:}, xdot_val_rf_ss-Xk};

%         % Add terminal constraint
%         if iter == N
%             cons   = {cons{:}, xdot_val_rf_ss-Xk};
%         else
%             % Add equality constraint
%             cons   = {cons{:}, Xk_end-Xk};
%         end
        
%         % bottom concentration
%         cons= {cons{:}, SumX1 + Xk(1,1)};
        
%         % active constraint
%         cons= {cons{:}, [Xk(84,1)]};
        
%         % average holdup constraints
%         cons = {cons{:}, [Xk(43,1) - AvgX43; Xk(83,1) - AvgX83] };
        
%         % equality constraints for level controllers
%         cons = {cons{:}, [Xk(43,1) - 0.5; Xk(83,1) - 0.5]};
        
%         % Add soft-constraint here
%         % introduce slack variable (pikirkan lagi!)
%         Sk = MX.sym(['S_' num2str((iter-1))], ns);
%         V  = {V{:}, Sk};
%         %cons = {cons{:}, [Xk(1,1) - Sk(1,1); Xk(84,1) - Sk(2,1)] };
%         %cons = {cons{:}, [Xk(1,1) - Sk(1,1)] };
%         cons = {cons{:}, [Xk(1,1) - Sk(1,1); Xk(83,1) + Sk(2,1); Xk(84,1) - Sk(3,1)] };
%         ssoftc = ssoftc + sum(Sk);

        Jecon  = (pf*F_0 + pV*Uk(2) - pB*Uk(5) - pD*Uk(4)) * delta_time;
        %Qmax(43,1) = 2*Qmax(43,1);
        %Qmax(83,1) = 2*Qmax(83,1);
        Jstate =(Qmax(1:nx,1).*(Xk - xdot_val_rf_ss))' * (Xk - xdot_val_rf_ss) * delta_time;
        
        % compute rotated cost function
        fm  = sf(Xk,Uk);
        % load Lagrange multipliers from steady-state optimization
        load LamdaCstrDist.mat; % lamda
        Jmodel = lamda'*fm;
        
        alpha  = 1;
        beta   = 1;
        gamma  = 1;
        %rho    = 1e3;
        rho    = 0.5;
        %rho    = 5e-2;
        
        %J = J + alpha*Jcontrol + gamma*Jstate + beta*Jecon;
        obj = obj + alpha*Jcontrol + gamma*Jstate + beta*Jecon;
        %obj = obj + alpha*Jcontrol + gamma*Jstate + beta*Jecon + Jmodel;
        %obj = obj + alpha*Jcontrol + gamma*Jstate + beta*Jecon + Jmodel + rho*ssoftc;
        %obj = obj + beta*Jecon;
        %obj = obj + alpha*Jcontrol*0.01 + gamma*Jstate*0.01 + beta*Jecon;
        %obj = obj + alpha*Jcontrol + gamma*Jstate;
        
        % add state soft-constraint
        %rho = 1e3;
        %J = J + rho*epsilon;

    end
end

