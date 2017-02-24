function [prob] = distColA_casadi_pn(p)
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
    nDual   = numel(y);
    nParam  = numel(p);
    % model parameters
    NT = 41;
    Uf = 1.2;              % Feeding rate F
    %P  = SX.sym('P');      % parameter used in path-following algorithm
    % invoke the model
    %[~,state,xdot,inputs] = DistColAP(Uf,p);
    [~,state,xdot,inputs] = DistColA(Uf);
    sf = Function('sf',{state,inputs}, {xdot});
    
    %     % soft-constraint for state variables
     epsilon = 0.1; 
    
    % dimensions
    global nx nu nk d tf;
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
    load u_opt_ss.mat;
    xdot_val_rf_ss = xf;
    
    % prices
    pf = 1; 
    pV = 0.0001;
    pB = 1; 
    pD = 2;
    F  = Uf;
    
    % compact price variable
    params.price.pf = pf;
    params.price.pV = pV;
    params.price.pB = pB;
    params.price.pD = pD;
    params.price.F  = F;
    
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
    
    % formulate the NLP
    Xk = X0;
    
%     load weights2.mat;
%     params.weights = weights;

    load Qmax.mat;
    params.Qmax = Qmax;
    
    % objective function at initial state value
    obj   = obj +( Qmax(1:nx,1).*(X0 - xdot_val_rf_ss))' * (X0 - xdot_val_rf_ss);
    
    for i=1:N
        [obj,cons,V,Xk,params] = iterateOnPredictionHorizon(Xk, V, cons, obj, params, i, epsilon);
    end
    
    V = vertcat(V{:});
    % Concatenate constraints
    cons  = vertcat(cons{:});
    
    % objective function and constraint functions
    f = Function('f', {V}, {obj}, char('V'), char('objective'));
    c = Function('c', {V}, {cons}, char('V'), char('constraint'));
    
    % construct Lagrangian
    lag_expr = obj + y'*cons;
    
    g    = f.gradient();
    lagr = Function('lagr', {V}, {lag_expr}, char('V'), char('lag_expr'));
    H    = Function(lagr.hessian('V','lag_expr'));
    Hobj = f.hessian('V','objective');
    J    = c.jacobian('V','constraint');
    
    f   = f(x);
    g   = g(x);
    H   = H(x);
    Lxp = H(1:nPrimal,1:nParam);
    J   = J(x);
    cp  = J(1:nDual,1:nParam);
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

function [J,g,w,Xk,params] = iterateOnPredictionHorizon(Xk, w, g, J, params, iter, epsilon)

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
   F  = params.price.F;
   
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
   
   global nx nu nk d;
   %weights = ones(nx+nu,1);
   %params.weights = [ones(nu,1);ones(nx,1)];  
   %load weights.mat;
   %load weights1.mat;

   count = 1;
   for k=0:nk-1
        % New NLP variable for the control
        %Uk  = MX.sym(['U_' num2str(k)], nu);
        Uk  = MX.sym(['U_' num2str((iter-1)*nk+k)], nu);
        w   = {w{:}, Uk};
        
        J   = J + (Qmax(nx+nu,1).*(Uk - u_opt))' * (Uk - u_opt);

        % State at collocation points
        Xkj = {};
        for j=1:d
            %Xkj{j} = MX.sym(['X_' num2str(k) '_' num2str(j)], nx);
            Xkj{j} = MX.sym(['X_' num2str((iter-1)*nk+k) '_' num2str(j)], nx);
            w      = {w{:}, Xkj{j}};

        end

        % Loop over collocation points
        Xk_end = D(1)*Xk; 
        for j=1:d
           % Expression for the state derivative at the collocation point
           xp = C(1,j+1)*Xk;
           for r=1:d
               xp = xp + C(r+1,j+1)*Xkj{r};
           end

           % Append collocation equations
           fj  = sf(Xkj{j},Uk);
           g   = {g{:}, h*fj - xp};

           % Add contribution to the end state
           Xk_end = Xk_end + D(j+1)*Xkj{j};

           % Add contribution to quadrature function
           %J = J + B(j+1)*qj*h;
                   
           % objective function
           MB = Xkj{j}(NT+1);  
           MD = Xkj{j}(2*NT);          % Actual reboiler and condenser holdup
           Do = Ds+(MD-MDs)*KcD;       % Distillate flow
           Bo = Bs+(MB-MBs)*KcB;       % Bottoms flow
           
           J = J + (pf*F + pV*Uk(2) - pB*Bo - pD*Do) * delta_time;
%            J = J + (Uk - u_opt)' *(Uk - u_opt) * delta_time;
%            J = J + (Xkj{j} - xdot_val_rf_ss)' *(Xkj{j} - xdot_val_rf_ss) * delta_time;
           
%            %params.weights(count:count+nx-1,1)
%            J = J + (params.weights(count:count+nx-1,1).*(Xkj{j} - xdot_val_rf_ss))' *(Xkj{j} - xdot_val_rf_ss)*delta_time;
%            count = count + nx;
%            %params.weights(count:count+nu-1,1)
%            J = J + (Uk - u_opt)' * (params.weights(count:count+nu-1,1).*(Uk - u_opt))*delta_time;
%            count = count + nu;
           J  = J + (Qmax(1:nx,1).*(Xkj{j} - xdot_val_rf_ss))' * (Xkj{j} - xdot_val_rf_ss) * delta_time;
        end    

        % New NLP variable for state at end of interval
        %Xk  = MX.sym(['X_' num2str(k+1)], nx);
        Xk  = MX.sym(['X_' num2str((iter-1)*nk+k)], nx);
        w   = {w{:}, Xk};

        % Add equality constraint
        g   = {g{:}, Xk_end-Xk};
        
        
        % objective function
        MB = Xk(NT+1);  
        MD = Xk(2*NT);              % Actual reboiler and condenser holdup
        Do = Ds+(MD-MDs)*KcD;       % Distillate flow
        Bo = Bs+(MB-MBs)*KcB;       % Bottoms flow
        
        J = J + (pf*F + pV*Uk(2) - pB*Bo - pD*Do) * delta_time; 
%         J = J + (Uk - u_opt)' * (Uk - u_opt) * delta_time;
%         J = J + (Xk - xdot_val_rf_ss)' * (Xk - xdot_val_rf_ss) * delta_time;
        
%         %params.weights(count:count+nx-1,1)
%         J = J + (params.weights(count:count+nx-1,1).*(Xk - xdot_val_rf_ss))' *(Xk - xdot_val_rf_ss)*delta_time;
%         count = count + nx;
%         %params.weights(count:count+nu-1,1)
%         J = J + (Uk - u_opt)' * (params.weights(count:count+nu-1,1).*(Uk - u_opt))*delta_time;
%         count = count + nu;

        J  = J + (Qmax(1:nx,1).*(Xkj{j} - xdot_val_rf_ss))' * (Xkj{j} - xdot_val_rf_ss) * delta_time;

        % add state soft-constraint
        %rho = 1e3;
        %J = J + rho*epsilon;

    end
end

