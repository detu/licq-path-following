function [prob] = distColACstr_pnNew(p)
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
    nDual   = numel(y.lam_g);
    nParam  = numel(p);
    % model parameters
    NT = 41;
    Uf = 0.3;                                 % Feeding rate F
    [~,state,xdot,inputs] = DistColACstr(Uf);
    sf = Function('sf',{state,inputs}, {xdot});
    
    % dimensions
    global nx nu nk d tf ns;
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
    load CstrDistXinit.mat;
    xf             = Xinit(1:84);
    u_opt          = Xinit(85:89);
    
    % prices
    pf = 1; 
    pV = 0.02;
    pB = 2; 
    pD = 0;
    
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
    X0   = MX.sym('X0', nx);
    V    = {V{:}, X0};
    cons = {cons{:}, X0 - x(1:nx,1)};
    cons_x0 = X0 - x(1:nx,1);
        
    % formulate the NLP
    Xk = X0;


    load Qmax.mat;
    params.Qmax = Qmax;

    ssoftc = 0;
    for i=1:N
        [obj,cons,V,Xk,params,ssoftc] = iterateOnPredictionHorizon(Xk, V, cons, obj, params, i,ssoftc);
    end
    
    V = vertcat(V{:});
    % Concatenate constraints
    cons  = vertcat(cons{:});
    
    % objective function and constraint functions
    f   = Function('f', {V}, {obj}, char('V'), char('objective'));
    c   = Function('c', {V}, {cons}, char('V'), char('constraint'));
    cx0 = Function('cx0', {X0}, {cons_x0}, char('X0'), char('constraint'));
    
    % construct Lagrangian
    lag_expr = obj + y.lam_g'*cons;

   
    [Hj,g] = hessian(obj,V);
    Hobj   = Function('Hj',{V},{Hj});
    gj     = Function('gj',{V},{g});
    Hl     = hessian(lag_expr,V);
    H      = Function('Hl',{V},{Hl});
    J      = Function('J',{V},{jacobian(cons,V),cons});
    Jp     = Function('Jp',{X0},{jacobian(cons_x0,X0),cons_x0});
    
    f                        = f(x);
    g                        = gj(x);
    H                        = H(x);
    Lxp                      = H(1:nPrimal,1:nParam);
    J                        = J(x);
    Jtemp                    = zeros(nDual,nParam);
    cp                       = Jp(x(1:nParam));
    Jtemp(1:nParam,1:nParam) = full(cp);
    cp                       = sparse(Jtemp);
    cst                      = c(x);
    
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
   
   global nx nu nk d ns N;

   count = 1;
   for k=0:nk-1
        % New NLP variable for the control
        Uk  = MX.sym(['U_' num2str((iter-1)*nk+k)], nu);
        V   = {V{:}, Uk};
        
        Jcontrol   = (Qmax(nx+1:nx+nu,1).*(Uk - u_opt))' * (Uk - u_opt);

        % State at collocation points
        Xkj = {};

        SumX1 = 0;
        for j=1:d
            Xkj{j} = MX.sym(['X_' num2str((iter-1)*nk+k) '_' num2str(j)], nx);
            V      = {V{:}, Xkj{j}};

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
           fj   = sf(Xkj{j},Uk);
           cons = {cons{:}, h*fj - xp};

           % Add contribution to the end state
           Xk_end = Xk_end + D(j+1)*Xkj{j};

        end    

        % New NLP variable for state at end of interval
        Xk  = MX.sym(['X_' num2str((iter-1)*nk+k)], nx);
        V   = {V{:}, Xk};

        % Add equality constraint
        cons= {cons{:}, Xk_end-Xk};


        Jecon  = (pf*F_0 + pV*Uk(2) - pB*Uk(5) - pD*Uk(4)) * delta_time;
        Jstate =(Qmax(1:nx,1).*(Xk - xdot_val_rf_ss))' * (Xk - xdot_val_rf_ss) * delta_time;
        
        % compute rotated cost function
        fm  = sf(Xk,Uk);
        % load Lagrange multipliers from steady-state optimization
        load LamdaCstrDist.mat; % lamda
        Jmodel = lamda'*fm;
        
        alpha  = 1;
        beta   = 1;
        gamma  = 1;
        rho    = 0.5;
        obj    = obj + alpha*Jcontrol + gamma*Jstate + beta*Jecon;


    end
end

