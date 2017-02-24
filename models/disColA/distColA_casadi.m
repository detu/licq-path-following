function [prob] = distColA_casadi(p)
%DISTCOLA_CASADI Summary of this function goes here
% 
% [OUTPUTARGS] = DISTCOLA_CASADI(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2016/02/01 17:45:24 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2016

import casadi.*

prob = struct('neq',{0},'niq',{0},'cin',{0},'ceq',{0},'dp_in',{0},'dp_eq',{0},'hess',{0},'lxp',{0},'x',0,'name',0);
prob.neq  = 2000;         % HARD-CODE !    % number of equality constraint
prob.niq  = 0;            % number of inequality constraint
prob.name = 'Distillation Column A model';
prob.x    = zeros(2,1);

prob.obj  = (@(x,y,p)(objective(x,y,p)));

end

function [f,g,H,Lxp,cst,J,cp,Jeq,dpe] = objective(x,y,p)
    import casadi.*
   
    % model parameters
    NT = 41;
    Uf = 1.2;              % Feeding rate F
    P  = SX.sym('P');      % parameter used in path-following algorithm
    % invoke the model
    %[~,state,xdot,inputs] = DistColAP(Uf,P);
    [~,state,xdot,inputs] = DistColAP(Uf,p);
    sf = Function('sf',{state,inputs}, {xdot});
    
    % dimensions
    nx = 82;
    nu = 2;

    nk = 10;       % control discretization
    tf = 10.0;     % end time
    h  = tf/nk;
    
    % preparing collocation matrices
    [~,C,D,d] = collocationSetup();
    
    % start with an empty NLP
%     % total number of variables 
%     NX  = nk*(d+1)*nx;            % collocated states 
%     NU  = nk*nu;                  % parameterized controls 
%     NXF = nx;                     % final state 
%     NV  = NX+NU+NXF;

    % NLP variable vector 
    %V = SX.sym('V',NV);
    V    = {};      % decision variables contain both control and state variables
    obj  = 0;       % objective function
    cons = {};      % nonlinear constraint
    
    delta_time = 1;
    alpha = 1;
    beta  = 1;
    gamma = 1;

    % Initial states and controls 
%     load cola_init.mat;
%     x_init = Xinit;
%     u_init = [2.784262975002144; 3.411499338095353];  % these are obtained from FMINCON steady-state optimization.
%     % test steady-state optimization
%     disA = Function('disA', {state,inputs}, {xdot}, char('x', 'u'), char('xdot'));
%     rf   = rootfinder('rf','newton', disA);
% 
%     % skip the steady-state optimization using IPOPT because needs to be
%     % supplied with FEASIBLE initial guess !
%     u_opt          = u_init;                           % from FMINCON 
%     xdot_val_rf_ss = full(rf(x_init, u_opt));
    load u_opt_ss.mat;
    xdot_val_rf_ss = xf;
    
    % prices
    pf = 1; 
    pV = 0.01; 
    pB = 1; 
    pD = 2;
    F  = Uf;
    % controller gains
    KcB = 10;  
    KcD = 10;
    % Nominal holdups - these are rather small 
    MDs = 0.5; 
    MBs = 0.5;         
    % Nominal flows
    Ds  = 0.5; 
    Bs  = 0.5;
    
    % "Lift" initial conditions
    X0   = SX.sym('X0', nx);
    V    = {V{:}, X0};
    cons = {cons{:}, X0 - xdot_val_rf_ss};
    
    % formulate the NLP
    Xk = X0;
    for k=0:nk-1
        % New NLP variable for the control
        Uk  = SX.sym(['U_' num2str(k)], nu);
        V   = {V{:}, Uk};
        
        % State at collocation points
        Xkj = {};
        for j=1:d
            Xkj{j} = SX.sym(['X_' num2str(k) '_' num2str(j)], nx);
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
            %[fj, qj] = easycall(f, Xkj{j},Uk);
            fj   = sf(Xkj{j},Uk);
            cons = {cons{:}, h*fj - xp};
            
            % Add contribution to the end state
            Xk_end = Xk_end + D(j+1)*Xkj{j};
            
            % Add contribution to quadrature function
            %J = J + B(j+1)*qj*h;
        end
        
        % New NLP variable for state at end of interval
        Xk   = SX.sym(['X_' num2str(k+1)], nx);
        V    = {V{:}, Xk};
        
        % Add equality constraint
        cons = {cons{:}, Xk_end-Xk};
        
        % objective function
        MB = Xk(NT+1);
        MD = Xk(2*NT);            % Actual reboiler and condenser holdup
        Do = Ds+(MD-MDs)*KcD;       % Distillate flow
        Bo = Bs+(MB-MBs)*KcB;       % Bottoms flow
        econ_term    = (pf*F + pV*Uk(2) - pB*Bo - pD*Do)*delta_time;            % economic term
        control_term = (Uk - u_opt)'*(Uk - u_opt)*delta_time;                   % tracking: control input term
        state_term   = (Xk - xdot_val_rf_ss)'*(Xk - xdot_val_rf_ss)*delta_time;  % tracking: states term
        obj          = obj + alpha*econ_term + beta*control_term + gamma*state_term;
        %obj            = obj + state_term;
        %obj            = obj + control_term;
    end

    
%     % choose collocation points 
%     tau_root = collocationPoints(3,'radau');
% 
%     % degree of interpolationg polynomial 
%     %d = len(tau_root) - 1;
%     d = size(tau_root,2)-1;
% 
%     % size of the finite elements 
%     h = tf/nk; 
% 
%     % coefficients of the collocation equation 
%     C = zeros(d+1,d+1);
% 
%     % coefficients of the continuity equation 
%     D = zeros(d+1);
% 
%     % dimensionless time inside one control interval 
%     tau = SX.sym('tau');
% 
%     % all collocation time points 
%     T = zeros(nk, d+1);        
%     for k=1:nk
%         for j=1:d+1
%             T(k,j) = h*(k + tau_root{j});
%         end
%     end
%     
%     % for all collocation points         
%     for j=1:d+1
%         % construct Lagrange polynomials to get the polynomials basis at the collocation point
%         L =1;
%         for r =1:d+1
%             if r ~= j
%                 L = L * ((tau - tau_root{r}) / (tau_root{j} - tau_root{r}));
%             end
%         end
%         lfcn = SXFunction('lfcn',{tau},{L});
% 
%         % evaluate the polynomial at the final time to get the coefficients of the continuity equation 
%         %D(j) = lfcn(1.0);
%         lfcn.setInput(1.0);
%         lfcn.evaluate();
%         D(j) = full(lfcn.getOutput());
% 
%         % evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
%         tfcn = lfcn.tangent();
%         for r=1:d+1
%             %[C(j,r),~] = tfcn(tau_root(r));
%             tfcn.setInput(tau_root{r});
%             tfcn.evaluate();
%             C(j,r) = full(tfcn.getOutput());
%         end
% 
%     end
% 
%     % total number of variables 
%     NX  = nk*(d+1)*nx;            % collocated states 
%     NU  = nk*nu;                  % parameterized controls 
%     NXF = nx;                     % final state 
%     NV  = NX+NU+NXF;
% 
%     % NLP variable vector 
%     %V = MX.sym('V',NV,1);
%     %V = cd.SX.sym('V',NV);
%     V = SX.sym('V',NV);
% 
%     % all variables with bounds and initial guess
%     %vars_lb   = zeros(NV,1);
%     %vars_ub   = inf*ones(NV,1);
%     %vars_init = zeros(NV,1);
%     offset    = 1;
% 
%     % get collocated states and parametrized control 
%     %X = MX.sym('X',((nk+1)*(d+1)*nx),1);
%     %X = MX.sym('X',(nk+1),(d+1));
%     %U = MX.sym('U',NU,1);
%     X{(nk+1),(d+1)} = [];
%     U{nk} = [];
%     for k=1:nk
%         % collocated states
%         for j=1:d+1
%             % get the expression for the state vector
%             %X(offset:offset+nx-1,1) =  V(offset:offset+nx-1,1);
%             X{k,j} = V(offset:offset+nx-1,1);
% 
%             % add the initial condition
%             %vars_init(offset:offset+nx-1,1) = x_init;
% 
%             % add bounds
%             %vars_lb(offset:offset+nx-1) = x_min;
%             %vars_ub(offset:offset+nx-1) = x_max;
% 
%             offset = offset + nx;
%         end
% 
%         % Parametrized controls
%         %U(((k-1)*nu)+1:(k*nu)) = V(offset:offset+nu-1);
%         U{k} = V(offset:offset+nu-1);
%         %vars_lb(offset:offset+nu-1)   = u_min;
%         %vars_ub(offset:offset+nu-1)   = u_max;
%         %vars_init(offset:offset+nu-1) = u_init;
%         offset = offset + nu;
%     end
%     % State at end time
%     %X(offset:offset+nx-1) = V(offset:offset+nx-1);
%     X{nk+1,1} = V(offset:offset+nx-1);
%     %vars_lb(offset:offset+nx-1)   = x_min;
%     %vars_ub(offset:offset+nx-1)   = x_max;
%     %vars_init(offset:offset+nx-1) = x_init;
%     offset = offset + nx;
% 
%     % Constraint function for the NLP
%     cons = [];
%     %lbg  = [];
%     %ubg  = [];
%     p_F  = 1.0;
%     p_V  = 0.01;
%     p_B  = 1.0;
%     p_D  = 2.0; 
%     alpha = 1;
%     beta  = 0.1;
%     gamma = 100;
%     
%     % load optimized control input and state values from steady-state
%     % optimization
%     %load sim_steady_state.mat;
%     load u_opt_ss.mat;
%     %xf = xf(1:41);
%     delta_time = 60; % [minute] convert second to minute
% 
%     F   = Uf;
% %     Bo  = 0.5;  % Fix value from the model
% %     Do  = 0.5;  % Fix value from the model
%     %Do = 0.627
%     %Bo = 0.573
%     KcB=10;  KcD=10;         % controller gains
%     MDs=0.5; MBs=0.5;        % Nominal holdups - these are rather small  
%     Ds=0.5; Bs=0.5;          % Nominal flows
%     obj = 0;
% 
%     % For all finite elements
%     for k=1:nk
%         % For all collocation points
%         for j=2:d+1
%             % Get an expression for the state derivative at the collocation point
%             xp_jk = 0;
%             for r=1:d+1
%                 xp_jk = xp_jk + C(r,j)*X{k,r};
%             end
%             % Add collocation equations to the NLP
%             fk   = sf.call({T(k,j), X{k,j}, U{k}});
%             cons = [cons;(h*fk{1} - xp_jk)];
%             %lbg  = [lbg;zeros(nx,1)];
%             %ubg  = [ubg;zeros(nx,1)];
%         end
% 
%         % Get an expression for the state at the end of the finite element
%         xf_k = 0;
%         for r=1:(d+1)
%             xf_k = xf_k + D(r)*X{k,r};
%         end
% 
%         % add continuity equation for NLP 
%         cons = [cons;X{k+1,1} - xf_k];
%         %lbg  = [lbg;zeros(nx,1)];
%         %ubg  = [ubg;zeros(nx,1)];
%         
%         % objective function
%         MB = xf_k(NT+1);
%         MD = xf_k(2*NT);            % Actual reboiler and condenser holdup
%         Do = Ds+(MD-MDs)*KcD;       % Distillate flow
%         Bo = Bs+(MB-MBs)*KcB;       % Bottoms flow
%         econ_term    = (p_F*F + p_V*U{k}(2) - p_B*Bo - p_D*Do)*delta_time;       % economic term
%         control_term = (U{k}(1) - u_opt(1)) + (U{k}(2) - u_opt(2))*delta_time;   % tracking: control input term
%         state_term   = (xf_k - xf)'*(xf_k - xf)*delta_time;                      % tracking: states term
%         obj          = obj + alpha*econ_term + beta*control_term + gamma*state_term;
%     end

    %obj = p_F*F + p_V*U{k}(2) - p_B*Bo - p_D*Do;    
    % Concatenate constraints
    %cons = vertcat(cons);
    cons  = vertcat(cons{:});
    
    % construct Lagrangian
    lag_expr = obj + y'*cons;
    
    % objective function related derivatives
    %opf = struct;
    %opf.input_scheme = char('V','P');
    %opf.output_scheme = char('obj');
    %f = SXFunction('f',{V,P},{obj},opf);
    %g = SXFunction(f.gradient);
    V = vertcat(V{:});
    f = Function('f', {V,P}, {obj}, char('V', 'P'), char('obj'));
    g = f.gradient();
    
    % Lagrangian and constraint related derivatives
    %op       = struct;
    %op.input_scheme  = char('V','P');
    %op.output_scheme = char('lag_expr');
    %lagr       = SXFunction('lagr',{V,P},{lag_expr},op);
    %H          = SXFunction(lagr.hessian('V','lag_expr'));
    %jac_lagr_x = SXFunction(lagr.jacobian('V','lag_expr'));
    %Lxp        = jac_lagr_x.jacobian(1,0);
    lagr       = Function('lagr', {V,P}, {lag_expr}, char('V','P'), char('lag_expr'));
    H          = Function(lagr.hessian('V','lag_expr'));
    jac_lagr_x = Function(lagr.jacobian('V','lag_expr'));
    Lxp        = jac_lagr_x.jacobian(1,0);
    
    % Inequality constraint
    %opc = struct;
    %opc.input_scheme  = char('V','P');
    %opc.output_scheme = char('c_expr');
    %cst = SXFunction('cst',{V,P},{cons},opc);
    %J   = cst.jacobian('V','c_expr');
    %cp  = cst.jacobian('P','c_expr');
    cst = Function('cst', {V,P}, {cons}, char('V','P'), char('c_expr'));
    J   = cst.jacobian('V','c_expr');
    cp  = cst.jacobian('P','c_expr');
 
    % evaluate and obtain their values
%     f.setInput(x,'V');
%     f.setInput(p,'P');
%     g.setInput(x,'V');
%     g.setInput(p,'P');
%     H.setInput(x,'V');
%     H.setInput(p,'P');
%     Lxp.setInput(x,'V');
%     Lxp.setInput(p,'P');
%     J.setInput(x,'V');
%     J.setInput(p,'P');
%     cp.setInput(x,'V');
%     cp.setInput(p,'P');
%     cst.setInput(x,'V');
%     cst.setInput(p,'P');
% 
%     f.evaluate();
%     g.evaluate();
%     H.evaluate();   
%     J.evaluate();
%     cp.evaluate();
%     Lxp.evaluate();
%     cst.evaluate();
    f   = f(x,p);
    g   = g(x,p);
    H   = H(x,p);
    Lxp = Lxp(x,p);
    J   = J(x,p);
    cp  = cp(x,p);
    cst = cst(x,p);
    
%     f   = full(f.getOutput());
%     g   = sparse(g.getOutput());
%     H   = sparse(H.getOutput());
%     Lxp = sparse(Lxp.getOutput());
%     J   = sparse(J.getOutput());
%     cp  = sparse(cp.getOutput());
%     cst = full(cst.getOutput());

    f   = full(f);
    g   = sparse(g);    
    H   = sparse(H);
    Lxp = sparse(Lxp);
    J   = sparse(J);
    cp  = sparse(cp);
    cst = full(cst);

    % Equality constraint
    %Jeq = [];
    %dpe = [];
    Jeq = J;
    dpe = cp;
end
