function [primal,dual] = optDistillationParam(param)
%OPTDISTILLATIONPARAM Summary of this function goes here
% 
% [OUTPUTARGS] = OPTDISTILLATIONPARAM(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2016/02/02 20:24:39 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2016

import casadi.* 
format long;

NT = 41;
Uf = 1.2;           % Feeding rate F
% invoke the model
%[t,state,xdot,inputs] = DistColA(Uf);
[~,state,xdot,inputs] = DistColAP(Uf,param);
% [t,state,xdot,inputs] = DistColAPCM(Uf,param);

%f = SXFunction('f', [t,state,inputs],[xdot]);
%f = SXFunction('f',daeIn('t',t, 'x',state, 'p', inputs), daeOut('ode', xdot));
%f = SXFunction('f',{t,state,inputs}, {xdot});
f = Function('f',{state,inputs}, {xdot});
% load x_init
% x_init = x_init(1:41);
load u_opt_ss.mat;
load cola_init.mat;
%x_init = Xinit;
%u_init = [2.70629; 3.20629];
%u_init = [2.827; 3.454];
%u_init = [2.784262975002144; 3.411499338095353];  % these are obtained from FMINCON steady-state optimization.

% bound constraints
%xB_min = 0.0100;
%xB_min = 0.0080;
xB_max = 0.008;
%xD_min = 0.9500;
xD_min = 0.9490;
V_max  = 4.008 ;

% State bounds and initial guess
x_min =  zeros(82,1);
x_max =  ones(82,1);
% x_min =  zeros(41,1);
% x_max =  ones(41,1);
%x_min(1)  = xB_min;
x_max(1)  = xB_max;
x_min(41) = xD_min;

% Control bounds
u_min = [0.1; 0.1];
u_max = [10; V_max];


%------------------------------------------------------------------------------------------------------#    
% Dynamic Optimization 
% Collocation method
%------------------------------------------------------------------------------------------------------#
    
% % dimensions
% nx = 82;
% % nx = 41;
% nu = 2;
%     
% nk = 10;       % control discretization
% tf = 10.0;     % end time
% 
% % choose collocation points 
% tau_root = collocationPoints(3,'radau');
% 
% % degree of interpolationg polynomial 
% %d = len(tau_root) - 1;
% d = size(tau_root,2)-1;
% 
% % size of the finite elements 
% h = tf/nk; 
% 
% % coefficients of the collocation equation 
% C = zeros(d+1,d+1);
% 
% % coefficients of the continuity equation 
% D = zeros(d+1);
% 
% % dimensionless time inside one control interval 
% tau = SX.sym('tau');
% 
% % all collocation time points 
% T = zeros(nk, d+1);        
% for k=1:nk
%     for j=1:d+1
%         T(k,j) = h*(k + tau_root{j});
%     end
% end
% 
% % for all collocation points         
% for j=1:d+1
%     % construct Lagrange polynomials to get the polynomials basis at the collocation point
%     L =1;
%     for r =1:d+1
%         if r ~= j
%             L = L * ((tau - tau_root{r}) / (tau_root{j} - tau_root{r}));
%         end
%     end
%     lfcn = SXFunction('lfcn',{tau},{L});
%     
%     % evaluate the polynomial at the final time to get the coefficients of the continuity equation 
%     %D(j) = lfcn(1.0);
%     lfcn.setInput(1.0);
%     lfcn.evaluate();
%     D(j) = full(lfcn.getOutput());
%     
%     % evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
%     tfcn = lfcn.tangent();
%     for r=1:d+1
%         %[C(j,r),~] = tfcn(tau_root(r));
%         tfcn.setInput(tau_root{r});
%         tfcn.evaluate();
%         C(j,r) = full(tfcn.getOutput());
%     end
%     
% end
%         
% % total number of variables 
% NX  = nk*(d+1)*nx;            % collocated states 
% NU  = nk*nu;                  % parameterized controls 
% NXF = nx;                     % final state 
% NV  = NX+NU+NXF;
% 
% % NLP variable vector 
% V = MX.sym('V',NV,1);
% %V = cd.SX.sym('V',NV);
% %V = SX.sym('V',NV);
% 
% % all variables with bounds and initial guess
% vars_lb   = zeros(NV,1);
% vars_ub   = inf*ones(NV,1);
% vars_init = zeros(NV,1);
% offset    = 1;
% 
% % get collocated states and parametrized control 
% %X = MX.sym('X',((nk+1)*(d+1)*nx),1);
% %X = MX.sym('X',(nk+1),(d+1));
% %U = MX.sym('U',NU,1);
% X{(nk+1),(d+1)} = [];
% U{nk} = [];
% for k=1:nk
%     % collocated states
%     for j=1:d+1
%         % get the expression for the state vector
%         %X(offset:offset+nx-1,1) =  V(offset:offset+nx-1,1);
%         X{k,j} = V(offset:offset+nx-1,1);
%         
%         % add the initial condition
%         vars_init(offset:offset+nx-1,1) = x_init;
%         
%         % add bounds
%         vars_lb(offset:offset+nx-1) = x_min;
%         vars_ub(offset:offset+nx-1) = x_max;
% 
%         offset = offset + nx;
%     end
%     
%     % Parametrized controls
%     %U(((k-1)*nu)+1:(k*nu)) = V(offset:offset+nu-1);
%     U{k} = V(offset:offset+nu-1);
%     vars_lb(offset:offset+nu-1)   = u_min;
%     vars_ub(offset:offset+nu-1)   = u_max;
%     vars_init(offset:offset+nu-1) = u_init;
%     offset = offset + nu;
% end
% % State at end time
% %X(offset:offset+nx-1) = V(offset:offset+nx-1);
% X{nk+1,1} = V(offset:offset+nx-1);
% vars_lb(offset:offset+nx-1)   = x_min;
% vars_ub(offset:offset+nx-1)   = x_max;
% vars_init(offset:offset+nx-1) = x_init;
% offset = offset + nx;
% 
% % Constraint function for the NLP
% g   = [];
% lbg = [];
% ubg = [];
% p_F = 1.0;
% p_V = 0.01;
% p_B = 1.0;
% p_D = 2.0; 
% alpha = 1;
% beta  = 0.1;
% gamma = 100;
% 
% % load optimized control input and state values from steady-state
% % optimization
% %load sim_steady_state.mat;
% %load u_opt_ss.mat;
% delta_time = 60; % [minute] convert second to minute
% 
% F   = Uf;
% % Bo  = 0.5;  % Fix value from the model
% % Do  = 0.5;  % Fix value from the model
% %Do = 0.627
% %Bo = 0.573
% KcB=10;  KcD=10;         % controller gains
% MDs=0.5; MBs=0.5;        % Nominal holdups - these are rather small  
% Ds=0.5; Bs=0.5;          % Nominal flows
% obj = 0;
% 
% % For all finite elements
% for k=1:nk
%     % For all collocation points
%     for j=2:d+1
%         % Get an expression for the state derivative at the collocation point
%         xp_jk = 0;
%         for r=1:d+1
%             xp_jk = xp_jk + C(r,j)*X{k,r};
%         end
%         % Add collocation equations to the NLP
%         fk  = f.call({T(k,j), X{k,j}, U{k}});
%         g   = [g;(h*fk{1} - xp_jk)];
%         lbg = [lbg;zeros(nx,1)];
%         ubg = [ubg;zeros(nx,1)];
%     end
%     
%     % Get an expression for the state at the end of the finite element
%     xf_k = 0;
%     for r=1:(d+1)
%         xf_k = xf_k + D(r)*X{k,r};
%     end
%     
%     % add continuity equation for NLP 
%     g   = [g;X{k+1,1} - xf_k];
%     lbg = [lbg;zeros(nx,1)];
%     ubg = [ubg;zeros(nx,1)];
%     
%     % objective function
%     MB = xf_k(NT+1);  
%     MD = xf_k(2*NT);            % Actual reboiler and condenser holdup
%     Do = Ds+(MD-MDs)*KcD;       % Distillate flow
%     Bo = Bs+(MB-MBs)*KcB;       % Bottoms flow
%     econ_term    = (p_F*F + p_V*U{k}(2) - p_B*Bo - p_D*Do)*delta_time;       % economic term
%     control_term = (U{k}(1) - u_opt(1)) + (U{k}(2) - u_opt(2))*delta_time;   % tracking: control input term
%     state_term   = (xf_k - xf)'*(xf_k - xf)*delta_time;                      % tracking: states term
%     obj          = obj + alpha*econ_term + beta*control_term + gamma*state_term;            
% end
%     
% %obj = p_F*F + p_V*U{k}(2) - p_B*Bo - p_D*Do;    
% % Concatenate constraints
% g = vertcat(g);
% fprintf('number of variables: %d \n', NV);
% fprintf('number of constraints:%d \n',size(g,1));
% 
% % NLP
% nlp = MXFunction('nlp', nlpIn('x',V),nlpOut('f',obj,'g',g));
% 
% 
% %% SOLVE THE NLP
% 
% % Set options
% opts = struct;
% %opts.expand = True;
% %opts["max_iter"] = 500
% %opts["linear_solver"] = 'ma27'
% opts.hessian_approximation = 'exact';
% opts.print_level = 0;
% 
% % Allocate an NLP solver
% solver = NlpSolver('solver', 'ipopt', nlp, opts);
% arg    = struct;
% 
% % Initial condition
% arg.x0 = vars_init;
% 
% % Bounds on x
% arg.lbx = vars_lb;
% arg.ubx = vars_ub;
% 
% % Bounds on g
% arg.lbg = lbg;
% arg.ubg = ubg;
% 
% 
% % Solve the problem
% startnlp = tic;
% res = solver(arg);
% elapsednlp = toc(startnlp);
% fprintf('NLP solver runtime:%f\n', elapsednlp);

% test steady-state optimization
%disA = Function('disA', {state,inputs}, {xdot}, char('x', 'u'), char('xdot'));
%rf   = rootfinder('rf','newton', disA);

% skip the steady-state optimization using IPOPT because needs to be
% supplied with FEASIBLE initial guess !
%u_opt          = u_init;                           % from FMINCON 
%xdot_val_rf_ss = full(rf(x_init, u_opt));

xdot_val_rf_ss = xf;  % from FMINCON 

% prices
pf = 1; 
pV = 0.01; 
pB = 1; 
pD = 2;
% 
F = Uf;
% %V = u{2};
% V = inputs{2};
% 
% controller gains
KcB = 10;  
KcD = 10;
% Nominal holdups - these are rather small 
MDs = 0.5; 
MBs = 0.5;         
% Nominal flows
Ds  = 0.5; 
Bs  = 0.5;

%% Dynamic Optimization
% create a function !   
% Collocation method
%------------------------------------------------------------------------------------------------------#
    
% dimensions
nx = 82;
nu = 2;
    
nk = 10;       % control discretization
tf = 10.0;     % end time
%nk = 100;
%tf = 50.0;
h  = tf/nk;

% preparing collocation matrices
[~,C,D,d] = collocationSetup();

% start with an empty NLP
w   = {};      % decision variables contain both control and state variables
w0  = [];      % initial guess
lbw = [];      % lower bound for decision variable
ubw = [];      % upper bound
J   = 0;       % objective function
g   = {};      % nonlinear constraint
lbg = [];      % lower bound for nonlinear constraint
ubg = [];      % upper bound

%delta_time = 60; % [minute] convert second to minute
delta_time = 1;
alpha = 1;
beta  = 1;
gamma = 1;

% "Lift" initial conditions
X0  = MX.sym('X0', nx);
w   = {w{:}, X0};
lbw = [lbw; x_min];
ubw = [ubw; x_max];
w0  = [w0; xdot_val_rf_ss];   % use initial guess from steady-state optimization results. this should give zero objective function (for tracking term!)
g   = {g{:}, X0 - xdot_val_rf_ss};
lbg = [lbg; zeros(nx,1)];
ubg = [ubg; zeros(nx,1)];

% formulate the NLP
Xk = X0;
for k=0:nk-1
    % New NLP variable for the control
    Uk  = MX.sym(['U_' num2str(k)], nu);
    w   = {w{:}, Uk};
    lbw = [lbw; u_min];
    ubw = [ubw; u_max];
    w0  = [w0;  u_opt];       % optimized results from steady-state optimization

    % State at collocation points
    Xkj = {};
    for j=1:d
        Xkj{j} = MX.sym(['X_' num2str(k) '_' num2str(j)], nx);
        w      = {w{:}, Xkj{j}};
        lbw    = [lbw; x_min];
        ubw    = [ubw; x_max];
        w0     = [w0; xdot_val_rf_ss];
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
       fj  = f(Xkj{j},Uk);
       g   = {g{:}, h*fj - xp};
       lbg = [lbg; zeros(nx,1)];
       ubg = [ubg; zeros(nx,1)];
       
       % Add contribution to the end state
       Xk_end = Xk_end + D(j+1)*Xkj{j};
  
       % Add contribution to quadrature function
       %J = J + B(j+1)*qj*h;
    end    
   
    % New NLP variable for state at end of interval
    Xk  = MX.sym(['X_' num2str(k+1)], nx);
    w   = {w{:}, Xk};
    lbw = [lbw; x_min];
    ubw = [ubw; x_max];
    w0  = [w0; xdot_val_rf_ss];

    % Add equality constraint
    g   = {g{:}, Xk_end-Xk};
    lbg = [lbg; zeros(nx,1)];
    ubg = [ubg; zeros(nx,1)];
    
    % objective function
    MB = Xk(NT+1);  
    MD = Xk(2*NT);            % Actual reboiler and condenser holdup
    Do = Ds+(MD-MDs)*KcD;       % Distillate flow
    Bo = Bs+(MB-MBs)*KcB;       % Bottoms flow
    econ_term    = (pf*F + pV*Uk(2) - pB*Bo - pD*Do)*delta_time;            % economic term
    control_term = (Uk - u_opt)'*(Uk - u_opt)*delta_time;                   % tracking: control input term
    state_term   = (Xk - xdot_val_rf_ss)'*(Xk - xdot_val_rf_ss)*delta_time;  % tracking: states term
    J            = J + alpha*econ_term + beta*control_term + gamma*state_term;
    %J            = J + state_term;
    %J            = J + control_term;
end

% Create an NLP solver
prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));
solver = nlpsol('solver', 'ipopt', prob);

% Solve the NLP
sol   = solver('x0', w0, 'lbx', lbw, 'ubx', ubw, 'lbg', lbg, 'ubg', ubg);
%sol   = solver('x0', w0, 'lbg', lbg, 'ubg', ubg);  % without bound constraint... drive system immediately to steady-state equilibrium! (with combination objective function, IPOPT doesn't converge quickly!)

% Print the optimal cost
fprintf('optimal cost: %f \n', full(sol.f));

% Retrieve the solution
primal = full(sol.x);
dual   = full(sol.lam_g);
% dual   = abs(full(sol.lam_g));

end
