function [prob] = distColA_casadi_sens(p)
%DISTCOLA_CASADI_SENS Summary of this function goes here
% 
% [OUTPUTARGS] = DISTCOLA_CASADI_SENS(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2016/04/04 12:05:14 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2016

import casadi.*

prob = struct('neq',{0},'niq',{0},'cin',{0},'ceq',{0},'dp_in',{0},'dp_eq',{0},'hess',{0},'lxp',{0},'x',0,'name',0);
prob.neq  = 2000;         % HARD-CODE !    % number of equality constraint
prob.niq  = 0;            % number of inequality constraint
prob.name = 'Distillation Column A model';
prob.x    = zeros(2,1);

%prob.obj  = (@(x,y,p,b)(objective(x,y,p,b)));
prob.obj  = (@(x,p,b)(objective(x,p,b)));

end

%function [f,g,H,Lxp,cst,J,cp,Jeq,dpe] = objective(x,y,p,b)
%function [fxx,fxw,fwx,fww,cx,cw] = objective(x,p,b)
function [fxx,cx] = objective(x,p,b)

    import casadi.*
   
    % change parameter from Uf to X0 (initial state)
    
    % model parameters
    NT = 41;
    Uf = 1.2;              % Feeding rate F
    %P  = SX.sym('P');      % parameter used in path-following algorithm
    % invoke the model
    %[~,state,xdot,inputs] = DistColAP(Uf,p);
    [~,state,xdot,inputs] = DistColA(Uf);
    sf = Function('sf',{state,inputs}, {xdot});
    
    % dimensions
    nx = 82;
    nu = 2;

    nk = 10;       % control discretization
    tf = 10.0;     % end time
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
    %P    = SX.sym('P', nx);   % initial state as parameter
    %P    = X0;
    %cons = {cons{:}, P - xdot_val_rf_ss};
    
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

    % define objective function with bound constraints
    V = vertcat(V{:});
    mu   = 1e-8;  % assumption. this should be obtained from IPOPT
    %mObj = obj + sum(log(V - b.lb)) + sum(log(b.ub - V));
    mObj = obj;
%     for i=1:numel(V)
%         mObj = mObj - mu*log(V(i) - b.lb(i)) - mu*log(b.ub(i) - V(i)); 
%     end

    % Concatenate constraints
    cons  = vertcat(cons{:});
    
    % objective function and constraint functions
    %f = Function('f', {V,P}, {mObj}, char('V', 'P'), char('objective'));
    %c = Function('c', {V,P}, {cons}, char('V', 'P'), char('constraint'));
    %f = Function('f', {V,X0}, {mObj}, char('V', 'X0'), char('objective'));
    %c = Function('c', {V,X0}, {cons}, char('V', 'X0'), char('constraint'));
    %f = Function('f', {V}, {mObj}, char('V'), char('objective'));
    %c = Function('c', {V}, {cons}, char('V'), char('constraint'));
    
    f = Function('f', {V,X0}, {mObj}, char('V', 'P'), char('objective'));
    c = Function('c', {V,X0}, {cons}, char('V', 'P'), char('constraint'));

    % take derivatives of f and c 
    %fx  = Function('fx' ,{V},{jacobian(f,V)});
    %fx  = f.jacobian(0,0);
    %fw  = Function('fw' ,{P},{jacobian(f,P)});
    %fw  = f.jacobian(1,0);
    %fxx = Function('fxx',{V},{hessian(f,V)});
    fxx = f.hessian(0,0);
    %fww = Function('fww',{P},{hessian(f,P)});
    %fww = f.hessian(1,0);
    %fxw = Function('fxw',{P},{jacobian(fx,P)});
    %fxw = fx.jacobian(1,0);
    %fwx = Function('fwx',{V},{jacobian(fw,V)});
    %fwx = fw.jacobian(0,0);
    %cx  = Function('cx' ,{V},{jacobian(c,V)});
    cx = c.jacobian(0,0);
    %cw  = Function('cw' ,{P},{jacobian(c,P)});
    %cw = c.jacobian(1,0);
    
    % evaluate and obtain their values
    %fxw   = full(fxw(x,p));
    %fww   = full(fww(x,p));
    %fwx   = full(fwx(x,p));
    %fxx   = full(fxx(x,p));
    %cx    = full(cx(x,p));
    %cw    = full(cw(x,p));
    fxx   = full(fxx(x));
    cx    = full(cx(x));
end

