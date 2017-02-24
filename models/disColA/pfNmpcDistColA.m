function pfNmpcDistColA
%PFNMPCDISTCOLA Summary of this function goes here
% 
% [OUTPUTARGS] = PFNMPCDISTCOLA(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2016/04/06 21:22:25 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2016

import casadi.* 
format long;

% IMPLEMENT PATH-FOLLOWING ALGORITHM !
% CHANGE THE PARAMETER TO INITIAL STATE ! NOT X0 but X1 !
% PLOT OBJECTIVE FUNCTION EVERY MPC INSTANCE, IS IT DECREASING ? 

% To Do: (look example in Lars Grune's book)
% 1. Create a loop running dynamic optimization (use collocation and
% path-following)
% 2. Apply first-input to integrator (ode15). Set the rest control for
% initial guess in the optimizer.
% 3. Run simulator (imitating real plant)
% 4. Get the latest measurement (or initial state value)
% 5. Continue to Step 1

% number of mpc iteration
%mpciterations = 20;
%mpciterations = 5;
mpciterations = 1;
% number of prediction horizon
N             = 10;
% sampling time
T             = 1;  % [minute]
% initial controls
u0            = [2.70629; 3.20629];   % [LT; VB]
%u0            = repmat(u0,N,1);      % multiply with prediction horizon
%u0            = [2.70629; 3.20629];   % [LT VB]
% get initial measurement (states) at time T = 0.
tmeasure      = 0.0;
load cola_init.mat;
%Xinit         = 0.5*ones(82,1);
xmeasure      = Xinit;

%[t, x, u, obj] = pfNmpc(@optProblem, @system, mpciterations, N, T, tmeasure, xmeasure, u0);
%[t, x, u] = pfNmpc(@optProblem, @system, mpciterations, N, T, tmeasure, xmeasure, u0);
[t, x, u] = convNmpc(@optProblem, @system, mpciterations, N, T, tmeasure, xmeasure, u0);
%[t, x, u] = iNmpc(@optProblem, @system, mpciterations, N, T, tmeasure, xmeasure, u0);

keyboard;

% PLOT THE RESULTS !
% check if the objective function is decreasing ???
%plot(obj);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of the NMPC functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = system(t, x, u, T)
    
    % QUESTION: how long this simulation ? 
    %[~,y]=ode15s('cola_lv',[0 20000],0.5*ones(1,82)');
    global uc;
    uc = u;
    %[~,x_out] = ode15s('cola_sim',[0 1],x);    % ADD cola_sim.m to SVN!!!
    %[~,x_out] = ode15s('cola_sim_1',[0 1],x);
    [~,x_out] = ode15s('cola_sim_1',[t t+T],x);
    lengthx   = size(x_out); 
    y         = x_out(lengthx(1),:)'; 
    
end

function [J,g,w0,w,lbg,ubg,lbw,ubw] = optProblem(x, u, N)   %add prediction horizon 
    import casadi.*
    % the model
    NT = 41;
    Uf = 1.2;           % Feeding rate F
    % invoke the model
    [~,state,xdot,inputs] = DistColA(Uf);
    
    f = Function('f',{state,inputs}, {xdot});
    
    %load u_opt_ss.mat;
    
    %load cola_init.mat;
    %x_init = Xinit;
    %x_init = x;
    
    % bound constraints
    xB_max = 0.01;
    xD_min = 0.95;
    V_max  = 4.008 ;

    % State bounds and initial guess
    x_min =  zeros(82,1);
    x_max =  ones(82,1);
    x_max(1)  = xB_max;
    x_min(41) = xD_min;

    % Control bounds
    u_min = [0.1; 0.1];
    u_max = [10; V_max];
    
    % Construct objective function
    load u_opt_ss.mat; %result from FMINCON
    xdot_val_rf_ss = xf;

    % prices
    pf = 1; 
    %pV = 0.01;
    pV = 0.0001;
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
    
    % dimensions
    global nx nu nk d tf;
    nx = 82;
    nu = 2;

    nk = 50;       % control discretization
    tf = 100;     % end time
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
    beta  = 0;
    %gamma = 1e3;
    gamma = 0;

    
    %% HERE SHOULD BE LOOP N-TIMES ACCORDING TO THE NUMBER OF PREDICTION HORIZON
    
    
    % "Lift" initial conditions
    X0  = MX.sym('X0', nx);
    w   = {w{:}, X0};
    lbw = [lbw; x_min];
    ubw = [ubw; x_max];
    %w0  = [w0; xdot_val_rf_ss];   % use initial guess from steady-state optimization results. this should give zero objective function (for tracking term!)
    w0  = [w0; x];
    %g   = {g{:}, X0 - xdot_val_rf_ss};
    %g   = {g{:}, X0 - x};
    %lbg = [lbg; zeros(nx,1)];
    %ubg = [ubg; zeros(nx,1)];

    % formulate the NLP
    Xk = X0;
    for k=0:nk-1
        % New NLP variable for the control
        Uk  = MX.sym(['U_' num2str(k)], nu);
        w   = {w{:}, Uk};
        lbw = [lbw; u_min];
        ubw = [ubw; u_max];
        %w0  = [w0;  u_opt];       % optimized results from steady-state optimization
        w0  = [w0;  u];            % check u dimension !

        % State at collocation points
        Xkj = {};
        for j=1:d
            Xkj{j} = MX.sym(['X_' num2str(k) '_' num2str(j)], nx);
            w      = {w{:}, Xkj{j}};
            lbw    = [lbw; x_min];
            ubw    = [ubw; x_max];
            %w0     = [w0; xdot_val_rf_ss];
            w0     = [w0; x];
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
        %w0  = [w0; xdot_val_rf_ss];
        w0  = [w0; x];

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
        %J            = J + econ_term;
    end
end

function J = evaluateObjective(x)
    
    % Construct objective function
    load u_opt_ss.mat; %result from FMINCON
    xdot_val_rf_ss = xf;

    % prices
    pf = 1; 
    pV = 0.01; 
    pB = 1; 
    pD = 2;

    % controller gains
    KcB = 10;  
    KcD = 10;
    % Nominal holdups - these are rather small 
    MDs = 0.5; 
    MBs = 0.5;         
    % Nominal flows
    Ds  = 0.5; 
    Bs  = 0.5;
    
    % dimensions
    nx = 82;
    nu = 2;
    
    % preparing collocation matrices
    [~,C,D,d] = collocationSetup();

    nk = 10;       % control discretization
    tf = 10.0;     % end time
    h  = tf/nk;

    %delta_time = 60; % [minute] convert second to minute
    delta_time = 1;
    alpha = 1;
    beta  = 1;
    gamma = 1;

    % formulate the NLP
    %Xk = X0;
    X0      = x(1:nx,1);
    x(1:nx) = [];
    x       = reshape(x,nu + d*nx + nx, nk);
    for k=0:nk-1
        % New NLP variable for the control
        %Uk  = MX.sym(['U_' num2str(k)], nu);
        %w   = {w{:}, Uk};
        Uk = x

        % State at collocation points
        Xkj = {};
        for j=1:d
            Xkj{j} = MX.sym(['X_' num2str(k) '_' num2str(j)], nx);
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
           %[fj, qj] = easycall(f, Xkj{j},Uk);
           fj  = f(Xkj{j},Uk);
           g   = {g{:}, h*fj - xp};

           % Add contribution to the end state
           Xk_end = Xk_end + D(j+1)*Xkj{j};

        end    

        % New NLP variable for state at end of interval
        Xk  = MX.sym(['X_' num2str(k+1)], nx);
        w   = {w{:}, Xk};

        % Add equality constraint
        g   = {g{:}, Xk_end-Xk};


        % objective function
        MB = Xk(NT+1);  
        MD = Xk(2*NT);            % Actual reboiler and condenser holdup
        Do = Ds+(MD-MDs)*KcD;       % Distillate flow
        Bo = Bs+(MB-MBs)*KcB;       % Bottoms flow
        econ_term    = (pf*F + pV*Uk(2) - pB*Bo - pD*Do)*delta_time;            % economic term
        control_term = (Uk - u_opt)'*(Uk - u_opt)*delta_time;                   % tracking: control input term
        state_term   = (Xk - xdot_val_rf_ss)'*(Xk - xdot_val_rf_ss)*delta_time;  % tracking: states term
        J            = J + alpha*econ_term + beta*control_term + gamma*state_term;
    end
end
