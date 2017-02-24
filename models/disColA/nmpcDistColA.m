function nmpcDistColA
%NMPCDISTCOLA Summary of this function goes here
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

% number of mpc iteration
mpciterations = 40;
%mpciterations = 20;
%mpciterations = 10;
%mpciterations = 5;  % WORK !
%mpciterations = 3;
%mpciterations = 1;
% number of prediction horizon
%N             = 8;  % WORK !
%N             = 10;
N             = 5;
%N             = 1;
% sampling time
T             = 1;  % [minute]
% initial controls
%u0            = [2.70629 3.20629];   % [LT VB]
%u0            = repmat(u0,N,1);      % multiply with prediction horizon
u0            = [2.70629; 3.20629];   % [LT VB]
%u0            = [3.3779; 4.0080];    % steady-state result for Greshgorin
u0            = repmat(u0,1,N);
% get initial measurement (states) at time T = 0.
tmeasure      = 0.0;
load xinit.mat;
xmeasure      = Xinit;

% call NMPC 
[~, xmeasureAll, uAll, obj, optRes, params] = iNmpc(@optProblem, @system, mpciterations, N, T, tmeasure, xmeasure, u0);
%[~, xmeasureAll_pf, uAll_pf, obj_pf, optRes, params] = pfNmpc(@optProblem, @system, mpciterations, N, T, tmeasure, xmeasure, u0);

% save iNMPC.mat xmeasureAll uAll obj;
% 
% keyboard;

load iNMPC.mat;     % load iNMPC results
load u_opt_ss.mat;       % load steady-state results

nu      = size(u0,1);
uAll    = reshape(uAll,nu,mpciterations);
uAll_pf = reshape(uAll_pf,nu,mpciterations);

% add initial control
uAll    = [u0(:,1) uAll];
uAll_pf = [u0(:,1) uAll_pf];

% add initial states
xmeasureAll    = horzcat(xmeasure,xmeasureAll);
xmeasureAll_pf = horzcat(xmeasure,xmeasureAll_pf);

global nk;
x = linspace(1,nk*N*T,mpciterations);
xi = [0 x];

figure(5);
clf;
figure('Units', 'pixels', ...
    'Position', [100 100 500 375]);
hold on;
hV_nlp = plot(x,obj,'LineWidth',6.0,'Color','g');
hold on;
hV_pf  = plot(x,obj_pf,'LineWidth',1.0,'Color','k');
hTitle  = title ('Objective functions comparison PF - NLP');
hXLabel = xlabel('Time [second]'                      );
hYLabel = ylabel('Objective function [-]'                  );
hLegend = legend( ...
  [hV_nlp, hV_pf],  ...
  'NLP: Obj. Func.' , ...
  'PF:  Obj. Func.' , ...
  'location', 'NorthWest' );

set( gca                       , ...
    'FontName'   , 'Helvetica' );
set([hTitle, hXLabel, hYLabel], ...
    'FontName'   , 'AvantGarde');
set([hLegend, gca]             , ...
    'FontSize'   , 8           );
set([hXLabel, hYLabel]  , ...
    'FontSize'   , 10          );
set( hTitle                    , ...
    'FontSize'   , 12          , ...
    'FontWeight' , 'bold'      );

set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1         );

set(gcf, 'PaperPositionMode', 'auto');
print -depsc2 ObjFunc.eps



figure(6);
clf;
plot(xi,xmeasureAll(1,:),'LineWidth',6.0,'Color','g');
hold on;
plot(xi,xmeasureAll_pf(1,:),'LineWidth',1.0,'Color','k');
hold on;
plot(xi, xf(1)*ones(1,mpciterations+1),'-r');


figure(7);
clf;
%plot(xmeasureAll(2,:),'LineWidth',2.5);
plot(xi,xmeasureAll(21,:),'LineWidth',6.0,'Color','g');
hold on;
plot(xi,xmeasureAll_pf(21,:),'LineWidth',1.0,'Color','k');
hold on;
plot(xi, xf(21)*ones(1,mpciterations+1),'-r');


figure(8);
clf;
%plot(xmeasureAll(3,:),'LineWidth',2.5);
plot(xi,xmeasureAll(42,:),'LineWidth',6.0,'Color','g');
hold on;
plot(xi,xmeasureAll_pf(42,:),'LineWidth',1.0,'Color','k');
hold on;
plot(xi, xf(42)*ones(1,mpciterations+1),'-r');

figure(9);
clf;
%plot(xmeasureAll(3,:),'LineWidth',2.5);
plot(xi,xmeasureAll(82,:),'LineWidth',6.0,'Color','g');
hold on;
plot(xi,xmeasureAll_pf(82,:),'LineWidth',1.0,'Color','k');
hold on;
plot(xi, xf(82)*ones(1,mpciterations+1),'-r');

figure(10);
clf;
plot(xi,uAll(1,:),'LineWidth',6.0,'Color','g');
hold on;
plot(xi,uAll_pf(1,:),'LineWidth',1.0,'Color','k');
hold on;
plot(xi, u_opt(1,:)*ones(1,mpciterations+1),'-r');

figure(11);
clf;
plot(xi,uAll(2,:),'LineWidth',6.0,'Color','g');
hold on;
plot(xi,uAll_pf(2,:),'LineWidth',1.0,'Color','k');
hold on;
plot(xi, u_opt(2,:)*ones(1,mpciterations+1),'-r');


% % EVALUATE ECONOMIC OBJECTIVE FUNCTION FROM OPTIMIZATION
% [economic_obj, track_u, track_x] = evaluateEconObj(optRes, params, N);
% fprintf('Economic objective value: %f\n', economic_obj);
% fprintf('Input tracking value: %f\n', sum(track_u));
% fprintf('State tracking value: %f\n', sum(track_x));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of the NMPC functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = system(t, x, u, T)
    
    global uc;
    uc = u;
    [~,x_out] = ode15s('cola_sim_1',[t t+T],x);
    lengthx   = size(x_out); 
    y         = x_out(lengthx(1),:)'; 
    
end

function [J,g,w0,w,lbg,ubg,lbw,ubw,params] = optProblem(x, u, N)   %add prediction horizon 
    import casadi.*
    
    % the model
    NT = 41;
    Uf = 1.2;           % Feeding rate F
    
    % invoke the model
    [~,state,xdot,inputs] = DistColA(Uf);
    f = Function('f',{state,inputs}, {xdot});
    
    % bound constraints
    xB_max = 0.01;
    xD_min = 0.95;
    V_max  = 4.008 ;

    % State bounds and initial guess
    x_min =  zeros(82,1);
    x_max =  ones(82,1);
    %x_min =  1e-5*ones(82,1);
    %x_max =  0.9999*ones(82,1);
    x_max(1)  = xB_max;
    x_min(41) = xD_min;
    
    % soft-constraint for state variables
    epsilon = 0.1;  
%     x_min   = x_min - epsilon*ones(82,1);
%     x_max   = x_max + epsilon*ones(82,1);

    % Control bounds
    u_min = [0.1; 0.1];
    u_max = [10; V_max];
    
    % compact bound variable for ease of function invocation 
    params.bound.x_min = x_min;
    params.bound.x_max = x_max;
    params.bound.u_min = u_min;
    params.bound.u_max = u_max;
    
    % Construct objective function
    load u_opt_ss.mat; %result from FMINCON
    %xdot_val_rf_ss = xf;
    
    % prices
    pf = 1; 
    %pV = 0.01;
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
    
    % dimensions
    global nx nu nk d tf;
    nx = 82;
    nu = 2;

    %nk = 50;      % control discretization
    %tf = 100;     % end time
    %nk = 1;
    nk = 2;
    tf = 5;        % in [minutes]
    h  = tf/nk;
    
    % compact model variable
    params.model.NT = NT;
    params.model.f  = f;
    params.model.xdot_val_rf_ss = xf;
    params.model.x  = x;
    %model.u  = u;
    params.model.u_opt = u_opt;
    % duplicate u0
    params.model.u  = repmat(u,1,nk); 


    % preparing collocation matrices
    [~,C,D,d] = collocationSetup();
    
    % compact collocation variable
    params.colloc.C = C;
    params.colloc.D = D;
    params.colloc.h = h;

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
    %gamma = 0;
    
    % compact weight variable
    params.weight.delta_time = delta_time;
    params.weight.alpha      = alpha;
    params.weight.beta       = beta;
    params.weight.gamma      = gamma;
    
    % "Lift" initial conditions
    X0  = MX.sym('X0', nx);
    w   = {w{:}, X0};
    lbw = [lbw; x_min];
    ubw = [ubw; x_max];
    w0  = [w0; x(1,:)'];     
    g   = {g{:}, X0 - x(1,:)'};
    lbg = [lbg; zeros(nx,1)];
    ubg = [ubg; zeros(nx,1)];

    % formulate the NLP
    Xk = X0;
    %load weights2.mat;
    %params.weights = weights;
    load Qmax.mat;
    params.Qmax = Qmax;
    
    % objective function at initial state value
    J   = J +( Qmax(1:nx,1).*(X0 - xf))' * (X0 - xf);
    
    % HERE SHOULD BE LOOP N-TIMES ACCORDING TO THE NUMBER OF PREDICTION HORIZON
    for i=1:N
        [J,g,w0,w,lbg,ubg,lbw,ubw,Xk,params] = iterateOnPredictionHorizon(Xk, w, w0, lbw, ubw, lbg, ubg, g, J, params, i, epsilon);
    end

    
end


function [J,g,w0,w,lbg,ubg,lbw,ubw,Xk,params] = iterateOnPredictionHorizon(Xk, w, w0, lbw, ubw, lbg, ubg, g, J, params, iter, epsilon)

   import casadi.*
   
   % extract compact variables
   x_min = params.bound.x_min;
   x_max = params.bound.x_max;
   u_min = params.bound.u_min;
   u_max = params.bound.u_max;
   
   NT = params.model.NT;
   f  = params.model.f;
   xdot_val_rf_ss = params.model.xdot_val_rf_ss;
   x = params.model.x;
   u = params.model.u;
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
   %load weights.mat;
   %load weights1.mat;
   
   %count = 1;
   for k=0:nk-1
        % New NLP variable for the control
        %Uk  = MX.sym(['U_' num2str(k)], nu);
        Uk     = MX.sym(['U_' num2str((iter-1)*nk+k)], nu);
        w      = {w{:}, Uk};
        lbw    = [lbw; u_min];
        ubw    = [ubw; u_max];
        indexU = (iter-1)*nk + (k+1);
        w0     = [w0;  u(:,indexU)];
        
        J   = J + (Qmax(nx+nu,1).*(Uk - u_opt))' * (Uk - u_opt);
        
%         % terminal set constraint
%         if k == (nk-1)
%             g   = {g{:}, u_opt-Uk};
%             lbg = [lbg; zeros(nu,1)];
%             ubg = [ubg; zeros(nu,1)];
%         end

        % State at collocation points
        Xkj = {};
        for j=1:d
            %Xkj{j} = MX.sym(['X_' num2str(k) '_' num2str(j)], nx);
            Xkj{j} = MX.sym(['X_' num2str((iter-1)*nk+k) '_' num2str(j)], nx);
            w      = {w{:}, Xkj{j}};
            lbw    = [lbw; x_min];
            ubw    = [ubw; x_max];
            w0     = [w0; x(iter+1,:)'];   
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
           fj  = f(Xkj{j},Uk);
           g   = {g{:}, h*fj - xp};
           lbg = [lbg; zeros(nx,1)];
           ubg = [ubg; zeros(nx,1)];

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
%            J = J + (Uk - u_opt)' * (Uk - u_opt) * delta_time;
%            J = J + (Xkj{j} - xdot_val_rf_ss)' * (Xkj{j} - xdot_val_rf_ss) * delta_time;
           
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
        lbw = [lbw; x_min];
        ubw = [ubw; x_max];
        %w0  = [w0; xdot_val_rf_ss];
        w0  = [w0; x(iter+1,:)'];

        % Add equality constraint
        g   = {g{:}, Xk_end-Xk};
        %g   = {g{:}, xdot_val_rf_ss-Xk};
        lbg = [lbg; zeros(nx,1)];
        ubg = [ubg; zeros(nx,1)];
        
        
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
%         rho = 1e3;
%         J = J + rho*epsilon;

    end
end

function [J, tu, tx] = evaluateEconObj(optRes, params, N)   %add prediction horizon 

    % Construct objective function
    load u_opt_ss.mat; %result from FMINCON
    
    % prices
    pf = params.price.pf; 
    pV = params.price.pV;
    pB = params.price.pB; 
    pD = params.price.pD;
    F  = params.price.F;

    % controller gains
    KcB = params.gain.KcB;  
    KcD = params.gain.KcD;
    MDs = params.gain.MDs; 
    MBs = params.gain.MBs;         
    Ds  = params.gain.Ds; 
    Bs  = params.gain.Bs;
    
    % dimensions
    global nx nu nk d;
    
    % compact model variable
    NT    = params.model.NT;
    
    % compact weight variable
    delta_time = params.weight.delta_time;

    % formulate the NLP
    X0           = optRes(1:nx);
    optRes(1:nx) = [];
    optRes       = reshape(optRes,(nu + nx*d + nx), N*nk);
    
    J  = 0;
    tu = zeros(nu,1);
    tx = zeros(nx,1);
    xdot_val_rf_ss = params.model.xdot_val_rf_ss;
    for i=1:N
        Xki = optRes(:,i);
        Uk  = optRes(1:nu);
        Xk  = Xki(end-nx+1:end);
        % objective function
        MB = Xk(NT+1);
        MD = Xk(2*NT);
        Do = Ds+(MD-MDs)*KcD;
        Bo = Bs+(MB-MBs)*KcB;
        J  = J + (pf*F + pV*Uk(2) - pB*Bo - pD*Do)*delta_time;
        
        for k=1:(nu+nx)
            if k<=nu
                tu(k) = params.weights(k) * ( Uk(k) - u_opt(k) )^2;
            else
                tx(k-nu) = params.weights(k) * ( Xk(k-nu) - xdot_val_rf_ss(k-nu) )^2;
            end
        end
    end
    
end
