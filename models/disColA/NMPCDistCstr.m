function NMPCDistCstr
%NMPCDISTCSTR Summary of this function goes here
% 
% NMPC for CSTR + Distillation Column A 
%
% [OUTPUTARGS] = NMPCDISTCSTR(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2016/06/30 11:42:49 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2016

import casadi.* 
format long;

global N;
% number of mpc iteration
mpciterations = 150;
% number of prediction horizon
N             = 30;  
% sampling time
T             = 1;  % [minute]
% initial controls (different initial conditions)
%load Xinit40.mat;  
%load Xinit32.mat;
load Xinit29.mat;  
%load Xinit28.mat;
%load Xinit30.mat
%load Xinit31.mat;
%u0            = Xinit40(85:89);
u0            = Xinit29(85:89);
%u0            = Xinit32(85:89);
%u0            = Xinit28(85:89);
u0            = repmat(u0,1,N);
% get initial measurement (states) at time T = 0.
tmeasure      = 0.0;
%xmeasure      = Xinit40(1:84);
xmeasure      = Xinit29(1:84);
%xmeasure      = Xinit32(1:84);
%xmeasure      = Xinit28(1:84);

% either call iNMPC 
%[~, xmeasureAll, uAll, obj, optRes, params, runtime] = iNmpc(@optProblem, @system, mpciterations, N, T, tmeasure, xmeasure, u0);

% or pf-NMPC
[~, xmeasureAll_pf, uAll_pf, obj_pf, optRes_pf, params_pf, runtime_pf] = pfNmpc(@optProblem, @system, mpciterations, N, T, tmeasure, xmeasure, u0);
 

keyboard;

%% THE CODE BELOW IS JUST FOR PLOTTING
% load CstrDistXinit.mat;
% xf    = Xinit(1:84);
% u_opt = Xinit(85:89);
% 
% nu      = size(u0,1);
% uAll    = reshape(uAll,nu,mpciterations);
% uAll_pf = reshape(uAll_pf,nu,mpciterations);
% 
% % add initial control
% uAll    = [u0(:,1) uAll];
% uAll_pf = [u0(:,1) uAll_pf];
% 
% % add initial states
% xmeasureAll    = horzcat(xmeasure,xmeasureAll);
% xmeasureAll_pf = horzcat(xmeasure,xmeasureAll_pf);
% 
% %global nk;
% %x = linspace(1,nk*N*T,mpciterations);
% x = linspace(1,mpciterations,mpciterations/T);
% xi = [0 x];
% 
% close all; % close all figures
% 
% figure(1);
% clf;
% % figure('Units', 'pixels', ...
% %     'Position', [100 100 500 375]);
% hold on;
% hV_nlp = plot(x,obj,'LineWidth',6.0,'Color','g');
% hold on;
% hV_pf  = plot(x,obj_pf,'LineWidth',1.0,'Color','k');
% hTitle  = title ('Objective functions comparison PF - NLP');
% hXLabel = xlabel('Number of MPC iteration [-]'             );
% hYLabel = ylabel('Objective function [-]'                  );
% hLegend = legend( ...
%   [hV_nlp, hV_pf],  ...
%   'NLP: Obj. Func.' , ...
%   'PF:  Obj. Func.' , ...
%   'location', 'NorthWest' );
% 
% set( gca                       , ...
%     'FontName'   , 'Helvetica' );
% set([hTitle, hXLabel, hYLabel], ...
%     'FontName'   , 'AvantGarde');
% set([hLegend, gca]             , ...
%     'FontSize'   , 8           );
% set([hXLabel, hYLabel]  , ...
%     'FontSize'   , 10          );
% set( hTitle                    , ...
%     'FontSize'   , 12          , ...
%     'FontWeight' , 'bold'      );
% 
% set(gca, ...
%   'Box'         , 'off'     , ...
%   'TickDir'     , 'out'     , ...
%   'TickLength'  , [.02 .02] , ...
%   'XMinorTick'  , 'on'      , ...
%   'YMinorTick'  , 'on'      , ...
%   'YGrid'       , 'on'      , ...
%   'XColor'      , [.3 .3 .3], ...
%   'YColor'      , [.3 .3 .3], ...
%   'LineWidth'   , 1         );
% 
% set(gcf, 'PaperPositionMode', 'auto');
% % print -depsc2 ObjFunc.eps
% 
% 
% 
% figure(2);
% clf;
% plot(xi,xmeasureAll(1,:),'o','LineWidth',6.0,'Color','g');
% hold on;
% plot(xi,xmeasureAll_pf(1,:),'*','LineWidth',1.0,'Color','k');
% hold on;
% plot(xi, xf(1)*ones(1,mpciterations+1),'-r');
% title ('x(1) comparison PF - NLP');
% xlabel('Time [minute]'                      );
% ylabel('Concentration [-]'                  );
% 
% 
% figure(3);
% clf;
% %plot(xmeasureAll(2,:),'LineWidth',2.5);
% plot(xi,xmeasureAll(21,:),'LineWidth',6.0,'Color','g');
% hold on;
% plot(xi,xmeasureAll_pf(21,:),'LineWidth',1.0,'Color','k');
% hold on;
% plot(xi, xf(21)*ones(1,mpciterations+1),'-r');
% title ('x(21) comparison PF - NLP');
% xlabel('Time [minute]'                      );
% ylabel('Concentration [-]'                  );
% 
% 
% figure(4);
% clf;
% plot(xi,xmeasureAll(41,:),'LineWidth',6.0,'Color','g');
% hold on;
% plot(xi,xmeasureAll_pf(41,:),'LineWidth',1.0,'Color','k');
% hold on;
% plot(xi, xf(41)*ones(1,mpciterations+1),'-r');
% title ('x(41) comparison PF - NLP');
% xlabel('Time [minute]'                      );
% ylabel('Concentration [-]'                  );
% 
% figure(5);
% clf;
% %plot(xmeasureAll(3,:),'LineWidth',2.5);
% plot(xi,xmeasureAll(42,:),'LineWidth',6.0,'Color','g');
% hold on;
% plot(xi,xmeasureAll_pf(42,:),'LineWidth',1.0,'Color','k');
% hold on;
% plot(xi, xf(42)*ones(1,mpciterations+1),'-r');
% title ('x(42) comparison PF - NLP');
% xlabel('Time [minute]'                      );
% ylabel('Concentration [-]'                  );
% 
% figure(6);
% clf;
% %plot(xmeasureAll(3,:),'LineWidth',2.5);
% plot(xi,xmeasureAll(82,:),'LineWidth',6.0,'Color','g');
% hold on;
% plot(xi,xmeasureAll_pf(82,:),'LineWidth',1.0,'Color','k');
% hold on;
% plot(xi, xf(82)*ones(1,mpciterations+1),'-r');
% title ('x(82) comparison PF - NLP');
% xlabel('Time [minute]'                      );
% ylabel('Hold-up [-]'                        );
% 
% figure(7);
% clf;
% plot(xi,uAll(1,:),'LineWidth',6.0,'Color','g');
% hold on;
% plot(xi,uAll_pf(1,:),'LineWidth',1.0,'Color','k');
% hold on;
% plot(xi, u_opt(1,:)*ones(1,mpciterations+1),'-r');
% title ('u(1)-LT: control input comparison PF - NLP');
% xlabel('Time [minute]'                          );
% ylabel('LT [m^3/minute]'                        );
% 
% figure(8);
% clf;
% plot(xi,uAll(2,:),'LineWidth',6.0,'Color','g');
% hold on;
% plot(xi,uAll_pf(2,:),'LineWidth',1.0,'Color','k');
% hold on;
% plot(xi, u_opt(2,:)*ones(1,mpciterations+1),'-r');
% title ('u(2)-VB: control input comparison PF - NLP');
% xlabel('Time [minute]'                          );
% ylabel('VB [kmol/minute]'                        );
% 
% figure(9);
% clf;
% plot(xi,uAll(3,:),'LineWidth',6.0,'Color','g');
% hold on;
% plot(xi,uAll_pf(3,:),'LineWidth',1.0,'Color','k');
% hold on;
% plot(xi, u_opt(3,:)*ones(1,mpciterations+1),'-r');
% title ('u(3)-F: control input comparison PF - NLP');
% xlabel('Time [minute]'                          );
% ylabel('F [kmol/minute]'                        );
% 
% figure(10);
% clf;
% plot(xi,uAll(4,:),'LineWidth',6.0,'Color','g');
% hold on;
% plot(xi,uAll_pf(4,:),'LineWidth',1.0,'Color','k');
% hold on;
% plot(xi, u_opt(4,:)*ones(1,mpciterations+1),'-r');
% title ('u(4)-D: control input comparison PF - NLP');
% xlabel('Time [minute]'                          );
% ylabel('D [kmol/minute]'                        );
% 
% figure(11);
% clf;
% plot(xi,uAll(5,:),'LineWidth',6.0,'Color','g');
% hold on;
% plot(xi,uAll_pf(5,:),'LineWidth',1.0,'Color','k');
% hold on;
% plot(xi, u_opt(5,:)*ones(1,mpciterations+1),'-r');
% title ('u(5)-B: control input comparison PF - NLP');
% xlabel('Time [minute]'                          );
% ylabel('B [kmol/minute]'                        );
% 
% keyboard;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of the NMPC functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = system(t, x, u, T)
    
    global uc;
    uc = u;
    [~,x_out] = ode15s('cola_lv_cstr',[t t+T], x);
    lengthx   = size(x_out); 
    y         = x_out(lengthx(1),:)'; 
    
end

function [J,g,w0,w,lbg,ubg,lbw,ubw,params] = optProblem(x, u, N, x0_measure)   %add prediction horizon 
    import casadi.*
    
    % the model
    NT = 41;
    Uf = 0.3;           % Feeding rate F_0
    
    % invoke the model
    [~,state,xdot,inputs] = DistColACstr(Uf);
    f = Function('f',{state,inputs}, {xdot});
    
    % bound constraints
    VB_max = 4.008;
    
    % State bounds and initial guess
    x_min =  zeros(84,1);  % try without epsilon here, later put epsilon
    x_max =  ones(84,1);

    % Control bounds
    u_min = [0.1; 0.1; 0.1; 0.1; 0.1];
    u_max = [10; VB_max; 10; 1.0; 1.0];
    
    % compact bound variable for ease of function invocation 
    params.bound.x_min = x_min;
    params.bound.x_max = x_max;
    params.bound.u_min = u_min;
    params.bound.u_max = u_max;
    
    % Construct objective function
    load CstrDistXinit.mat;
    xf    = Xinit(1:84);
    u_opt = Xinit(85:89);
    
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
    
    % dimensions
    global nx nu nk d tf ns;
    nx = 84;   % CSTR + Distillation Column A
    nu = 5;    % LT, VB, F, D, B
    nk = 1;
    tf = 1;      % in [minutes]
    h  = tf/nk;
    ns = 0;
    
    % compact model variable
    params.model.NT = NT;
    params.model.f  = f;
    params.model.xdot_val_rf_ss = xf;
    params.model.x  = x;
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
    w0  = [w0; x(1,1:nx)'];
    g   = {g{:}, X0 - x0_measure};  % USE MEASUREMENT HERE !
    lbg = [lbg; zeros(nx,1)];
    ubg = [ubg; zeros(nx,1)];

    % formulate the NLP
    Xk = X0;

    
    load Qmax.mat;
    params.Qmax = Qmax;
    
    % HERE SHOULD BE LOOP N-TIMES ACCORDING TO THE NUMBER OF PREDICTION HORIZON
    count  = 2; % counter for state variable as initial guess
    ssoftc = 0;
    for i=1:N
        [J,g,w0,w,lbg,ubg,lbw,ubw,Xk,params,count,ssoftc] = iterateOnPredictionHorizon(Xk, w, w0, lbw, ubw, lbg, ubg, g, J, params, i, count,ssoftc);
    end

    
end


function [J,g,w0,w,lbg,ubg,lbw,ubw,Xk,params,count,ssoftc] = iterateOnPredictionHorizon(Xk, w, w0, lbw, ubw, lbg, ubg, g, J, params, iter, count, ssoftc)

   import casadi.*
   global N;
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
   
   global nx nu nk d ns;

   for k=0:nk-1
        % New NLP variable for the control
        %Uk  = MX.sym(['U_' num2str(k)], nu);
        Uk     = MX.sym(['U_' num2str((iter-1)*nk+k)], nu);
        w      = {w{:}, Uk};
        lbw    = [lbw; u_min];
        ubw    = [ubw; u_max];
        indexU = (iter-1)*nk + (k+1);
        w0     = [w0;  u(:,indexU)];
        
        Jcontrol   = (Qmax(nx+1:nx+nu,1).*(Uk - u_opt))' * (Uk - u_opt);
      
        % State at collocation points
        Xkj   = {};
        SumX1 = 0;
        for j=1:d
            Xkj{j} = MX.sym(['X_' num2str((iter-1)*nk+k) '_' num2str(j)], nx);
            w      = {w{:}, Xkj{j}};
            lbw    = [lbw; x_min];
            ubw    = [ubw; x_max];   
            w0     = [w0; x(iter+1,:)'];  
            count  = count + 1;
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
           
        end    

        % New NLP variable for state at end of interval
        Xk  = MX.sym(['X_' num2str((iter-1)*nk+k)], nx);
        w   = {w{:}, Xk};
        lbw = [lbw; x_min];
        x_maxEnd       =  ones(84,1);
        x_maxEnd(1,1)  = 0.1;
        x_maxEnd(84,1) = 0.7;
        ubw = [ubw; x_maxEnd];
        w0  = [w0; x(iter+1,:)'];
        count  = count + 1;

        % Add equality constraint
        g   = {g{:}, (Xk_end-Xk)};
        lbg = [lbg; zeros(nx,1)];
        ubg = [ubg; zeros(nx,1)];
               
        Jecon  = (pf*F_0 + pV*Uk(2) - pB*Uk(5) - pD*Uk(4)) * delta_time;
        Jstate =(Qmax(1:nx,1).*(Xk - xdot_val_rf_ss))' * (Xk - xdot_val_rf_ss) * delta_time;
        
        alpha  = 1;
        beta   = 1;
        gamma  = 1;

        
        % compute rotated cost function
        fm  = f(Xk,Uk);
        % load Lagrange multipliers from steady-state optimization
        load LamdaCstrDist.mat; % lamda
        Jmodel = lamda'*fm;
        
        J = J + alpha*Jcontrol + gamma*Jstate + beta*Jecon;
    end
end


