function NMPCGasLift
%NMPCGASLIFT Summary of this function goes here
% 
% [OUTPUTARGS] = NMPCGASLIFT(INPUTARGS) Explain usage here
%
% NMPC for Gas Lift
%
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2017/10/24 21:18:07 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2017

import casadi.* 
format long;

global N;
mpciterations = 60;
% number of prediction horizon
N             = 60;  
% sampling time
T             = 1;  % [minute]

%% First initial states and controls !
% All the parameter values are defined inside this function
param       = GasLiftRiser_Param;
param.sv    = 0;        % use slack variable? (1=Yes; 0=No)
param.N     = 24;           % no. of samples in prediction horizon
param.T     = 2*3600;       % 2 hrs predition horizon
param.tf    = param.T/param.N;  % each sample is 5 min
param.tSim  = 1;                % Shorter simulation time for EKF
param.nIter = 60;               % number of closed loop iterations (60-->5h of CL simulation)


PlotResult  = 1;        % plot prediction horizon at each sampling interval
EKF         = 1;        % use EKF? (1=Yes; 0=No)

param.sv          = 0;        % use slack variable? (1=Yes; 0=No)
param.slip_real   = 1.0;      
param.GOR_real    = param.GOR;

% import initial conditions
%[dx0,z0,u0,lbx,lbz,lbu,ubx,ubz,ubu] = GasLiftRiser_Initialization_bounds(param);
[dx0,z0,u0,param.lbx,param.lbz,param.lbu,param.ubx,param.ubz,param.ubu] = GasLiftRiser_Initialization_bounds(param);

% initial states and controls
load xuInit.mat;
tmeasure = 0.0;
xmeasure = [xInit(1:8);xInit(9:end)];

% either call iNMPC 
%[~, xmeasureAll, uAll, obj, optRes, params, runtime] = iNmpc(@optProblem, @system, mpciterations, N, T, tmeasure, xmeasure, u0);
%[~, xmeasureAll, uAll, obj, optRes, params, runtime] = OutputINmpc(@optProblem, @system, mpciterations, N, T, tmeasure, xmeasure, u0, param);

% or pf-NMPC
%[~, xmeasureAll_pf, uAll_pf, obj_pf, optRes_pf, params_pf, runtime_pf] = pfNmpc(@optProblem, @system, mpciterations, N, T, tmeasure, xmeasure, u0);
[~, xmeasureAll_pf, uAll_pf, obj_pf, optRes_pf, params_pf, runtime_pf] = OutputPfNmpc(@optProblem, @system, mpciterations, N, T, tmeasure, xmeasure, u0, param);
 

keyboard;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of the NMPC functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = system(t, x, u, T, par)
    % real plant simulation using IDAS integrator
    F  = wellModel(par);
    nu = par.n_w;
    %nz = 30;
    nx = 8;
    
    % initialization Kalman Filter 
    EKF  = 1;
    dx0  = x(1:nx);
    u0   = u(1:nu);
    z0   = x(nx+1:end);
    %u_in = u;
    u_in    = [u;par.GOR];
    
    [f_EKF,JacFx,h_EKF,JacHx,z_EKF,yIndex,nxEKF,~,~,~,~,~] = EKF_integrating(par);
    nyEKF   = length(yIndex);
    xk_hat  = [dx0];
    uEKF    = u0;
    
    Pk = 1e3.*eye(nxEKF);
    Qk = 1e3.*eye(nxEKF);
    Rk = 1e0.*eye(nyEKF);
    P = Pk;
    
    xf      = dx0;
    zf      = z0;
    Cost    = 0;
    
    for EKF_k = 1:par.tf/par.tSim

        Fk = F('x0',xf,'z0',zf,'p',u_in);
        xf = full(Fk.xf);
        zf = full(Fk.zf);
        %J_real(sim_k) = full(Fk.qf);

        ymeas = zf(yIndex) + (randn(nyEKF,1).*0.0);

%         meas.p_wh(:,(sim_k-1)*300+EKF_k)    = ymeas(3:4);
%         meas.p_bh(:,(sim_k-1)*300+EKF_k)    = ymeas(5:6);
%         meas.p_rh(:,(sim_k-1)*300+EKF_k)    = ymeas(7);
%         meas.p_m(:,(sim_k-1)*300+EKF_k)     = ymeas(8);
%         meas.w_gl(:,(sim_k-1)*300+EKF_k)    = u_in(1:2);
%         meas.w_to((sim_k-1)*300+EKF_k)      = ymeas(9);
%         meas.w_tg((sim_k-1)*300+EKF_k)      = ymeas(10);

        Cost         = Cost +  ymeas(9);
        %iCost(sim_k) = Cost;

        if EKF
            % Extended Kalman filter for state estimation

            Fj = full(JacFx(xk_hat,uEKF));

            if max(max(isnan(Fj)))
                disp('NaN in x EKF')
            end

            xk_hat_1 =  full(f_EKF(xk_hat,uEKF));
            Pk_1 = Fj*Pk*Fj' + Qk;

            %uEKF    = [u_in_1;u_in_2];
            uEKF    = u(1:nu);

            Hj      = full(JacHx(xk_hat_1,uEKF));
            ek      = full(ymeas - h_EKF(xk_hat_1,uEKF));
            Sk      = Hj*Pk_1*Hj' + Rk;
            Kk      = (Pk_1*Hj')/(Sk);
            xk_hat  = xk_hat_1 + Kk*ek;
            Pk      = (eye(nxEKF) - Kk*Hj)*Pk_1;
            zk_hat  = full(z_EKF(xk_hat,uEKF));

            x_hat   = xk_hat(1:8);
            %             GOR_hat = xk_hat(9:10);
            %             GOR_est(:,sim_k) = GOR_hat;

            if isnan(xk_hat)
                disp('NaN in y EKF')
            end

            % Estimation Error
            xEstErr = abs(full(x_hat) - xf);
            zEstErr = abs(full(zk_hat) - zf);

            % set new initial values for the next iteration
            dx0     =  full(x_hat) + 0.0*randn(1);
            z0      = full(zk_hat) + 0.0*randn(1);
            %u0      = [u_in_1;u_in_2];

        else
            % set new initial values for the next iteration
            dx0     =  xf;
            z0      = zf;
            %u0      = [u_in_1;u_in_2];
            %J_real(sim_k) = full(Fk.qf);
        end
    end
    y = [dx0;z0];
end

function [J,g,w0,w,p,lbg,ubg,lbw,ubw,par] = optProblem(x, u, N, x0_measure, par)   %add prediction horizon
% TO DO: - add equality constraint for initial state variable values

import casadi.*
global nu nk tf ns;
nk = par.N;
tf = 1;      % in [minutes]
ns = 0;
nu = 2;
nz = 30;
nx = 8;
% extract u0, z0, and dx0
u0  = u(1:nu,1);
dx0 = x(1:nx);
z0  = x(nx+1:end);
% bounds
lbx = par.lbx;
lbz = par.lbz;
lbu = par.lbu;
ubx = par.ubx;
ubz = par.ubz;
ubu = par.ubu;

% nomimal model
f = nomModel(par);

%% Direct Collocation
[B,C,D,d] = collocationSetup();

%% Build NLP solver
% empty nlp
w   = {};
w0  = [];
lbw = [];
ubw = [];
J   = 0;

g   = {};
lbg = [];
ubg = [];

% initial conditions for each scenario
X0  = MX.sym('X0',nx);
Z0  = MX.sym('Z0',nz);
w   = {w{:}, X0,Z0};
lbw = [lbw;lbx;lbz];
ubw = [ubw;ubx;ubz];
w0  = [w0; dx0;z0];

% Formulate NLP
Xk  = X0;
Xkj = {};
Zkj = {};
% Uk_prev = U0;

p   = MX.sym('p',nu);
g   = {g{:},[X0;Z0]-x0_measure};
lbg = [lbg;zeros(nx+nz,1)];
ubg = [ubg;zeros(nx+nz,1)];
load QmaxGL.mat;
load xuSS.mat;

for k = 0:par.N-1
    
    Uk  = MX.sym(['U_' num2str(k)],nu);
    w   = {w{:},Uk};
    lbw = [lbw;lbu];
    ubw = [ubw;ubu];
    w0  = [w0;u0];
    
    % regularization term for control input
    Jcontrol   = (QmaxGL(nx+nz+1:nx+nz+nu,1).*(Uk - uSS))' * (Uk - uSS);
    
    Xkj = {};
    Zkj = {};
    
    for j = 1:d
        Xkj{j} = MX.sym(['X_' num2str(k) '_' num2str(j)],nx);
        Zkj{j} = MX.sym(['Z_' num2str(k) '_' num2str(j)],nz);
        w   = {w{:},Xkj{j},Zkj{j}};
        lbw = [lbw;lbx;lbz];
        ubw = [ubw;ubx;ubz];
        w0  = [w0; dx0;z0 ];
    end
    
    % Loop over collocation points
    Xk_end  = D(1)*Xk;
    
    for j = 1:d
        % Expression for the state derivative at the collocation point
        xp  = C(1,j+1)*Xk;  % helper state
        for r = 1:d
            xp = xp + C(r+1,j+1)*Xkj{r};
        end
        [fj,zj,qj] =  f(Xkj{j},Zkj{j},vertcat(Uk,p));
        
        g   = {g{:},par.tf*fj-xp,zj};  % dynamics and algebraic constraints
        lbg = [lbg;zeros(nx,1);zeros(nz,1)];
        ubg = [ubg;zeros(nx,1);zeros(nz,1)];
        
        % Add contribution to the end states
        Xk_end  = Xk_end + D(j+1)*Xkj{j};
        
    end
    
    
    % New NLP variable for state at end of interval
    Xk      = MX.sym(['X_' num2str(k+1) ], nx);
    w       = {w{:},Xk};
    lbw     = [lbw;lbx];
    ubw     = [ubw;ubx];
    w0      = [w0; dx0];
 
    % Shooting Gap constraint
    g   = {g{:},Xk_end-Xk};
    lbg = [lbg;zeros(nx,1)];
    ubg = [ubg;zeros(nx,1)];
        
%     % Gas capacity constraints
%     g   = {g{:},sum(Zkj{j}(17:18))};
%     lbg = [lbg;0];
%     ubg = [ubg;par.QgMax];
    
%     % Gas Lift constraint
%     g   = {g{:},sum(Uk)};
%     lbg = [lbg;0];
%     ubg = [ubg;par.qGLMax];
    
    % regularization term for state variable
    Jstate =(QmaxGL(1:nx+nz,1).*([Xk;Zkj{j}] - xSS))' * ([Xk;Zkj{j}] - xSS);
    % economic objective state cost
    Jecon  = -Zkj{j}(29) + sum(Uk);
    % stage cost    
    J = J + Jcontrol + Jstate + Jecon;
    
end

end