clear 
clc

% Import CasADi
addpath ('C:\Users\dineshk\CasADi\casadi-matlabR2014b-v3.1.0-rc1')
import casadi.*
addpath('C:\matlab')

%% Set parameters, initial values, upper and lower bounds

% All the parameter values are defined inside this function
par = GasLiftRiser_Param;
n_w = par.n_w; % no. of wells;

[dx0,z0,u0,lbx,lbz,lbu,ubx,ubz,ubu] = GasLiftRiser_Initialization_bounds(par);
GOR_val = par.GOR;
PI_val = par.PI;

%% Simulation Setup
par.N       = 24;%60;
par.T       = 2*3600;
par.tf      = par.T/par.N;
par.tSim    = 1;
par.nIter   = 144;

EKF             = 1;
par.CostMeas    = 0;

par.sv          = 1;
par.QgMax       = 10;
par.qGLMax      = 50;

par.PI_real     = par.PI;
par.RN          = 0*[0.75;0.5];
par.slip_real   = 1.0;
par.GOR_real    = par.GOR + par.RN.*par.GOR_var;

sv = par.sv;
nIter = par.nIter;
%% Nominal Model

R       = par.R;
Mw      = par.Mw;
L_w     = par.L_w;
H_w     = par.H_w;
D_w     = par.D_w;

L_bh    = par.L_bh;
H_bh    = par.H_bh;
D_bh    = par.D_bh;

L_a     = par.L_a;
H_a     = par.H_a;
D_a     = par.D_a;

L_r     = par.L_r;
H_r     = par.H_r;
D_r     = par.D_r;

rho_o   = par.rho_o;
C_iv    = par.C_iv;
C_pc    = par.C_pc;
C_pr    = par.C_pr;
rho_ro  = sum(rho_o)/2;

mu_oil  = 1*0.001; % 1cP oil viscosity

A_w     = pi.*(D_w/2).^2;
A_bh    = pi.*(D_bh/2).^2;
V_a     = L_a.*(pi.*(D_a/2).^2 - pi.*(D_w/2).^2);
A_r     = pi.*(D_r/2).^2;

% differential states
m_ga    = MX.sym('m_ga',n_w); % 1-2
m_gt    = MX.sym('m_gt',n_w); % 3-4
m_ot    = MX.sym('m_ot',n_w); % 5-6
m_gr    = MX.sym('m_gr',1);   % 7
m_or    = MX.sym('m_or',1);   % 8

% control input
w_gl    = MX.sym('w_gl',n_w);

% parameters
p_res   = MX.sym('p_res',n_w);
PI      = MX.sym('PI',n_w);
GOR     = MX.sym('GOR',n_w);
T_a     = MX.sym('T_a',n_w);
T_w     = MX.sym('T_w',n_w);
T_r     = MX.sym('T_r',1);
p_s     = MX.sym('p_s',1);

slip    = 1.0;

xGwH    = (m_gt.*1e3./max(1e-3,(m_gt.*1e3+m_ot.*1e3)));
xOwH    = (m_ot.*1e3./max(1e-3,(m_gt.*1e3+m_ot.*1e3)));
xGrH    = (m_gr.*1e3./(m_gr.*1e3+m_or.*1e3));
xOrH    = (m_or.*1e3./(m_gr.*1e3+m_or.*1e3));

xGw     = slip.*xGwH./(1 + (slip-1).*xGwH);
xOw     = 1 - xGw;
xGr     = slip.*xGrH./(1 + (slip-1).*xGrH);
xOr     = 1 - xGr;

d1      = MX.sym('d1',n_w); % 1-2
d2      = MX.sym('d2',n_w); % 3-4
d3      = MX.sym('d3',n_w); % 5-6
d4      = MX.sym('d4',1);   % 7
d5      = MX.sym('d5',1);   % 8

% algebraic equations used for substitution in the ODE model
p_rh    = 1e-5.*(((R.*par.T_r./Mw).*(m_gr.*1e3./(L_r.*A_r))) - ((m_gr.*1e3+m_or.*1e3 )./(L_r.*A_r)).*9.81.*H_r/2); %27
rho_r   = 1e-2.*((m_gr.*1e3 + m_or.*1e3).*p_rh.*1e5.*Mw.*rho_ro)./(m_or.*1e3.*p_rh.*1e5.*Mw + rho_ro.*R.*par.T_r.*m_gr.*1e3); % 28
w_pr    = C_pr.*sqrt(rho_r.*1e2.*(p_rh.*1e5 - par.p_s.*1e5)); %30
p_m     = 1e-5.*(p_rh.*1e5 + 9.81./(A_r.*L_r).*(m_or.*1e3+m_gr.*1e3).*H_r + 128.*mu_oil.*L_r.*w_pr./(3.14.*D_r.^4.*((m_gr.*1e3 + m_or.*1e3).*p_rh.*1e5.*Mw.*rho_ro)./(m_or.*1e3.*p_rh.*1e5.*Mw + rho_ro.*R.*par.T_r.*m_gr.*1e3))) ;%29
p_ai    = 1e-5.*(((R.*par.T_a./(V_a.*Mw) + 9.81.*H_a./V_a).*m_ga.*1e3) + (Mw./(R.*par.T_a).*((R.*par.T_a./(V_a.*Mw) + 9.81.*H_a./V_a).*m_ga.*1e3)).*9.81.*H_a);
p_wh    = 1e-5.*(((R.*par.T_w./Mw).*(m_gt.*1e3./(L_w.*A_w + L_bh.*A_bh - m_ot.*1e3./rho_o))) - ((m_gt.*1e3+m_ot.*1e3 )./(L_w.*A_w)).*9.81.*H_w/2);
rho_ai  = 1e-2.*(Mw./(R.*par.T_a).*p_ai.*1e5);
rho_m   = 1e-2.*(((m_gt.*1e3 + m_ot.*1e3).*p_wh.*1e5.*Mw.*rho_o)./(m_ot.*1e3.*p_wh.*1e5.*Mw + rho_o.*R.*par.T_w.*m_gt.*1e3));
w_pc    = C_pc.*sqrt(rho_m.*1e2.*(p_wh.*1e5 - p_m.*1e5));
w_pg    = xGw.*w_pc;
w_po    = xOr.*w_pc;
p_wi    = 1e-5.*((p_wh.*1e5 + 9.81./(A_w.*L_w).*max(0,(m_ot.*1e3+m_gt.*1e3-rho_o.*L_bh.*A_bh)).*H_w + 128.*mu_oil.*L_w.*w_pc./(3.14.*D_w.^4.*((m_gt.*1e3 + m_ot.*1e3).*p_wh.*1e5.*Mw.*rho_o)./(m_ot.*1e3.*p_wh.*1e5.*Mw + rho_o.*R.*par.T_w.*m_gt.*1e3))));
p_bh    = 1e-5.*(p_wi.*1e5 + rho_o.*9.81.*H_bh + 128.*mu_oil.*L_bh.*w_po./(3.14.*D_bh.^4.*rho_o));
w_iv    = C_iv.*sqrt(rho_ai.*1e2.*(p_ai.*1e5 - p_wi.*1e5));
w_ro    = (PI).*1e-6.*(par.p_res.*1e5 - p_bh.*1e5);% w_ro = (-IPR.a + sqrt(IPR.a.^2+4*IPR.b.*(p_res-p_bh).*1e5))./(2.*IPR.b);
w_rg    = 1e1.*GOR.*w_ro;
w_to    = xOr.*w_pr; %29
w_tg    = xGr.*w_pr; %30

% differential equations
df1     = (w_gl - w_iv).*1e-3 + d1;
df2     = (w_iv + w_rg.*1e-1 - w_pg).*1e-3 + d2;
df3     = (w_ro - w_po).*1e-3 + d3;
df4     = (sum(w_pg) - w_tg).*1e-3 + d4 ;
df5     = (sum(w_po) - w_to).*1e-3 + d5 ;

% Form the DAE system
diff    = vertcat(df1,df2,df3,df4,df5);

% give parameter values
diff    = substitute(diff,p_res,par.p_res);
diff    = substitute(diff,p_s,par.p_s);
diff    = substitute(diff,T_a,par.T_a);
diff    = substitute(diff,T_w,par.T_w);
diff    = substitute(diff,T_r,par.T_r);

% concatenate the differential and algebraic states
x_var   = vertcat(m_ga,m_gt,m_ot,m_gr,m_or);
p_var   = vertcat(w_gl,GOR,PI,d1,d2,d3,d4,d5);

% Objective function
L       = -(w_to.^2) + 0.5.*sum((w_gl).^2);

z_vec   = vertcat(p_ai,p_wh,p_wi,p_bh,rho_ai,rho_m,w_iv,w_pc,w_pg,w_po,w_ro,w_rg,p_rh,rho_r,p_m,w_pr,w_to,w_tg);
zF      = Function('zF',{x_var,p_var},{z_vec});

w_tg_eval = Function('w_tg_eval',{x_var,p_var},{w_tg});
%% Simulator model

[F,w_gl_SP,~] 	= WellSimulator(par);

%initialize simulator
xf      = dx0;
zf      = z0;
Cost    = 0;

disp([' True Steady state optimum'])
disp(' ')
disp(['Well 1: ' num2str(w_gl_SP(1)) ' [kg/s]'])
disp(['Well 2: ' num2str(w_gl_SP(2)) ' [kg/s]'])
disp(['Total : ' num2str(sum(w_gl_SP)) ' [kg/s]'])
disp(' ------------------------------------------- ')
disp(' ------------------------------------------- ')

%% For extended Kalman Filter
if EKF
    
    [f_EKF,JacFx,JacFw,h_EKF,JacHx,z_EKF,yIndex,nxEKF,JacAx,JacBu,JacJx,JacJu,J_EKF] = EKF_integrating(par);
    
    nyEKF   = length(yIndex);
    d0      = zeros(8,1);
    xk_hat  = [dx0;par.GOR;par.PI;d0];
    uEKF    = u0;
    
    PIestimate = 0;
    IntDist    = 1;
    
    Px  = 1e3.*eye(12);
    Pd  = IntDist.*1e3.*eye(8);
    Qx  = 1e3.*eye(12);
    Qd  = IntDist.*1e1.*eye(8);
    if ~PIestimate
        Qx(11,11) = 0;
        Qx(12,12) = 0;
        Px(11,11) = 0;
        Px(12,12) = 0;
    end
    Pk  = [Px,zeros(12,8);zeros(8,12),Pd];
    Qk  = [Qx,zeros(12,8);zeros(8,12),Qd];
    Rk  = 1e3.*eye(nyEKF);
    P   = Pk;
    K   = 0;
    
    GOR_hat = par.GOR;
    PI_hat  = par.PI;
    d_hat   = d0;
end
%% Setpoint tracking NMPC

MPCinit.dx0      = dx0;
MPCinit.z0       = z0;
MPCinit.u0       = u0;
MPCinit.u_in     = u0;
[solverNMPC,MPC] = SetpointNMPC(par,MPCinit);

%%
nu = par.n_w;
nz = 30;
nx = 8;
nd = nx;
np = nu;

%% Build Static NLP

% decision variables
w   = {};
w0  = [];
lbw = [];
ubw = [];

% constraints
g   = {};
lbg = [];
ubg = [];

w   = {w{:},x_var,p_var};
lbw = [lbw;lbx;lbu;par.GOR;par.PI;d0];
ubw = [ubw;ubx;ubu;par.GOR;par.PI;d0];
w0  = [w0;dx0;u0;par.GOR;par.PI;d0];

J   = L; % Economic objective

%Add the system model as constraints
g   = {g{:},vertcat(diff)};
lbg = [lbg;zeros(8,1)];
ubg = [ubg;zeros(8,1)];

% Add gas lift and gas capacity constraints
g   = {g{:},sum(p_var(1:2)),w_tg_eval(x_var,p_var)};
lbg = [lbg;0;0];
ubg = [ubg;par.qGLMax;par.QgMax];

nlp = struct('x',vertcat(w{:}),'f',J,'g',vertcat(g{:}));

solver = nlpsol('solver','ipopt',nlp);


%% Simulation

h = waitbar(0,'Simulation in Progress...');
for sim_k = 1:nIter
    waitbar(sim_k /nIter,h,sprintf('Time: %0.0f min',sim_k*5))
    
    if rem(sim_k,17) == 0
        par.RN(1) = rand - 0.5;
    end
    if rem(sim_k,23) ==0 || sim_k == 7
        par.RN(2) = rand - 0.5;
    end
    
    RN(:,sim_k)         = par.RN;    
    GOR_real            = par.GOR + par.RN.*par.GOR_var;
    GOR_true(:,sim_k)   = GOR_real;
    
    tic;
    sol         = solver('x0',w0,'lbx',lbw,'ubx',ubw,'lbg',lbg,'ubg',ubg);
    t(sim_k)    = toc;
    NMPC.t      = t;
    
    w_opt       = full(sol.x);
    J_MPC(sim_k)= full(sol.f);
    
    %%  extract solution
    
    u_opt1 = w_opt(nx+1);
    u_opt2 = w_opt(nx+2);
    
    u_in_RTO1(sim_k)  = u_opt1;
    u_in_RTO2(sim_k)  = u_opt2;
    
    NMPC.u_in_RTO = [u_in_RTO1;u_in_RTO2];
    
    %% SetpointTracking NMPC
    
    u_SP = [u_opt1;u_opt2];
    
    MPC.w0 = [];
    MPC.lbw = [];
    MPC.ubw = [];
    
    MPC.w0  = [MPC.w0;dx0;z0];
    MPC.lbw = [MPC.lbw;dx0;z0];
    MPC.ubw = [MPC.ubw;dx0;z0];
    MPC.w0  = [MPC.w0;u0];
    MPC.lbw = [MPC.lbw;u0];
    MPC.ubw = [MPC.ubw;u0];
    for i = 1:par.N
        MPC.w0  = [MPC.w0;u0;u_SP;GOR_hat;PI_hat;d_hat];
        MPC.lbw = [MPC.lbw;lbu;u_SP;GOR_hat;PI_hat;d_hat];
        MPC.ubw = [MPC.ubw;ubu;u_SP;GOR_hat;PI_hat;d_hat];
        for d = 1:3
            MPC.w0  = [MPC.w0;dx0;z0];
            MPC.lbw = [MPC.lbw;lbx;lbz];
            MPC.ubw = [MPC.ubw;ubx;ubz];
        end
        for d = 1:3
            if sv
                % slack  variable
                MPC.w0  = [MPC.w0;0];
                MPC.lbw = [MPC.lbw;0];
                MPC.ubw = [MPC.ubw;2];
            end
        end
        MPC.w0  = [MPC.w0;dx0];
        MPC.lbw = [MPC.lbw;lbx];
        MPC.ubw = [MPC.ubw;ubx];
    end
    
    solNMPC   = solverNMPC('x0',MPC.w0,'lbx',MPC.lbw,'ubx',MPC.ubw,'lbg',MPC.lbg,'ubg',MPC.ubg);
    w_optMPC  = full(solNMPC.x);
    
    n_w_i     = nx+nz+nu+ (nu+nu+nu+np+nd + (nx+nz+sv)*d+nx)*(par.N) ;
    
    u_opt1MPC = [w_optMPC((nx+nz+nu+1):d*(nx+nz+sv)+nx+2*nu+nu+np+nd:n_w_i);NaN];
    u_opt2MPC = [w_optMPC((nx+nz+nu+2):d*(nx+nz+sv)+nx+2*nu+nu+np+nd:n_w_i);NaN];
    
    u_in      = [u_opt1MPC(1,1); u_opt2MPC(1,1);GOR_real];
    
    u_in_MPC1(sim_k)    = u_opt1MPC(1,1);
    u_in_MPC2(sim_k)    = u_opt2MPC(1,1);
    NMPC.u_in_MPC       = [u_in_MPC1;u_in_MPC2];
    %% Simulator using IDAS integrator
    
    for EKF_k = 1:par.tf/par.tSim
        
        Fk = F('x0',xf,'z0',zf,'p',u_in);
        xf = full(Fk.xf);
        zf = full(Fk.zf);
        J_real(sim_k) = full(Fk.qf);
        
        ymeas = zf(yIndex) + (randn(nyEKF,1).*0.0);
        
        meas.p_wh(:,(sim_k-1)*300+EKF_k)    = ymeas(3:4);
        meas.p_bh(:,(sim_k-1)*300+EKF_k)    = ymeas(5:6);
        meas.p_rh(:,(sim_k-1)*300+EKF_k)    = ymeas(9:10);
        meas.p_m(:,(sim_k-1)*300+EKF_k)     = ymeas(11:12);
        meas.w_gl(:,(sim_k-1)*300+EKF_k)    = u_in(1:2);
        meas.w_to((sim_k-1)*300+EKF_k)      = ymeas(11);
        meas.w_tg((sim_k-1)*300+EKF_k)      = ymeas(12);
        
        Cost         = Cost +  ymeas(11);
        iCost(sim_k) = Cost;        
        
        if EKF
            % Extended Kalman filter for state estimation
            wp = 0.00000001.*randn(4+nd,1);
            
            Fj = full(JacFx(xk_hat,uEKF,wp));
            
            if max(max(isnan(Fj)))
                disp('NaN in x EKF')
            end
            
            xk_hat_1 =  full(f_EKF(xk_hat,uEKF,wp));
            Pk_1 = Fj*Pk*Fj' + Qk;
            
            Hj      = full(JacHx(xk_hat_1,uEKF));
            ek      = full(ymeas - h_EKF(xk_hat_1,uEKF));
            Sk      = Hj*Pk_1*Hj' + Rk;
            Kk      = (Pk_1*Hj')/(Sk);
            xk_hat  = xk_hat_1 + Kk*ek;            
            Pk      = (eye(nxEKF) - Kk*Hj)*Pk_1;
            zk_hat  = full(z_EKF(xk_hat,uEKF));
            
            x_hat   = xk_hat(1:8);
            GOR_hat = xk_hat(9:10);
            PI_hat  = xk_hat(11:12);
            d_hat   = xk_hat(13:20);
            
            GOR_est(:,sim_k) = GOR_hat;
            if isnan(xk_hat)
                disp('NaN in y EKF')
            end
            
            % Estimation Error
            xEstErr = abs(full(x_hat) - xf);
            zEstErr = abs(full(zk_hat) - zf);
            
            % set new initial values for the next iteration
            dx0     =  full(x_hat);
            u0      = u_in(1:2);
            z0      = full(zk_hat);
            uEKF    = u_in(1:2);
            
            A       = full(JacAx(xk_hat,uEKF));
            B       = full(JacBu(xk_hat,uEKF));
            C       = full(JacJx(xk_hat,uEKF));
            D       = full(JacJu(xk_hat,uEKF));
            Ju_hat  = -C*(A\B) + D;
            Ju(:,sim_k) = Ju_hat;
            
        else
            % set new initial values for the next iteration
            dx0     =  xf;
            u0      = u_in(1:2);
            z0      = zf;
            J_real(sim_k) = full(Fk.qf);
            
            GOR_hat = par.GOR;
            PI_hat  = par.PI;
            d_hat   = d0;
            
        end
    end
    
    %% Set new initial values for next iteration
    
    w0  = [];
    lbw = [];
    ubw = [];
    
    lbw = [lbw;lbx;lbu;GOR_hat;PI_hat;d_hat];
    ubw = [ubw;ubx;ubu;GOR_hat;PI_hat;d_hat];
    w0  = [w0;dx0;u0;GOR_hat;PI_hat;d_hat];
    
end
close(h)

h1 = waitbar(0,'Simulating Ideal Case...');
for sim_k = 1:nIter
    waitbar(sim_k /nIter,h1)
    par.slip_real   = 1.0;
    par.GOR_real    = par.GOR + RN(:,sim_k).*par.GOR_var;
    
    [F,w_gl_SP,~]   = WellSimulator(par);
    u_opt_ideal(:,sim_k) = w_gl_SP;
end
close(h1)

%%

SRTO.GOR_true   = GOR_true;
SRTO.GOR_est    = GOR_est;
SRTO.par        = par;
SRTO.meas       = meas;
SRTO.NMPC       = NMPC;
SRTO.sim_k      = sim_k;
SRTO.nIter      = nIter;
SRTO.iCost      = iCost;
SRTO.RN         = RN;
SRTO.u_opt_ideal= u_opt_ideal;
save('SRTO','SRTO')

GOR_RN = SRTO.RN;
save('GOR_RN','GOR_RN')

%%

load('SRTO.mat')

tH = (1:SRTO.sim_k).*SRTO.par.tf./3600; 
meas.tH = (1:SRTO.sim_k*SRTO.par.tf).*SRTO.par.tSim./3600; 

figure(124)
clf
subplot(221)
stairs(meas.tH,SRTO.meas.p_wh')
hold all
grid on
ylabel('wellhead pressure[bar]','Interpreter','latex')


subplot(222)
hold all
plot(meas.tH,SRTO.meas.w_tg)
grid on
ylabel('Oil rate [kg/s]','Interpreter','latex')

subplot(224)
hold all
plot(tH,SRTO.GOR_true')
plot(tH,SRTO.GOR_est','--')
grid on
ylabel('GOR [kg/kg]','Interpreter','latex')

subplot(223)
plot(meas.tH,SRTO.meas.p_bh')
grid on
ylabel('bottom hole pressure [bar]','Interpreter','latex')


figure(125)
clf
subplot(221)
stairs(tH,SRTO.NMPC.u_in_MPC(1,:))
hold all
stairs(tH,SRTO.NMPC.u_in_MPC(2,:))
stairs(tH,SRTO.NMPC.u_in_RTO(1,:))
stairs(tH,SRTO.NMPC.u_in_RTO(2,:))
stairs(tH,SRTO.u_opt_ideal','--')
grid on
ylabel('GasLift rate rate [kg/s]','Interpreter','latex')

subplot(222)
hold all
plot(meas.tH,SRTO.meas.w_tg)
grid on
ylabel('Oil rate [kg/s]','Interpreter','latex')

subplot(223)
hold all
plot(tH,SRTO.GOR_true')
plot(tH,SRTO.GOR_est','--')
grid on
ylabel('GOR [kg/kg]','Interpreter','latex')

subplot(224)
hold all
plot(meas.tH,SRTO.meas.w_to)
grid on
ylabel('Oil rate [kg/s]','Interpreter','latex')

%%



figure(221)
clf
subplot(321)
hold all
% stairs(tH,SRTO.NMPC.u_in_RTO(1,:))
% stairs(tH,SRTO.NMPC.u_in_RTO(2,:))
stairs(tH,SRTO.NMPC.u_in_MPC(1,:))
stairs(tH,SRTO.NMPC.u_in_MPC(2,:))
% stairs(tH,SRTO.u_opt_ideal','--')
grid on
ylabel('GasLift rate rate [kg/s]','Interpreter','latex')
xlabel('time [h]','Interpreter','Latex')
leg = legend('well 1 ','well 2','well 1 - MPC','well 2 - MPC');
set(leg,'Interpreter','latex','Orientation','Horizontal','Location','Best')
axs = gca;
axs.TickLabelInterpreter = 'latex';
xlim([0,12])
ylim([0,3.2])
box on

subplot(322)
hold all
plot(meas.tH,10.*ones(size(SRTO.meas.w_tg)),'k:')
plot(meas.tH,SRTO.meas.w_tg)
grid on
ylabel('Gas rate [kg/s]','Interpreter','latex')
xlabel('time [h]','Interpreter','Latex')
axs = gca;
axs.TickLabelInterpreter = 'latex';
xlim([0,12])
box on

subplot(323)
hold all
plot(tH,SRTO.GOR_true')
plot(tH,SRTO.GOR_est','--')
grid on
ylabel('GOR [kg/kg]','Interpreter','latex')
xlabel('time [h]','Interpreter','Latex')
leg = legend('well 1','well 2 ');
set(leg,'Interpreter','latex','Orientation','Horizontal','Location','Best')
axs = gca;
axs.TickLabelInterpreter = 'latex';
xlim([0,12])
ylim([0.07,0.13])
box on

subplot(324)
hold all
plot(meas.tH,SRTO.meas.w_to)
grid on
ylabel('Oil rate [kg/s]','Interpreter','latex')
xlabel('time [h]','Interpreter','Latex')
axs = gca;
axs.TickLabelInterpreter = 'latex';
xlim([0,12])
ylim([60,75])
box on

subplot(325)
stairs(meas.tH,SRTO.meas.p_wh')
hold all
grid on
ylabel('$p_{wh}$ [bar]','Interpreter','latex')
xlabel('time [h]','Interpreter','Latex')
leg = legend('well 1','well2');
set(leg,'Interpreter','latex','Orientation','Horizontal','Location','Best')
axs = gca;
axs.TickLabelInterpreter = 'latex';
xlim([0,12])
ylim([46,54])
box on


subplot(326)
plot(meas.tH,SRTO.meas.p_bh')
grid on
ylabel('$p_{bh}$ [bar]','Interpreter','latex')
xlabel('time [h]','Interpreter','Latex')
leg = legend('well 1','well2');
set(leg,'Interpreter','latex','Orientation','Horizontal','Location','Best')
axs = gca;
axs.TickLabelInterpreter = 'latex';
xlim([0,12])
box on

%%
figure(222)
clf
hold all
stairs(tH,SRTO.u_opt_ideal','--')
stairs(tH,SRTO.NMPC.u_in_RTO(1,:))
stairs(tH,SRTO.NMPC.u_in_RTO(2,:))
stairs(tH,SRTO.NMPC.u_in_MPC(1,:))
stairs(tH,SRTO.NMPC.u_in_MPC(2,:))
grid on
ylabel('GasLift rate rate [kg/s]','Interpreter','latex')
xlabel('time [h]','Interpreter','Latex')
leg = legend('well 1 - ideal','well 2 - ideal','well 1 - SRTO','well 2 - SRTO','well 1 - MPC','well 2 - MPC');
set(leg,'Interpreter','latex','Orientation','Vertical','Location','Best')
axs = gca;
axs.TickLabelInterpreter = 'latex';
xlim([0,12])
% ylim([0,3.0])
box on

%%
figure(95)
plot(tH,abs(SRTO.GOR_true'- SRTO.GOR_est'))
grid on
ylabel('GOR estimation Error','Interpreter','latex')
xlabel('time [h]','Interpreter','Latex')
leg = legend('well 1 ','well 2 ','well 1 - SRTO','well 2 - SRTO');
set(leg,'Interpreter','latex','Orientation','Vertical','Location','Best')
axs = gca;
axs.TickLabelInterpreter = 'latex';


%%




figure(85)
clf
subplot(312)
plot(meas.tH,SRTO.meas.w_to,'k')
grid on
ylabel('Oil rate [kg/s]','Interpreter','latex')
xlabel('time [h]','Interpreter','Latex')
axs = gca;
axs.TickLabelInterpreter = 'latex';
xlim([0,12])
ylim([60,75])
box on

subplot(311)
hold all
plot(meas.tH,10.*ones(size(SRTO.meas.w_tg)),'k:')
plot(meas.tH,SRTO.meas.w_tg,'k')
grid on
ylabel('Gas rate [kg/s]','Interpreter','latex')
xlabel('time [h]','Interpreter','Latex')
axs = gca;
axs.TickLabelInterpreter = 'latex';
xlim([0,12])
ylim([7,12])
box on

subplot(313)
plot(meas.tH,SRTO.meas.p_bh')
grid on
ylabel(' bottom hole pressure [bar]','Interpreter','latex')
xlabel('time [h]','Interpreter','Latex')
leg = legend('well 1','well 2');
set(leg,'Interpreter','latex')
axs = gca;
axs.TickLabelInterpreter = 'latex';
xlim([0,12])
box on

%%
figure(23)
clf
subplot(211)
stairs(meas.tH,SRTO.meas.p_rh(1,:))
grid on
ylabel('well head pressure [bar]','Interpreter','latex')
xlabel('time [h]','Interpreter','Latex')
leg = legend('well 1','well 2');
set(leg,'Interpreter','latex')
axs = gca;
axs.TickLabelInterpreter = 'latex';
xlim([0,12])
% ylim([46,54])
box on

% 


%%

figure()
hold all
plot(tH,SRTO.GOR_true')
plot(tH,SRTO.GOR_est','--')
grid on
ylabel('GOR [kg/kg]','Interpreter','latex')
xlabel('time [h]','Interpreter','Latex')
leg = legend('well 1 true','well 2 true', 'well 1 Est', 'well 2 Est');
set(leg,'Interpreter','latex','Orientation','Horizontal','Location','Best')
axs = gca;
axs.TickLabelInterpreter = 'latex';
xlim([0,12])
ylim([0.07,0.13])
box on