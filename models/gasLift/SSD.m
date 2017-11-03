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


%%
[f,w_gl_0,~]    = NominalModel(par);

[F,w_gl_SP,~]   = WellSimulator(par);

%initialize simulator
xf      = dx0;
zf      = z0;
Cost    = 0;
%%
yIndex          = [1:4,7,8,23,24,25,27,29,30];

load('SRTO.mat')
h = waitbar(0,'Simulation in Progress...');

for sim_k = 1:nIter
    
    waitbar(sim_k /nIter,h,sprintf('Time: %0.0f min',sim_k*5))
    GOR_real = par.GOR + SRTO.RN(:,sim_k).*par.GOR_var;
    GOR_true(:,sim_k) = GOR_real;
    
    u_in = [w_gl_0;GOR_real];
    for EKF_k = 1:par.tf/par.tSim
        Fk = F('x0',xf,'z0',zf,'p',u_in);
        
        xf = full(Fk.xf);
        zf = full(Fk.zf);
        J_real(sim_k) = full(Fk.qf);
        
        ymeas = zf(yIndex) ;
        
        meas.p_wh(:,(sim_k-1)*300+EKF_k) = ymeas(3:4);
        meas.p_bh(:,(sim_k-1)*300+EKF_k) = ymeas(5:6);
        meas.p_rh(:,(sim_k-1)*300+EKF_k) = ymeas(9:10);
        meas.p_m(:,(sim_k-1)*300+EKF_k) = ymeas(11:12);
        meas.w_gl(:,(sim_k-1)*300+EKF_k) = u_in(1:2);
        meas.w_to((sim_k-1)*300+EKF_k) = ymeas(11);
        meas.w_tg((sim_k-1)*300+EKF_k) = ymeas(12);
        
        if sim_k>3
            pbh = meas.p_bh(1,(sim_k-1-3)*300+EKF_k : (sim_k-1)*300+EKF_k);
            pwh = meas.p_wh(1,(sim_k-1-3)*300+EKF_k : (sim_k-1)*300+EKF_k);
            wto = meas.w_to(1,(sim_k-1-3)*300+EKF_k : (sim_k-1)*300+EKF_k);
            n = length(wto);
            pbh0 = mean(pbh);
            pwh0 = mean(pwh);
            wto0 = mean(wto);
            pbhS = 0;
            pwhS = 0;
            wtoS = 0;
            pbhSd = 0;
            pwhSd = 0;
            wtoSd = 0;
            for jn = 1:n
                pbhS = pbhS + (pbh(jn) - pbh0).^2;
                pwhS = pwhS + (pwh(jn) - pwh0).^2;
                wtoS = wtoS + (wto(jn) - wto0)^2;
                if jn >1
                    pbhSd = pbhSd + (pbh(jn) - pbh(jn-1)).^2;
                    pwhSd = pwhSd + (pwh(jn) - pwh(jn-1)).^2;
                    wtoSd = wtoSd + (wto(jn) - wto(jn-1))^2;
                end
            end
            pbh_s = (1/(n-1))*pbhS;
            pwh_s = (1/(n-1))*pwhS;
            wto_s = (1/(n-1))*wtoS;
            
            pbh_sd = (1/(n-1))*pbhSd;
            pwh_sd = (1/(n-1))*pwhSd;
            wto_sd = (1/(n-1))*wtoSd;
            
            p_bhR = pbh_sd./pbh_s;
            p_whR = pwh_sd./pwh_s;
            w_toR = wto_sd./wto_s;
            
            p_bhC(:,sim_k) = pbhS;%1- 0.5.*p_bhR;
            p_whC(:,sim_k) = pwhS; %1- 0.5.*p_whR;
            w_toC(sim_k) = wtoS;%1- 0.5.*w_toR;
            
            if w_toC(sim_k) <0.001
                Var(sim_k) = 1;
            else
                Var(sim_k) = 0;
            end
            
        end
    end
end


close(h)
%%

tH = (1:SRTO.sim_k).*SRTO.par.tf./3600;
meas.tH = (1:SRTO.sim_k*SRTO.par.tf).*SRTO.par.tSim./3600;

Var(22) = 0;
Var(91) = 0;

figure(121)
clf
subplot(321)
plot(meas.tH ,meas.p_wh')
hold all
grid on
% xlim([0,SRTO.nIter])
ylabel('$p_{wh}$ [bar]','Interpreter','latex')
xlabel('time [h]','Interpreter','Latex')
leg = legend('well 1','well2');
set(leg,'Interpreter','latex','Orientation','Horizontal','Location','Best')
axs = gca;
axs.TickLabelInterpreter = 'latex';
xlim([0,12])
box on

subplot(322)
hold all
plot(meas.tH ,meas.w_tg)
plot(meas.tH,10.*ones(size(meas.w_tg)),'k:')
grid on
ylabel('Gas rate [kg/s]','Interpreter','latex')
xlabel('time [h]','Interpreter','Latex')
axs = gca;
axs.TickLabelInterpreter = 'latex';
xlim([0,12])
box on

subplot(326)
hold all
bar(tH,~Var,'k')
ylabel('SSD [0-1]','Interpreter','latex')
xlabel('time [h]','Interpreter','Latex')
axs = gca;
axs.TickLabelInterpreter = 'latex';
xlim([0,12])
box on

subplot(323)
plot(meas.tH ,meas.p_bh')
grid on
% xlim([0,SRTO.nIter])
ylabel('$p_{bh}$ [bar]','Interpreter','latex')
xlabel('time [h]','Interpreter','Latex')
leg = legend('well 1','well2');
set(leg,'Interpreter','latex','Orientation','Horizontal','Location','Best')
axs = gca;
axs.TickLabelInterpreter = 'latex';
xlim([0,12])
box on

subplot(325)
hold all
plot(tH,GOR_true')
grid on
ylabel('GOR [kg/kg]','Interpreter','latex')
xlabel('time [h]','Interpreter','Latex')
leg = legend('well 1','well2');
set(leg,'Interpreter','latex','Orientation','Horizontal','Location','Best')
axs = gca;
axs.TickLabelInterpreter = 'latex';
xlim([0,12])
box on

subplot(324)
hold all
plot(meas.tH ,meas.w_to)
grid on
ylabel('Oil rate [kg/s]','Interpreter','latex')
xlabel('time [h]','Interpreter','Latex')
axs = gca;
axs.TickLabelInterpreter = 'latex';
xlim([0,12])
box on


figure(122)
bar(tH,~Var,'k')
ylabel('SSD [0-1]','Interpreter','latex')
xlabel('time [h]','Interpreter','Latex')
axs = gca;
axs.TickLabelInterpreter = 'latex';
xlim([0,12])
box on