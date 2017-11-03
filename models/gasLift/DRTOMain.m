clear 
clc

% Import CasADi
%addpath ('C:\Users\dineshk\CasADi\casadi-matlabR2014b-v3.1.0-rc1')
import casadi.*


%% Set parameters, initial values, upper and lower bounds

% All the parameter values are defined inside this function
par = GasLiftRiser_Param;
n_w = par.n_w; % no. of wells;

% import initial conditions
[dx0,z0,u0,lbx,lbz,lbu,ubx,ubz,ubu] = GasLiftRiser_Initialization_bounds(par);


%% Simulation Setup
par.N       = 24;           % no. of samples in prediction horizon
par.T       = 2*3600;       % 2 hrs predition horizon
par.tf      = par.T/par.N;  % each sample is 5 min
par.tSim    = 1;            % Shorter simulation time for EKF
par.nIter   = 60;           % number of closed loop iterations (60-->5h of CL simulation)

PlotResult      = 0;        % plot prediction horizon at each sampling interval
EKF             = 1;        % use EKF? (1=Yes; 0=No)

par.sv          = 0;        % use slack variable? (1=Yes; 0=No)

par.slip_real   = 1.0;      
par.GOR_real    = par.GOR;

sv = par.sv;
nIter = par.nIter;

%%
[f,w_gl_0,~] = NominalModel(par);

[F,w_gl_SP,~]  = WellSimulator(par);

nu = par.n_w;
nz = 30;
nx = 8;

%initialize simulator
xf      = dx0;
zf      = z0;
Cost    = 0;

disp(' ------------------------------------------- ')
disp(' ------------------------------------------- ')
disp([' Nominal Steady state optimum'])
disp(' ')
disp(['Well 1: ' num2str(w_gl_0(1)) ' [kg/s]'])
disp(['Well 2: ' num2str(w_gl_0(2)) ' [kg/s]'])
disp(['Total : ' num2str(sum(w_gl_0)) ' [kg/s]'])
disp(' ------------------------------------------- ')
disp([' True Steady state optimum'])
disp(' ')
disp(['Well 1: ' num2str(w_gl_SP(1)) ' [kg/s]'])
disp(['Well 2: ' num2str(w_gl_SP(2)) ' [kg/s]'])
disp(['Total : ' num2str(sum(w_gl_SP)) ' [kg/s]'])
disp(' ------------------------------------------- ')
disp(' ------------------------------------------- ')

%% For extended Kalman Filter
if EKF
    
    [f_EKF,JacFx,h_EKF,JacHx,z_EKF,yIndex,nxEKF,~,~,~,~,~] = EKF_integrating(par);
    
    nyEKF   = length(yIndex);
    xk_hat  = [dx0];
    uEKF    = u0;
    
    Pk = 1e3.*eye(nxEKF);
    Qk = 1e3.*eye(nxEKF);
    Rk = 1e0.*eye(nyEKF);
    P = Pk;
    K = 0;
    
end
%% Direct Collocation

% Degree of interpolating polynomial
d = 3;

% Get collocation points
tau_root = [0, collocation_points(d, 'radau')];

% Coefficients of the collocation equation
C = zeros(d+1,d+1);

% Coefficients of the continuity equation
D = zeros(d+1, 1);

% Coefficients of the quadrature function
B = zeros(d+1, 1);

% Construct polynomial basis
for j=1:d+1
    % Construct Lagrange polynomials to get the polynomial basis at the collocation point
    coeff = 1;
    for r=1:d+1
        if r ~= j
            coeff = conv(coeff, [1, -tau_root(r)]);
            coeff = coeff / (tau_root(j)-tau_root(r));
        end
    end
    % Evaluate the polynomial at the final time to get the coefficients of the continuity equation
    D(j) = polyval(coeff, 1.0);
    
    % Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
    pder = polyder(coeff);
    for r=1:d+1
        C(j,r) = polyval(pder, tau_root(r));
    end
    
    % Evaluate the integral of the polynomial to get the coefficients of the quadrature function
    pint = polyint(coeff);
    B(j) = polyval(pint, 1.0);
end


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
lbw = [lbw; dx0;z0];
ubw = [ubw; dx0;z0];
w0  = [w0; dx0;z0];

U0  = MX.sym('U0',2);
w   = {w{:}, U0};
lbw = [lbw; u0];
ubw = [ubw; u0];
w0  = [w0; u0];

% Formulate NLP
Xk  = X0;
Xkj = {};
Zkj = {};
Uk_prev = U0;

for k = 0:par.N-1
    
    Uk      = MX.sym(['U_' num2str(k)],nu);
    GOR_k   = MX.sym(['GOR_' num2str(k)],nu);   
    Upar    = vertcat(Uk,GOR_k);
    
    w   = {w{:},Upar};
    lbw = [lbw;lbu;par.GOR];
    ubw = [ubw;ubu;par.GOR];
    w0  = [w0;u0;par.GOR];
    
    Xkj = {};
    Zkj = {};
    
    for j = 1:d
        Xkj{j} = MX.sym(['X_' num2str(k) '_' num2str(j)],nx);
        Zkj{j} = MX.sym(['Z_' num2str(k) '_' num2str(j)],nz);
        if par.sv
            s{j} = MX.sym(['s_' num2str(k) '_' num2str(j)],1);
        end
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
        [fj,zj,qj] =  f(Xkj{j},Zkj{j},Upar);
        
        g   = {g{:},par.tf*fj-xp,zj};  % dynamics and algebraic constraints
        lbg = [lbg;zeros(nx,1);zeros(nz,1)];
        ubg = [ubg;zeros(nx,1);zeros(nz,1)];
        
        % Gas capacity constraints on all the collocation points
        if par.sv
            g   = {g{:},Zkj{j}(30)-s{j}};  %
        else
            g   = {g{:},Zkj{j}(30)};
        end
        lbg     = [lbg;0];
        ubg     = [ubg;par.QgMax];
        
        % Slack variables for gas capacity constraints
        if par.sv
            w   = {w{:},s{j}};
            lbw = [lbw;0];
            ubw = [ubw;1];
            w0  = [w0;0];
        end
        
        % Add contribution to the end states
        Xk_end  = Xk_end + D(j+1)*Xkj{j};
        
        if par.sv
            J   = J + (B(j+1)*qj*par.tf ) + 100000*s{j} + 10.*sum((Uk_prev - Uk).^2) ;
        else
            J   = J + (B(j+1)*qj*par.tf ) + 10.*sum((Uk_prev - Uk).^2) ; %
        end
    end
    
    Uk_prev = MX.sym(['Uprev_' num2str(k+1)],nu);
    Uk_prev = Uk;
    
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
    
    % Gas Lift constraint
    g   = {g{:},sum(Uk)};
    lbg = [lbg;0];
    ubg = [ubg;par.qGLMax];
    
    % MV Rate Of Change Constraints
    g = {g{:},(Uk_prev - Uk)};
    lbg = [lbg;-par.qGLROC];
    ubg = [ubg;par.qGLROC];
    
end

% create and solve NLP solver
nlp     = struct('x',vertcat(w{:}),'f',J,'g',vertcat(g{:}));
solver  = nlpsol('solver','ipopt',nlp);


%% Simulation
h = waitbar(0,'Simulation in Progress...');

for sim_k = 1:nIter
    waitbar(sim_k /nIter,h,sprintf('Time: %0.0f min',sim_k*5))
    
    GOR_true(:,sim_k)   = par.GOR;
    
    tic;
    sol         = solver('x0',w0,'lbx',lbw,'ubx',ubw,'lbg',lbg,'ubg',ubg);
    t(sim_k)    = toc;
    NMPC.t      = t; % save computation time for each iteration
    
    w_opt           = full(sol.x);
    J_MPC(sim_k)    = full(sol.f);
    
    %%  extract solution
    
    n_w_i       = nx+nz+nu+ (nu+nu + (nx+nz+sv)*d+nx)*(par.N) ;    
    u_opt1      = [w_opt((nx+nz+nu+1):d*(nx+nz+sv)+nx+2*nu:n_w_i);NaN];
    u_opt2      = [w_opt((nx+nz+nu+2):d*(nx+nz+sv)+nx+2*nu:n_w_i);NaN];    
    w_po1       = w_opt([(nx+nz-11),nu+(nx+nz-11)+(nx+nz+nu+nu):d*(nx+nz+sv)+nx+nu+nu:n_w_i]);
    w_po2       = w_opt([(nx+nz-10),nu+(nx+nz-10)+(nx+nz+nu+nu):d*(nx+nz+sv)+nx+nu+nu:n_w_i]);
    w_oEx       = w_opt([(nx+nz-1),nu+(nx+nz-1)+(nx+nz+nu+nu):d*(nx+nz+sv)+nx+nu+nu:n_w_i]);
    s_opt       = w_opt([nx+nz+nu+(nu+nu)+(nx+nz)*d+1:d*(nx+nz+sv)+nx+nu+nu:n_w_i]);    
    w_pg1       = w_opt([(nx+nz-13),nu+(nx+nz-13)+(nx+nz+nu+nu):d*(nx+nz+sv)+nx+nu+nu:n_w_i]);
    w_pg2       = w_opt([(nx+nz-12),nu+(nx+nz-12)+(nx+nz+nu+nu):d*(nx+nz+sv)+nx+nu+nu:n_w_i]);
    w_gEx       = w_opt([(nx+nz-0),nu+(nx+nz-0)+(nx+nz+nu+nu):d*(nx+nz+sv)+nx+nu+nu:n_w_i]);    
    p_bh1       = w_opt([(nx+nz-23),nu+(nx+nz-23)+(nx+nz+nu+nu):d*(nx+nz+sv)+nx+nu+nu:n_w_i]);
    p_bh2       = w_opt([(nx+nz-22),nu+(nx+nz-22)+(nx+nz+nu+nu):d*(nx+nz+sv)+nx+nu+nu:n_w_i]);
    
    % implement the first sample on the simulator
 
    u_in_1  = u_opt1(1) ;
    u_in_2  = u_opt2(1) ;
    u_in    = [u_in_1;u_in_2;par.GOR];
    
    u_in_MPC1(sim_k)  = u_in_1;
    u_in_MPC2(sim_k)  = u_in_2;
    NMPC.u_in_MPC     = [u_in_MPC1;u_in_MPC2];
    
    
    %% Simulator using IDAS integrator
    
    for EKF_k = 1:par.tf/par.tSim
        
        Fk = F('x0',xf,'z0',zf,'p',u_in);
        xf = full(Fk.xf);
        zf = full(Fk.zf);
        J_real(sim_k) = full(Fk.qf);
        
        ymeas = zf(yIndex) + (randn(nyEKF,1).*0.0);
        
        meas.p_wh(:,(sim_k-1)*300+EKF_k)    = ymeas(3:4);
        meas.p_bh(:,(sim_k-1)*300+EKF_k)    = ymeas(5:6);
        meas.p_rh(:,(sim_k-1)*300+EKF_k)    = ymeas(7);
        meas.p_m(:,(sim_k-1)*300+EKF_k)     = ymeas(8);
        meas.w_gl(:,(sim_k-1)*300+EKF_k)    = u_in(1:2);
        meas.w_to((sim_k-1)*300+EKF_k)      = ymeas(9);
        meas.w_tg((sim_k-1)*300+EKF_k)      = ymeas(10);
        
        Cost         = Cost +  ymeas(9);
        iCost(sim_k) = Cost;
        
        if EKF          
            % Extended Kalman filter for state estimation
            
            Fj = full(JacFx(xk_hat,uEKF));
            
            if max(max(isnan(Fj)))
                disp('NaN in x EKF')
            end
            
            xk_hat_1 =  full(f_EKF(xk_hat,uEKF));
            Pk_1 = Fj*Pk*Fj' + Qk;
            
            uEKF    = [u_in_1;u_in_2];
            
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
            u0      = [u_in_1;u_in_2];
            
        else
            % set new initial values for the next iteration
            dx0     =  xf;
            z0      = zf;
            u0      = [u_in_1;u_in_2];
            J_real(sim_k) = full(Fk.qf);                   
        end       
    end
    
    %% Set new initial values for next iteration
     
    w0 = [];
    lbw = [];
    ubw = [];
        w0 = [w0;dx0;z0];
        lbw = [lbw;dx0;z0];
        ubw = [ubw;dx0;z0];
        w0 = [w0;u0];
        lbw = [lbw;u0];
        ubw = [ubw;u0];
        for i = 1:par.N
            w0 = [w0;u0;par.GOR];
            lbw = [lbw;lbu;par.GOR];
            ubw = [ubw;ubu;par.GOR];
            for d = 1:3
                w0 = [w0;dx0;z0];
                lbw = [lbw;lbx;lbz];
                ubw = [ubw;ubx;ubz];
            end
            for d = 1:3
                if sv
                    w0 = [w0;0];
                    lbw = [lbw;0];
                    ubw = [ubw;1];
                end
            end
            w0 = [w0;dx0];
            lbw = [lbw;lbx];
            ubw = [ubw;ubx];
        end
        
    %% Plot results
    
    w_po_MPC1(sim_k)    = full(Fk.zf(19));
    w_po_MPC2(sim_k)    = full(Fk.zf(20));
    w_oEx_MPC(sim_k)    = full(Fk.zf(29));
    NMPC.w_po_MPC       = [w_po_MPC1;w_po_MPC2;w_oEx_MPC];
    
    w_pg_MPC1(sim_k)    = full(Fk.zf(17));
    w_pg_MPC2(sim_k)    = full(Fk.zf(18));
    w_gEx_MPC(sim_k)    = full(Fk.zf(30));
    NMPC.w_pg_MPC       = [w_pg_MPC1;w_pg_MPC2;w_gEx_MPC];
    
    tgrid = linspace(0,par.T,par.N+1)./60;
    tpast = linspace(-(sim_k-1)*300,0,length(u_in_MPC1))./60;
    

    
    if PlotResult
        
        figure(14)
        clf
        subplot(211)
        hold all
        plot(tpast,w_pg_MPC1,'b')
        plot(tpast,w_pg_MPC2,'r')
        plot(tpast,w_pg_MPC1+w_pg_MPC2,'k')
        plot(tgrid,w_pg1(:,1),'b-.')
        plot(tgrid,w_pg2(:,1),'r-.')
        %         plot(tgrid,w_pg1(:,1)+w_pg2(:,1),'k-.')
        plot(tgrid,w_gEx(:,1),'k-.')
        plot(tgrid,par.QgMax.*ones(size(tgrid)),'k')
        plot([0,0],[0,11],'k:')
        xlim([-60*600,18000]./60)
        grid on
        leg = legend('well 1','well 2','Total Gas rate');
        set(leg,'Interpreter','latex')
        ylabel ('Gas rate [kg/s]','Interpreter','LaTex')
        xlabel ('sample instant N','Interpreter','LaTex')
        %         title ([OptCase ' Optimization'])
        
        subplot(212)
        hold all
        stairs(tpast,u_in_MPC1,'b')
        stairs(tpast,u_in_MPC2,'r')
        plot(tpast,u_in_MPC1+u_in_MPC2,'k')
        stairs(tgrid,u_opt1(:,1),'b-.')
        stairs(tgrid,u_opt2(:,1),'r-.')
        plot(tgrid,u_opt1(:,1)+u_opt2(:,1),'k-.')
        plot([0,0],[0,7],'k:')
        xlim([-60*600,18000]./60)
%         ylim([1,7])
        grid on
        leg = legend('well 1','well 2','Total GL rate');
        set(leg,'Interpreter','latex')
        ylabel ('Gas Lift rate [kg/s]','Interpreter','LaTex')
        xlabel ('sample instant N','Interpreter','LaTex')
        
      
        
    end
    
end
close(h)

%%

DRTO.GOR_true   = GOR_true;
DRTO.par        = par;
DRTO.meas       = meas;
DRTO.NMPC       = NMPC;
DRTO.sim_k      = sim_k;
DRTO.nIter      = nIter;
DRTO.iCost      = iCost;
save('DRTO','DRTO')


%% Plot final results

load('DRTO.mat')

tH = (1:DRTO.sim_k).*DRTO.par.tf./3600; 
meas.tH = (1:DRTO.sim_k*DRTO.par.tf).*DRTO.par.tSim./3600; 

figure(126)
clf
subplot(221)
hold all
plot(meas.tH,DRTO.meas.w_tg,'k')
grid on
ylabel('Gas rate [kg/s]','Interpreter','latex')
xlim([0,5])
xlabel('time [h]','Interpreter','latex')


subplot(222)
hold all
plot(meas.tH,DRTO.meas.w_to,'k')
grid on
ylabel('Oil rate [kg/s]','Interpreter','latex')
xlim([0,5])
xlabel('time [h]','Interpreter','latex')

subplot(223)
stairs(tH,DRTO.NMPC.u_in_MPC(1,:))
hold all
stairs(tH,DRTO.NMPC.u_in_MPC(2,:))
grid on
ylabel('GasLift rate rate [kg/s]','Interpreter','latex')
xlabel('time [h]','Interpreter','latex')
xlim([0,5])

subplot(224)
hold all
plot(tH,DRTO.GOR_true')
grid on
ylabel('GOR [kg/kg]','Interpreter','latex')
xlim([0,5])
xlabel('time [h]','Interpreter','latex')
