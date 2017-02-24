%OPTDISTCOLACSTR Summary of this function goes here
% 
% [OUTPUTARGS] = OPTDISTCOLACSTR(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2016/03/14 10:12:15 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2016

import casadi.* 
format long;

%% the model
% model parameters
NT = 41;
Uf = 1.2;           % Feeding rate F

% invoke the model
%[t,state,xdot,inputs] = DistColA(Uf);
[t,state,xdot,inputs] = DistColACstr(Uf);

f = Function('f',{state,inputs}, {xdot});

% bound constraints
xB_max = 0.008;
xD_min = 0.9500;
%xD_min = 0.9490;
%xD_min = 0.900;
V_max  = 4.008;

% State bounds and initial guess
x_min     = zeros(84,1);
x_max     = ones(84,1);
x_max(1)  = xB_max;
x_min(41) = xD_min;

% Control bounds
u_min = [0.1; 0.1];
u_max = [10; V_max];

%% Steady-state optimization

load u_opt_cstr.mat %result from FMINCON
xdot_val_rf_ss = x_opt_ss;

%load cstr_init.mat;
%x0 = Xinit;
Xinit = x_opt_ss;

% prices
pf = 1; 
pV = 0.01; 
pB = 1; 
pD = 2;
F = Uf;

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
nx = 84;
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
%w0  = [w0; xdot_val_rf_ss];   % use initial guess from steady-state optimization results. this should give zero objective function (for tracking term!)
w0  = [w0; Xinit];
%g   = {g{:}, X0 - xdot_val_rf_ss};
g   = {g{:}, X0 - Xinit};
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
    w0  = [w0;  u_opt_ss];       % optimized results from steady-state optimization

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
    %MB = Xk(NT+1);  
    %MD = Xk(2*NT);              % Actual reboiler and condenser holdup
    MB = Xk(NT+2);
    MD = Xk(2*NT+1);
    Do = Ds+(MD-MDs)*KcD;       % Distillate flow
    Bo = Bs+(MB-MBs)*KcB;       % Bottoms flow
    econ_term    = (pf*F + pV*Uk(2) - pB*Bo - pD*Do)*delta_time;             % economic term                 
    control_term = (Uk - u_opt_ss)'*(Uk - u_opt_ss)*delta_time;                        % tracking: control input term
    state_term   = (Xk - xdot_val_rf_ss)'*(Xk - xdot_val_rf_ss)*delta_time;  % tracking: states term
    %J            = J + alpha*econ_term + beta*control_term + gamma*state_term;
    J            = J + state_term;
    %J            = J + control_term;
    %J            = J + econ_term;
end

% Create an NLP solver
prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));
solver = nlpsol('solver', 'ipopt', prob);

% Solve the NLP
sol   = solver('x0', w0, 'lbx', lbw, 'ubx', ubw, 'lbg', lbg, 'ubg', ubg);
%sol   = solver('x0', w0, 'lbg', lbg, 'ubg', ubg);  % without bound constraint... drive system immediately to steady-state equilibrium! (with combination objective function, IPOPT doesn't converge quickly!)

v_opt     = full(sol.x);
lamda_opt = full(sol.lam_g);

%keyboard


% %===============================================================================
% % SAVE RESULTS FOR PATH-FOLLOWING ALGORITHM INPUTS
% %===============================================================================
% X_init, Lamda AND upper bound and lower bound !
% x_opt=v_opt; y_opt=lamda_opt; lb = vars_lb; ub = vars_ub;
x_opt=v_opt; y_opt=lamda_opt; lb = lbw; ub = ubw;
save nlp0_cstr.mat x_opt y_opt lb ub;

%===============================================================================
% PLOTTING
%===============================================================================
% Get values at the beginning of each finite element
x1_opt  = v_opt(1:(nx+nu)+nx*d:end);
x2_opt  = v_opt(2:(nx+nu)+nx*d:end);
x3_opt  = v_opt(3:(nx+nu)+nx*d:end);
x4_opt  = v_opt(4:(nx+nu)+nx*d:end);
x5_opt  = v_opt(5:(nx+nu)+nx*d:end);
x6_opt  = v_opt(6:(nx+nu)+nx*d:end);
x7_opt  = v_opt(7:(nx+nu)+nx*d:end);
x8_opt  = v_opt(8:(nx+nu)+nx*d:end);
x9_opt  = v_opt(9:(nx+nu)+nx*d:end);
x10_opt = v_opt(10:(nx+nu)+nx*d:end);
x11_opt = v_opt(11:(nx+nu)+nx*d:end);
x12_opt = v_opt(12:(nx+nu)+nx*d:end);
x13_opt = v_opt(13:(nx+nu)+nx*d:end);
x14_opt = v_opt(14:(nx+nu)+nx*d:end);
x15_opt = v_opt(15:(nx+nu)+nx*d:end);
x16_opt = v_opt(16:(nx+nu)+nx*d:end);
x17_opt = v_opt(17:(nx+nu)+nx*d:end);
x18_opt = v_opt(18:(nx+nu)+nx*d:end);
x19_opt = v_opt(19:(nx+nu)+nx*d:end);
x20_opt = v_opt(20:(nx+nu)+nx*d:end);
x21_opt = v_opt(21:(nx+nu)+nx*d:end);
x22_opt = v_opt(22:(nx+nu)+nx*d:end);
x23_opt = v_opt(23:(nx+nu)+nx*d:end);
x24_opt = v_opt(24:(nx+nu)+nx*d:end);
x25_opt = v_opt(25:(nx+nu)+nx*d:end);
x26_opt = v_opt(26:(nx+nu)+nx*d:end);
x27_opt = v_opt(27:(nx+nu)+nx*d:end);
x28_opt = v_opt(28:(nx+nu)+nx*d:end);
x29_opt = v_opt(29:(nx+nu)+nx*d:end);
x30_opt = v_opt(30:(nx+nu)+nx*d:end);
x31_opt = v_opt(31:(nx+nu)+nx*d:end);
x32_opt = v_opt(32:(nx+nu)+nx*d:end);
x33_opt = v_opt(33:(nx+nu)+nx*d:end);
x34_opt = v_opt(34:(nx+nu)+nx*d:end);
x35_opt = v_opt(35:(nx+nu)+nx*d:end);
x36_opt = v_opt(36:(nx+nu)+nx*d:end);
x37_opt = v_opt(37:(nx+nu)+nx*d:end);
x38_opt = v_opt(38:(nx+nu)+nx*d:end);
x39_opt = v_opt(39:(nx+nu)+nx*d:end);
x40_opt = v_opt(40:(nx+nu)+nx*d:end);
x41_opt = v_opt(41:(nx+nu)+nx*d:end);
x42_opt = v_opt(42:(nx+nu)+nx*d:end);
x43_opt = v_opt(43:(nx+nu)+nx*d:end);
x44_opt = v_opt(44:(nx+nu)+nx*d:end);
x45_opt = v_opt(45:(nx+nu)+nx*d:end);
x46_opt = v_opt(46:(nx+nu)+nx*d:end);
x47_opt = v_opt(47:(nx+nu)+nx*d:end);
x48_opt = v_opt(48:(nx+nu)+nx*d:end);
x49_opt = v_opt(49:(nx+nu)+nx*d:end);
x50_opt = v_opt(50:(nx+nu)+nx*d:end);
x51_opt = v_opt(51:(nx+nu)+nx*d:end);
x52_opt = v_opt(52:(nx+nu)+nx*d:end);
x53_opt = v_opt(53:(nx+nu)+nx*d:end);
x54_opt = v_opt(54:(nx+nu)+nx*d:end);
x55_opt = v_opt(55:(nx+nu)+nx*d:end);
x56_opt = v_opt(56:(nx+nu)+nx*d:end);
x57_opt = v_opt(57:(nx+nu)+nx*d:end);
x58_opt = v_opt(58:(nx+nu)+nx*d:end);
x59_opt = v_opt(59:(nx+nu)+nx*d:end);
x60_opt = v_opt(60:(nx+nu)+nx*d:end);
x61_opt = v_opt(61:(nx+nu)+nx*d:end);
x62_opt = v_opt(62:(nx+nu)+nx*d:end);
x63_opt = v_opt(63:(nx+nu)+nx*d:end);
x64_opt = v_opt(64:(nx+nu)+nx*d:end);
x65_opt = v_opt(65:(nx+nu)+nx*d:end);
x66_opt = v_opt(66:(nx+nu)+nx*d:end);
x67_opt = v_opt(67:(nx+nu)+nx*d:end);
x68_opt = v_opt(68:(nx+nu)+nx*d:end);
x69_opt = v_opt(69:(nx+nu)+nx*d:end);
x70_opt = v_opt(70:(nx+nu)+nx*d:end);
x71_opt = v_opt(71:(nx+nu)+nx*d:end);
x72_opt = v_opt(72:(nx+nu)+nx*d:end);
x73_opt = v_opt(73:(nx+nu)+nx*d:end);
x74_opt = v_opt(74:(nx+nu)+nx*d:end);
x75_opt = v_opt(75:(nx+nu)+nx*d:end);
x76_opt = v_opt(76:(nx+nu)+nx*d:end);
x77_opt = v_opt(77:(nx+nu)+nx*d:end);
x78_opt = v_opt(78:(nx+nu)+nx*d:end);
x79_opt = v_opt(79:(nx+nu)+nx*d:end);
x80_opt = v_opt(80:(nx+nu)+nx*d:end);
x81_opt = v_opt(81:(nx+nu)+nx*d:end);
x82_opt = v_opt(82:(nx+nu)+nx*d:end);
x83_opt = v_opt(83:(nx+nu)+nx*d:end);
x84_opt = v_opt(84:(nx+nu)+nx*d:end);


%u1_opt  = v_opt((d+1)*nx:(d+1)*nx+nu:NV);
%u2_opt  = v_opt((d+1)*nx+1:(d+1)*nx+nu:NV);
u1_opt  = v_opt(85:(nx+nu)+nx*d:end);
u2_opt  = v_opt(86:(nx+nu)+nx*d:end);
tgrid   = linspace(0,tf,nk+1);
tgrid_u = linspace(0,tf,nk);


figure(1)
clf()
plot(tgrid,x1_opt, 'LineWidth',2.5); hold on;
plot(tgrid,x2_opt, 'LineWidth',2.5); hold on;
plot(tgrid,x21_opt, 'LineWidth',2.5); hold on;
plot(tgrid,x40_opt, 'LineWidth',2.5); hold on;
plot(tgrid,x41_opt, 'LineWidth',2.5); hold on;
plot(tgrid,x42_opt, 'LineWidth',2.5);
legend('x[1]-reboiler','x[2]','x[21]-feeder','x[40]', 'x[41]-condenser', 'x[42]-CSTR');
xlabel('time [minutes]');
ylabel('state variables');

figure(2)
plot(tgrid_u,u1_opt, 'LineWidth',2.5); hold on;
plot(tgrid_u,u2_opt, 'LineWidth',2.5);
legend('u[1]-LT','u[2]-VB');
ylabel('control inputs');
xlabel('time [minutes]');

% figure(3)
% plot(tgrid,x1_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x2_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x3_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x4_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x5_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x6_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x7_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x8_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x9_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x10_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x11_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x12_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x13_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x14_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x15_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x16_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x17_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x18_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x19_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x20_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x21_opt, 'LineWidth',2.5); 
% legend('x[1]','x[2]','x[3]','x[4]', 'x[5]', 'x[6]','x[7]','x[8]','x[9]', 'x[10]', 'x[11]','x[12]','x[13]','x[14]', 'x[15]', 'x[16]','x[17]','x[18]','x[19]', 'x[20]', 'x[21]');
% xlabel('time [minutes]');
% ylabel('state variables');
% title('Distillation column A states: x[1] - x[21]');
% 
% figure(4)
% plot(tgrid,x22_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x23_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x24_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x25_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x26_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x27_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x28_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x29_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x30_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x31_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x32_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x33_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x34_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x35_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x36_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x37_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x38_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x39_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x40_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x41_opt, 'LineWidth',2.5); 
% legend('x[22]','x[23]','x[24]', 'x[25]', 'x[26]','x[27]','x[28]','x[29]', 'x[30]', 'x[31]','x[32]','x[33]','x[34]', 'x[35]', 'x[36]','x[37]','x[38]','x[39]', 'x[40]', 'x[41]');
% xlabel('time [minutes]');
% ylabel('state variables');
% title('Distillation column A states: x[22] - x[41]');
% 
% figure(5)
% plot(tgrid,x42_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x43_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x44_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x45_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x46_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x47_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x48_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x49_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x50_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x51_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x52_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x53_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x54_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x55_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x56_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x57_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x58_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x59_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x60_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x61_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x62_opt, 'LineWidth',2.5); 
% legend('x[42]','x[43]','x[44]', 'x[45]', 'x[46]','x[47]','x[48]','x[49]', 'x[50]', 'x[51]','x[52]','x[53]','x[54]', 'x[55]', 'x[56]','x[57]','x[58]','x[59]', 'x[60]', 'x[61]', 'x[62]');
% xlabel('time [minutes]');
% ylabel('state variables');
% title('Distillation column A states: x[42] - x[62]');
% 
% figure(6)
% plot(tgrid,x63_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x64_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x65_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x66_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x67_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x68_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x69_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x70_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x71_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x72_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x73_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x74_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x75_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x76_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x77_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x78_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x79_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x80_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x81_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x82_opt, 'LineWidth',2.5);
% legend('x[63]','x[64]', 'x[65]', 'x[66]','x[67]','x[68]','x[69]', 'x[70]', 'x[71]','x[72]','x[73]','x[74]', 'x[75]', 'x[76]','x[77]','x[78]','x[79]', 'x[80]', 'x[81]', 'x[82]');
% xlabel('time [minutes]');
% ylabel('state variables');
% title('Distillation column A states: x[63] - x[82]');

x_sol    = v_opt(end-83:end);
% x_sol    = v_opt(end-40:end);
%diff_sol = xf - x_sol;
diff_sol = xdot_val_rf_ss - x_sol;
err_sol  = norm(diff_sol,inf);
fprintf('The states error: %f \n', err_sol);

