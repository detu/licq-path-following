function [solver,MPC] = SetpointNMPC(par,MPCinit)
addpath ('C:\Users\dineshk\CasADi\casadi-matlabR2014b-v3.1.0-rc1')
import casadi.*


[~,~,~,lbx,lbz,lbu,ubx,ubz,ubu] = GasLiftRiser_Initialization_bounds(par);

dx0 = MPCinit.dx0;
z0 = MPCinit.z0;
u0 = MPCinit.u0;
u_in = MPCinit.u_in;

GOR_val = par.GOR;
PI_val = par.PI;
d0 = zeros(8,1);

n_w = par.n_w;
R = par.R;
Mw = par.Mw;



%% Modelling

L_w = par.L_w;
H_w = par.H_w;
D_w = par.D_w;

L_bh = par.L_bh;
H_bh = par.H_bh;
D_bh = par.D_bh;

L_a = par.L_a;
H_a = par.H_a;
D_a = par.D_a;

L_r = par.L_r;
H_r = par.H_r;
D_r = par.D_r;

rho_o = par.rho_o;
C_iv = par.C_iv;
C_pc = par.C_pc;
C_pr = par.C_pr;
rho_ro = sum(rho_o)/2;

mu_oil = 1*0.001; % 1cP oil viscosity

A_w = pi.*(D_w/2).^2;
A_bh = pi.*(D_bh/2).^2;
V_a = L_a.*(pi.*(D_a/2).^2 - pi.*(D_w/2).^2);
A_r = pi.*(D_r/2).^2;

% differential states
m_ga = MX.sym('m_ga',n_w); % 1-2
m_gt = MX.sym('m_gt',n_w); % 3-4
m_ot = MX.sym('m_ot',n_w); % 5-6
m_gr = MX.sym('m_gr',1);   % 7
m_or = MX.sym('m_or',1);   % 8

% Algebraic states
p_ai = MX.sym('p_ai',n_w);      % 1-2
p_wh = MX.sym('p_wh',n_w);      % 3-4
p_wi = MX.sym('p_wi',n_w);      % 5-6
p_bh = MX.sym('p_bh',n_w);      % 7-8
rho_ai = MX.sym('rho_ai',n_w);  % 9-10
rho_m = MX.sym('rho_m',n_w);    % 11-12
w_iv = MX.sym('w_iv',n_w);      % 13-14
w_pc = MX.sym('w_pc',n_w);      % 15-16
w_pg = MX.sym('w_pg',n_w);      % 17-18
w_po = MX.sym('w_po',n_w);      % 19-20
w_ro = MX.sym('w_ro',n_w);      % 21-22
w_rg = MX.sym('w_rg',n_w);      % 23-24
p_rh = MX.sym('p_rh',1);        % 25
rho_r =  MX.sym('rho_r',1);     % 26
p_m = MX.sym('p_m',1);          % 27
w_pr = MX.sym('w_pr',1);        % 28
w_to = MX.sym('w_to',1);        % 29
w_tg = MX.sym('w_tg',1);        % 30

% control input
w_gl = MX.sym('w_gl',n_w);
w_gl_SP = MX.sym('w_gl_SP',n_w);

% parameters
p_res = MX.sym('p_res',n_w);
PI = MX.sym('PI',n_w);
GOR = MX.sym('GOR',n_w);
T_a = MX.sym('T_a',n_w);
T_w = MX.sym('T_w',n_w);
T_r = MX.sym('T_r',1);
p_s = MX.sym('p_s',1);

% hold-up gass and oil mass fractions in well and riser
xGwH = MX.sym('xGwH',n_w);
xOwH = MX.sym('xOwH',n_w);
xGrH = MX.sym('xGrH',1);
xOrH = MX.sym('xOrH',1);
% Flowing mass fraction in well and riser
xGw = MX.sym('xGw',n_w);
xOw = MX.sym('xOw',n_w);
xGr = MX.sym('xGr',1);
xOr = MX.sym('xOr',1);

slip = 1.0;

xGwH = (m_gt.*1e3./max(1e-3,(m_gt.*1e3+m_ot.*1e3)));
xOwH = (m_ot.*1e3./max(1e-3,(m_gt.*1e3+m_ot.*1e3)));
xGrH = (m_gr.*1e3./(m_gr.*1e3+m_or.*1e3));
xOrH = (m_or.*1e3./(m_gr.*1e3+m_or.*1e3));

xGw = slip.*xGwH./(1 + (slip-1).*xGwH);
xOw = 1 - xGw;
xGr = slip.*xGrH./(1 + (slip-1).*xGrH);
xOr = 1 - xGr;

d1 = MX.sym('d1',n_w); % 1-2
d2 = MX.sym('d2',n_w); % 3-4
d3 = MX.sym('d3',n_w); % 5-6
d4 = MX.sym('d4',1);   % 7
d5 = MX.sym('d5',1);   % 8

% algebraic equations
f1 = -p_ai.*1e5 + ((R.*T_a./(V_a.*Mw) + 9.81.*H_a./V_a).*m_ga.*1e3) + (Mw./(R.*T_a).*((R.*T_a./(V_a.*Mw) + 9.81.*H_a./V_a).*m_ga.*1e3)).*9.81.*H_a;
f2 = -p_wh.*1e5 + ((R.*T_w./Mw).*(m_gt.*1e3./(L_w.*A_w + L_bh.*A_bh - m_ot.*1e3./rho_o))) - ((m_gt.*1e3+m_ot.*1e3 )./(L_w.*A_w)).*9.81.*H_w/2;
f3 = -p_wi.*1e5 + (p_wh.*1e5 + 9.81./(A_w.*L_w).*max(0,(m_ot.*1e3+m_gt.*1e3-rho_o.*L_bh.*A_bh)).*H_w + 128.*mu_oil.*L_w.*w_pc./(3.14.*D_w.^4.*((m_gt.*1e3 + m_ot.*1e3).*p_wh.*1e5.*Mw.*rho_o)./(m_ot.*1e3.*p_wh.*1e5.*Mw + rho_o.*R.*T_w.*m_gt.*1e3)));
f4 = -p_bh.*1e5 + (p_wi.*1e5 + rho_o.*9.81.*H_bh + 128.*mu_oil.*L_bh.*w_ro./(3.14.*D_bh.^4.*rho_o));
f5 = -rho_ai.*1e2 +(Mw./(R.*T_a).*p_ai.*1e5);
% f6 = -rho_m.*1e2 + (m_gt.*1e3+m_ot.*1e3 )./(L_w.*A_w);
f6 = -rho_m.*1e2 + ((m_gt.*1e3 + m_ot.*1e3).*p_wh.*1e5.*Mw.*rho_o)./(m_ot.*1e3.*p_wh.*1e5.*Mw + rho_o.*R.*T_w.*m_gt.*1e3);
f7 = -w_iv + C_iv.*sqrt(rho_ai.*1e2.*(p_ai.*1e5 - p_wi.*1e5));
f8 = -w_pc + 1.*C_pc.*sqrt(rho_m.*1e2.*(p_wh.*1e5 - p_m.*1e5));
f9 = -w_pg + xGw.*w_pc;
f10 = -w_po + xOw.*w_pc;
f11 = -w_ro + par.PI.*1e-6.*(p_res.*1e5 - p_bh.*1e5);
f12 = -w_rg.*1e-1 + GOR.*w_ro;  % 23 - 24

f13 = -p_rh.*1e5 + ((R.*T_r./Mw).*(m_gr.*1e3./(L_r.*A_r))) - ((m_gr.*1e3+m_or.*1e3 )./(L_r.*A_r)).*9.81.*H_r/2; %25
% f14 = -rho_r.*1e2 + (m_gr.*1e3+m_or.*1e3 )./(L_r.*A_r); % 26
f14 = -rho_r.*1e2 + ((m_gr.*1e3 + m_or.*1e3).*p_rh.*1e5.*Mw.*rho_ro)./(m_or.*1e3.*p_rh.*1e5.*Mw + rho_ro.*R.*T_r.*m_gr.*1e3);
f15 = -p_m.*1e5 + (p_rh.*1e5 + 9.81./(A_r.*L_r).*(m_or.*1e3+m_gr.*1e3).*H_r + 128.*mu_oil.*L_r.*w_pr./(3.14.*D_r.^4.*((m_gr.*1e3 + m_or.*1e3).*p_rh.*1e5.*Mw.*rho_ro)./(m_or.*1e3.*p_rh.*1e5.*Mw + rho_ro.*R.*T_r.*m_gr.*1e3)));%27
f16 = -w_pr + 1.*C_pr.*sqrt(rho_r.*1e2.*(p_rh.*1e5 - p_s.*1e5)); %28
f17 = -w_to + xOr.*w_pr; % 29
f18 = -w_tg + xGr.*w_pr; %30

% differential equations
df1 = (w_gl - w_iv).*1e-3 + d1;
df2 = (w_iv + w_rg.*1e-1 - w_pg).*1e-3 + d2;
df3 = (w_ro - w_po).*1e-3 + d3;
df4 = (sum(w_pg) - w_tg).*1e-3 + d4 ;
df5 = (sum(w_po) - w_to).*1e-3 + d5 ;

% Form the DAE system
diff = vertcat(df1,df2,df3,df4,df5);
alg = vertcat(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18);

% give parameter values
alg = substitute(alg,p_res,par.p_res);
alg = substitute(alg,p_s,par.p_s);
alg = substitute(alg,T_a,par.T_a);
alg = substitute(alg,T_w,par.T_w);
alg = substitute(alg,T_r,par.T_r);

% concatenate the differential and algebraic states
x_var = vertcat(m_ga,m_gt,m_ot,m_gr,m_or);
z_var = vertcat(p_ai,p_wh,p_wi,p_bh,rho_ai,rho_m,w_iv,w_pc,w_pg,w_po,w_ro,w_rg,p_rh,rho_r,p_m,w_pr,w_to,w_tg);
p_var = vertcat(w_gl,w_gl_SP,GOR,PI,d1,d2,d3,d4,d5);

L = sum((w_gl-w_gl_SP).^2); % Penalize setpoint deviation

f = Function('f',{x_var,z_var,p_var},{diff,alg,L},{'x','z','p'},{'xdot','zeval','qj'});

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




nu = 2;%length(p_var);
nz = 30;
nx = 8;
np = nu;
nd = nx;

%% Build NLP solver

% empty nlp
w = {};
w0 = [];
lbw = [];
ubw = [];
J = 0;

g = {};
lbg = [];
ubg = [];

% initial conditions for each scenario
X0 = MX.sym('X0',nx);
Z0 = MX.sym('Z0',nz);
w = {w{:}, X0,Z0};
lbw = [lbw; dx0;z0];
ubw = [ubw; dx0;z0];
w0 = [w0; dx0;z0];

U0 = MX.sym('U0',2);
w = {w{:}, U0};
lbw = [lbw; u0];
ubw = [ubw; u0];
w0 = [w0; u0];

% Formulate NLP
Xk = X0;
Xkj = {};
Zkj = {};
Uk_prev = U0;
js = 1;
for k = 0:par.N-1
    
    Uk = MX.sym(['U_' num2str(k) '_' num2str(js)],nu);
    Usp_k = MX.sym(['Usp_' num2str(k) '_' num2str(js)],nu);
    GOR_k = MX.sym(['GOR_' num2str(k) '_' num2str(js)],nu);
    PI_k  = MX.sym(['PI_' num2str(k) '_' num2str(js)],nu);
    Udk   = MX.sym(['Ud_' num2str(k) '_' num2str(js)],nx);
    Upar = vertcat(Uk,Usp_k,GOR_k,PI_k,Udk);
    
    w = {w{:},Upar};
    lbw = [lbw;lbu;u_in;GOR_val;PI_val;d0];
    ubw = [ubw;ubu;u_in;GOR_val;PI_val;d0];
    w0 = [w0;u0;u_in;GOR_val;PI_val;d0];
    
    Xkj = {};
    Zkj = {};
    
    for j = 1:d
        Xkj{j} = MX.sym(['X_' num2str(k) '_' num2str(j) '_' num2str(js)],nx);
        Zkj{j} = MX.sym(['Z_' num2str(k) '_' num2str(j) '_' num2str(js)],nz);
        if par.sv
            s{j} = MX.sym(['s_' num2str(k) '_' num2str(j) '_' num2str(js)],1);
        end
        w = {w{:},Xkj{j},Zkj{j}};
        lbw = [lbw;lbx;lbz];
        ubw = [ubw;ubx;ubz];
        w0 = [w0; dx0;z0];
    end
    
    % Loop over collocation points
    Xk_end = D(1)*Xk;
    
    for j = 1:d
        % Expression for the state derivative at the collocation point
        xp = C(1,j+1)*Xk;  % helper state
        for r = 1:d
            xp = xp + C(r+1,j+1)*Xkj{r};
        end
        [fj,zj,qj] =  f(Xkj{j},Zkj{j},Upar);
        
        g = {g{:},par.tf*fj-xp,zj};  % dynamics and algebraic constraints
        lbg = [lbg;zeros(nx,1);zeros(nz,1)];
        ubg = [ubg;zeros(nx,1);zeros(nz,1)];
        
        % Gas capacity constraints on all the collocation points
        if par.sv
            g = {g{:},Zkj{j}(30)-s{j}};  %
        else
            g = {g{:},Zkj{j}(30)};
        end
        lbg = [lbg;0];
        ubg = [ubg;par.QgMax];
        
        % Slack variables for gas capacity constraints
        if par.sv
            w = {w{:},s{j}};
            lbw = [lbw;0];
            ubw = [ubw;2];
            w0 = [w0;0];
        end
        
        % Add contribution to the end states
        Xk_end = Xk_end + D(j+1)*Xkj{j};
        
        if par.sv
            J = J + (B(j+1)*qj*par.tf ) + 100*s{j} + 10.*sum((Uk_prev - Uk).^2) ;
        else
            J = J + (B(j+1)*qj*par.tf ) + 10.*sum((Uk_prev - Uk).^2) ; %
        end
    end
    
    Uk_prev = MX.sym(['Uprev_' num2str(k+1)],nu);
    Uk_prev = Uk;
    
    % New NLP variable for state at end of interval
    Xk = MX.sym(['X_' num2str(k+1) '_' num2str(js)], nx);
    w = {w{:},Xk};
    lbw = [lbw;lbx];
    ubw = [ubw;ubx];
    w0 = [w0; dx0];
    
    % Shooting Gap constraint
    g = {g{:},Xk_end-Xk};
    lbg = [lbg;zeros(nx,1)];
    ubg = [ubg;zeros(nx,1)];
    
    g = {g{:},sum(Uk)};
    lbg = [lbg;0];
    ubg = [ubg;par.qGLMax];
    
    
end

% create and solve NLP solver

nlp = struct('x',vertcat(w{:}),'f',J,'g',vertcat(g{:}));

solver = nlpsol('solver','ipopt',nlp);

MPC.lbg = lbg;
MPC.ubg = ubg;
