function [F,w_gl_SP,Cost,w_opt_SS] = WellSimulator(par)

import casadi.*

n_w     = par.n_w; % no. of wells
R       = par.R;
Mw      = par.Mw;

%% Modelling

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

% Algebraic states
p_ai    = MX.sym('p_ai',n_w);      % 1-2
p_wh    = MX.sym('p_wh',n_w);      % 3-4
p_wi    = MX.sym('p_wi',n_w);      % 5-6
p_bh 	= MX.sym('p_bh',n_w);      % 7-8
rho_ai  = MX.sym('rho_ai',n_w);    % 9-10
rho_m   = MX.sym('rho_m',n_w);     % 11-12
w_iv    = MX.sym('w_iv',n_w);      % 13-14
w_pc    = MX.sym('w_pc',n_w);      % 15-16
w_pg    = MX.sym('w_pg',n_w);      % 17-18
w_po    = MX.sym('w_po',n_w);      % 19-20
w_ro    = MX.sym('w_ro',n_w);      % 21-22
w_rg    = MX.sym('w_rg',n_w);      % 23-24
p_rh    = MX.sym('p_rh',1);        % 25
rho_r   =  MX.sym('rho_r',1);      % 26
p_m     = MX.sym('p_m',1);         % 27
w_pr    = MX.sym('w_pr',1);        % 28
w_to    = MX.sym('w_to',1);        % 29
w_tg    = MX.sym('w_tg',1);        % 30

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

% hold-up gass and oil mass fractions in well and riser
xGwH    = MX.sym('xGwH',n_w);
xOwH    = MX.sym('xOwH',n_w);
xGrH    = MX.sym('xGrH',1);
xOrH    = MX.sym('xOrH',1);
% Flowing mass fraction in well and riser
xGw     = MX.sym('xGw',n_w);
xOw     = MX.sym('xOw',n_w);
xGr     = MX.sym('xGr',1);
xOr     = MX.sym('xOr',1);

slip    = par.slip_real;

xGwH    = (m_gt.*1e3./max(1e-3,(m_gt.*1e3+m_ot.*1e3)));
xOwH    = (m_ot.*1e3./max(1e-3,(m_gt.*1e3+m_ot.*1e3)));
xGrH    = (m_gr.*1e3./(m_gr.*1e3+m_or.*1e3));
xOrH    = (m_or.*1e3./(m_gr.*1e3+m_or.*1e3));

xGw     = slip.*xGwH./(1 + (slip-1).*xGwH);
xOw     = 1 - xGw;
xGr     = slip.*xGrH./(1 + (slip-1).*xGrH);
xOr     = 1 - xGr;

% algebraic equations
f1      = -p_ai.*1e5 + ((R.*T_a./(V_a.*Mw) + 9.81.*H_a./V_a).*m_ga.*1e3) + (Mw./(R.*T_a).*((R.*T_a./(V_a.*Mw) + 9.81.*H_a./V_a).*m_ga.*1e3)).*9.81.*H_a;
f2      = -p_wh.*1e5 + ((R.*T_w./Mw).*(m_gt.*1e3./(L_w.*A_w + L_bh.*A_bh - m_ot.*1e3./rho_o))) - ((m_gt.*1e3+m_ot.*1e3 )./(L_w.*A_w)).*9.81.*H_w/2;
f3      = -p_wi.*1e5 + (p_wh.*1e5 + 9.81./(A_w.*L_w).*max(0,(m_ot.*1e3+m_gt.*1e3-rho_o.*L_bh.*A_bh)).*H_w + 128.*mu_oil.*L_w.*w_pc./(3.14.*D_w.^4.*((m_gt.*1e3 + m_ot.*1e3).*p_wh.*1e5.*Mw.*rho_o)./(m_ot.*1e3.*p_wh.*1e5.*Mw + rho_o.*R.*T_w.*m_gt.*1e3)));
f4      = -p_bh.*1e5 + (p_wi.*1e5 + rho_o.*9.81.*H_bh + 128.*mu_oil.*L_bh.*w_ro./(3.14.*D_bh.^4.*rho_o));
f5      = -rho_ai.*1e2 +(Mw./(R.*T_a).*p_ai.*1e5);
f6      = -rho_m.*1e2 + ((m_gt.*1e3 + m_ot.*1e3).*p_wh.*1e5.*Mw.*rho_o)./(m_ot.*1e3.*p_wh.*1e5.*Mw + rho_o.*R.*T_w.*m_gt.*1e3);
f7      = -w_iv + C_iv.*sqrt(rho_ai.*1e2.*(p_ai.*1e5 - p_wi.*1e5));
f8      = -w_pc + 1.*C_pc.*sqrt(rho_m.*1e2.*(p_wh.*1e5 - p_m.*1e5));
f9      = -w_pg + xGw.*w_pc;
f10     = -w_po + xOw.*w_pc;
f11     = -w_ro + PI.*1e-6.*(p_res.*1e5 - p_bh.*1e5);
f12     = -w_rg.*1e-1 + GOR.*w_ro;  % 23 - 24
f13     = -p_rh.*1e5 + ((R.*T_r./Mw).*(m_gr.*1e3./(L_r.*A_r))) - ((m_gr.*1e3+m_or.*1e3 )./(L_r.*A_r)).*9.81.*H_r/2; %25
f14     = -rho_r.*1e2 + ((m_gr.*1e3 + m_or.*1e3).*p_rh.*1e5.*Mw.*rho_ro)./(m_or.*1e3.*p_rh.*1e5.*Mw + rho_ro.*R.*T_r.*m_gr.*1e3); % 26
f15     = -p_m.*1e5 + (p_rh.*1e5 + 9.81./(A_r.*L_r).*(m_or.*1e3+m_gr.*1e3).*H_r + 128.*mu_oil.*L_r.*w_pr./(3.14.*D_r.^4.*((m_gr.*1e3 + m_or.*1e3).*p_rh.*1e5.*Mw.*rho_ro)./(m_or.*1e3.*p_rh.*1e5.*Mw + rho_ro.*R.*T_r.*m_gr.*1e3)));%27
f16     = -w_pr + 1.*C_pr.*sqrt(rho_r.*1e2.*(p_rh.*1e5 - p_s.*1e5)); %28
f17     = -w_to + xOr.*w_pr; % 29
f18     = -w_tg + xGr.*w_pr; %30

% differential equations
df1     = (w_gl - w_iv).*1e-3;
df2     = (w_iv + w_rg.*1e-1 - w_pg).*1e-3;
df3     = (w_ro - w_po).*1e-3;
df4     = (sum(w_pg) - w_tg).*1e-3 ;
df5     = (sum(w_po) - w_to).*1e-3 ;

% Form the DAE system
diff    = vertcat(df1,df2,df3,df4,df5);
alg     = vertcat(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18);

% give parameter values
alg     = substitute(alg,p_res,par.p_res);
alg     = substitute(alg,PI,par.PI);
alg     = substitute(alg,p_s,par.p_s);
alg     = substitute(alg,T_a,par.T_a);
alg     = substitute(alg,T_w,par.T_w);
alg     = substitute(alg,T_r,par.T_r);

% concatenate the differential and algebraic states
x_var   = vertcat(m_ga,m_gt,m_ot,m_gr,m_or);
z_var   = vertcat(p_ai,p_wh,p_wi,p_bh,rho_ai,rho_m,w_iv,w_pc,w_pg,w_po,w_ro,w_rg,p_rh,rho_r,p_m,w_pr,w_to,w_tg);
p_var   = vertcat(w_gl,GOR);

% stage cost
%L       = -(w_to.^2) + 0.5.*sum((w_gl).^2);
L       = -w_to + sum(w_gl);  %EKA

% Build DAE system
dae     = struct('x',x_var,'z',z_var,'p',p_var,'ode',diff,'alg',alg,'quad',L);  % der(m_tot) = w_in - w_out;
opts    = struct('tf',par.tSim);

% create IDAS integrator
F       = integrator('F','idas',dae,opts);


%% Steady state optimization
[dx0,z0,u0,lbx,lbz,lbu,ubx,ubz,ubu] = GasLiftRiser_Initialization_bounds(par);

% decision variables
w   = {};
w0  = [];
lbw = [];
ubw = [];

% constraints
g   = {};
lbg = [];
ubg = [];

w   = {w{:},x_var,z_var,p_var};
lbw = [lbw;lbx;lbz;lbu;par.GOR];
ubw = [ubw;ubx;ubz;ubu;par.GOR];
w0  = [w0;dx0;z0;u0;par.GOR];

J = L; % Economic objective

%Add the system model as constraints
g = {g{:},vertcat(diff,alg)};
lbg = [lbg;zeros(8+30,1)];
ubg = [ubg;zeros(8+30,1)];

% Add gas lift and gas capacity constraints
% g = {g{:},sum(p_var(1:2)),sum(z_var(17:18))};
% lbg = [lbg;0;0];
% ubg = [ubg;par.qGLMax;par.QgMax];
%g = {g{:},sum(p_var(1:2)),(sum(z_var(17:18))-par.QgMax)};
% g   = {g{:},sum(p_var(1:2)),sum(z_var(17:18))- z_var(30)};
% lbg = [lbg;0;0];
% ubg = [ubg;par.qGLMax;0];
% g   = {g{:},sum(z_var(17:18))};
% lbg = [lbg;0];
% ubg = [ubg;par.QgMax];

% formalize it into an NLP problem
nlp = struct('x',vertcat(w{:}),'f',J,'g',vertcat(g{:}));

% Assign solver
solver = nlpsol('solver','ipopt',nlp);

% Solve
sol = solver('x0',w0,'lbx',lbw,'ubx',ubw,'lbg',lbg,'ubg',ubg);

% Extract Solution
w_opt_SS = full(sol.x);
w_gl_SP = w_opt_SS(end-3:end-2);
Cost = full(sol.f);

% % Add Greshgorin Bound
% ceq  = vertcat(g{:});
% lamda  = full(sol.lam_g);
% numX   = size(w_opt_SS,1);
% numY   = size(lamda,1);
% lambda.eqnonlin = lamda;
% l = {};
% for i=1:numY       
%    l{i} = MX.sym(['l_' num2str(i)],1);
% end
% 
% l         = vertcat(l{:});
% w         = vertcat(w{:});
% Lagrange  = J + l'*ceq;
% Lagr      = Function('Lagr', {w,l}, {Lagrange}, char('w','l'), char('Lagr'));
% H         = Function(Lagr.hessian('w','Lagr'));
% cons      = Function('Const', {w}, {ceq}, char('w'), char('cons'));
% Jcon      = Function(cons.jacobian('w','cons'));
% 
% eqVal = cons(w_opt_SS);
% Hx   = H(w_opt_SS,lambda.eqnonlin);
% Hx   = full(Hx);
% Jac  = Jcon(w_opt_SS);
% Jac  = full(Jac);
% % display nullspace of the constraint and its eigenvalue
% rH   = null(Jac)'*Hx*null(Jac);
% erH  = eig(rH);
% 
% [Hxxl,QmaxGL]   = Greshgorin(Hx);
% save QmaxGL.mat QmaxGL;
end

function [H,Q] = Greshgorin(H)
numH    = size(H,1);
Q       = zeros(numH,numH);
delta   = 1e-2;      % normal case
%delta   = 2.5;      % with measurement noise 1 percent
for i=1:numH  % iterate all row of Hessian
    sumRow = 0;
    for j=1:numH
        if j ~= i
            sumRow = sumRow + abs(H(i,j));
        end
    end
    
    if H(i,i) <= sumRow   % include equality 
        Q(i,i) = sumRow - H(i,i) + delta;
    end
end
Q = diag(Q);
end