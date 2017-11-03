function [f_EKF,JacFx,h_EKF,JacHx,z_EKF,yIndex,nxEKF,JacAx,JacBu,JacJx,JacJu,J_EKF] = EKF_integrating(par)

% Import CasADi
%addpath ('C:\Users\dineshk\CasADi\casadi-matlabR2014b-v3.1.0-rc1')
import casadi.*

n_w = par.n_w; % no. of wells;
Mw      = par.Mw;
R       = par.R;

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
m_ga = MX.sym('m_ga',n_w); % 1-2
m_gt = MX.sym('m_gt',n_w); % 3-4
m_ot = MX.sym('m_ot',n_w); % 5-6
m_gr = MX.sym('m_gr',1);   % 7
m_or = MX.sym('m_or',1);   % 8

% GOR = MX.sym('GOR',n_w);   %9-10
% wp1 = MX.sym('wp1',n_w); % small artifical Noise for estimating GOR

% control input
w_gl = MX.sym('w_gl',n_w);

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
w_pg    = (m_gt.*1e3./(m_gt.*1e3+m_ot.*1e3)).*w_pc;
w_po    = (m_ot.*1e3./(m_gt.*1e3+m_ot.*1e3)).*w_pc;
p_wi    = 1e-5.*((p_wh.*1e5 + 9.81./(A_w.*L_w).*max(0,(m_ot.*1e3+m_gt.*1e3-rho_o.*L_bh.*A_bh)).*H_w + 128.*mu_oil.*L_w.*w_pc./(3.14.*D_w.^4.*((m_gt.*1e3 + m_ot.*1e3).*p_wh.*1e5.*Mw.*rho_o)./(m_ot.*1e3.*p_wh.*1e5.*Mw + rho_o.*R.*par.T_w.*m_gt.*1e3))));
p_bh    = 1e-5.*(p_wi.*1e5 + rho_o.*9.81.*H_bh + 128.*mu_oil.*L_bh.*w_po./(3.14.*D_bh.^4.*rho_o));
w_iv    = C_iv.*sqrt(rho_ai.*1e2.*(p_ai.*1e5 - p_wi.*1e5));
w_ro    = par.PI.*1e-6.*(par.p_res.*1e5 - p_bh.*1e5);% w_ro = (-IPR.a + sqrt(IPR.a.^2+4*IPR.b.*(p_res-p_bh).*1e5))./(2.*IPR.b);
w_rg    = 1e1.*par.GOR.*w_ro;
w_to    = (m_or.*1e3./(m_gr.*1e3+m_or.*1e3)).*w_pr; %29
w_tg    = (m_gr.*1e3./(m_gr.*1e3+m_or.*1e3)).*w_pr; %30

% differential equations in discrete time
df1 = m_ga + par.tSim.*(w_gl - w_iv).*1e-3;
df2 = m_gt + par.tSim.*(w_iv + w_rg.*1e-1 - w_pg).*1e-3;
df3 = m_ot + par.tSim.*(w_ro - w_po).*1e-3;
df4 = m_gr + par.tSim.*(w_pg(1)+w_pg(2) - w_tg).*1e-3;
df5 = m_or + par.tSim.*(w_po(1)+w_po(2) - w_to).*1e-3;


% Concatenate the differential and algebraic equations
diff = vertcat(df1,df2,df3,df4,df5);

% ------- not used -------- %
% differential equations in continous time
df1C = (w_gl - w_iv).*1e-3;
df2C = (w_iv + w_rg.*1e-1 - w_pg).*1e-3;
df3C = (w_ro - w_po).*1e-3;
df4C = (sum(w_pg) - w_tg).*1e-3;
df5C = (sum(w_po) - w_to).*1e-3;

% Concatenate the differential and algebraic equations
diffC = vertcat(df1C,df2C,df3C,df4C,df5C); % continous time model
% ------- not used -------- %

% concatenate the differential and algebraic states
x_EKF = vertcat(m_ga,m_gt,m_ot,m_gr,m_or);
p_EKF = vertcat(w_gl);

% stage cost
L = -(w_to.^2) + 0.5.*sum(w_gl).^2;

nxEKF = length(diff);
%%

f_EKF = Function('f_EKF',{x_EKF,p_EKF},{diff},{'x','p'},{'xdot'});
JacFx = Function('JacFx',{x_EKF,p_EKF},{jacobian(diff,x_EKF)});

y_model = vertcat(p_ai,p_wh,p_bh,p_rh,p_m,w_to,w_tg);

h_EKF = Function('h_EKF',{x_EKF,p_EKF},{y_model});
JacHx = Function('JacHx',{x_EKF,p_EKF},{jacobian(y_model,x_EKF)});
yIndex = [1:4,7,8,25,27,29,30];

z_vec = vertcat(p_ai,p_wh,p_wi,p_bh,rho_ai,rho_m,w_iv,w_pc,w_pg,w_po,w_ro,w_rg,p_rh,rho_r,p_m,w_pr,w_to,w_tg);
z_EKF = Function('z_EKF',{x_EKF,p_EKF},{z_vec});

% ------- not used -------- %
JacAx = Function('JacAx',{x_EKF,p_EKF},{jacobian(diffC,x_EKF)});
JacBu = Function('JacBu',{x_EKF,p_EKF},{jacobian(diffC,p_EKF)});
JacJx = Function('JacJx',{x_EKF,p_EKF},{jacobian(L,x_EKF)});
JacJu = Function('JacJu',{x_EKF,p_EKF},{jacobian(L,p_EKF)});
J_EKF = Function('J_EKF',{x_EKF,p_EKF},{L});
% ------- not used -------- %

