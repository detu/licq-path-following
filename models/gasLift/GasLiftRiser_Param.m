function  par = GasLiftRiser_Param

par.n_w     = 2;


par.L_w     = [1500;1500];
par.H_w     = [1000;1000];
par.D_w     = [0.121;0.121];

par.L_bh    = [500;500];
par.H_bh    = [500;500];
par.D_bh    = [0.121;0.121];

par.L_a     = par.L_w;
par.H_a     = par.H_w;
par.D_a     = [0.189;0.189];

par.L_r     = 500;
par.H_r     = 500;
par.D_r     = 0.121;

par.rho_o   = [8;8].*1e2;
par.C_iv    = [0.1e-3;0.1e-3];
par.C_pc    = [2e-3;2e-3];
par.C_pr    = [10e-3];

par.GOR     = [0.1;0.12];
par.p_m     = [20;20];
par.p_res   = [150;155];
par.PI      = [7;7];
par.T_a     = [28+273;28+273];
par.T_w     = [32+273;32+273];
par.T_r     = 30+273;
par.p_s     = 20;

par.GOR_var = [0.05;0.02];
par.rho_var = [150;25].*0;

par.wS      = [1;1]; % weights for the different scenarios

par.R       = 8.314;
par.Mw      = 20e-3;

par.nIter   = 60;
%par.qGLMax  = 40;
par.qGLMax  = 5;
%par.QgMax   = 12;
par.QgMax   = 9.5;
%par.QgMax   = 9;  % EKA, For initial condition setting.
%par.QgMax   = 8; % 
par.qGLROC  = [0.5;0.5];

