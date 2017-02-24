function [C, Ceq]= col_constr(u1)

%global d u NT
global u;
%F  = 1.0;
F  = 1.2;
% D  = 0.627000;
% B  = 0.373000;
% D  = 0.5;
% B  = 0.5;
zF = 0.65;
% u  = [u1; D; B; F; zF];
u  = [u1; F; zF];
%u = [u1 ; F zF D B] ; %L, V, F, zF, qF
%u = [u1 ; d(1:3); d(5)]; %L, V, F, zF, qF, F_0

%u = u1;  % u = [LT VB D B F zF]

%x0 = 0.5*ones(84,1);
load cstr_init.mat;
x0 = Xinit;
options = optimset('Display','none','TolFun',1e-10) ;
%x = fsolve(@cola_lv,x0,options) ;
x = fsolve(@cola_lv_cstr,x0,options) ;
NT = 41;
xB = x(1); 
xD = x(NT);

% C = [xB-0.08;
%     0.95-xD];
C = [xB-0.008;
    0.95-xD];
% C = [xB-0.10;   % relax the constraint.. 
%      0.95-xD];
%C = 0.95-xD;
Ceq=[];

