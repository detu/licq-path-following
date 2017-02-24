%function [Jcol, D, B, xD, xB, x] = col_obj(u1)
function [Jcol] = col_obj(u1)
% global d u NT
% 
% %u = [u1 ; d(1:3)] ; %L, V, F, zF, qF
% u = [u1 ; d(1:3); d(5)]; %L, V, F, zF, qF, F_0
global u;
%F  = 1.0;
F  = 1.2;
zF = 0.65;
u  = [u1; F; zF];
NT = 41;
load cstr_init.mat;
x0 = Xinit;
options = optimset('Display','none','TolFun',1e-10);

x = fsolve(@cola_lv_cstr,x0,options);

KcB=10;  KcD=10;         % controller gains
% KcB=1e-2;  KcD=1e-2;         % controller gains
%KcB=1;  KcD=1;
MDs=0.5; MBs=0.5;        % Nominal holdups - these are rather small  
Ds=0.5; Bs=0.5;          % Nominal flows
%MB=x(NT+1);  MD=x(2*NT); % Actual reboiler and condenser holdup
MB=x(NT+2);  MD=x(2*NT+1);
D=Ds+(MD-MDs)*KcD;       % Distillate flow
B=Bs+(MB-MBs)*KcB;       % Bottoms flow     

% D = 0.627000;
% B = 0.373000;
% D = 0.5;
% B = 0.5;
 
% Jcol = pf*F + pV*V - pB*B - pD*D ;
% prices
pf = 1; 
pV = 0.01; 
pB = 1; 
pD = 2;

%Jcol = pf*u1(5) + pV*u1(2) - pB*u1(4) - pD*u1(3);
Jcol = pf*F + pV*u1(2) - pB*B - pD*D;

end
