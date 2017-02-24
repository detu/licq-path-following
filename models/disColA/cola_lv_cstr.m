function xprime=cola_lv_cstr(t,X) 
% sample usage:   [t,x]=ode15s('cola_lv',[0 5000],0.5*ones(1,82));
%
% cola_lv - Subroutine for simulation with LV-configuration.
%           It calls the model colamod, and 
%           includes control of condenser and reboiler level 
%           using two P-controllers with the LV-configuration. 
%
%            Inputs are reflux (LT) and boilup (VB). Disturbances
%            are feedrate and feed composition. These are set by directly
%            altering 'cola_lv.m'. Outputs are liquid composition and
%            liquid hold up for stages 1 through NT, given in x. 

% Number of stages in the column


% global u NT %Make perturbed inputs/disturbances available to model
global uc;
NT = 41;
 
% LT  = 2.827;                        % Reflux
% VB  = 3.454;                        % Boilup
LT  = uc(1);
VB  = uc(2);
% D   = 0.5;                          % Distillate
% B   = 0.5;                          % Bottoms
%F   = 1.5;
F   = uc(3);
F_0 = 0.3;                          % CSTR feedrate 1.1 after 1.0 
%F_0 = 0.4;
%F_0 = uc(3);
%zF  = 0.45;
zF  = 1.0;
% qF  = 1.0;                          % Feed liquid fraction
% Uinit = [LT VB D B F zF qF F_0]';
% u1 = Uinit(1:2);
% pV = 0.01;
% d  = [F zF qF pV]';
% uc = [u1 ; d(1:3)];


% % P-Controllers for control of reboiler and condenser hold up.
% KcB = 10;    % controller gains
% KcD = 10;         
% MDs = 0.5;   % Nominal holdups - these are rather small
% MBs = 0.5;          
% Ds  = 0.5;   % Nominal flows
% Bs  = 0.5;          
% %Ds  = uc(4);   % Nominal flows
% %Bs  = uc(5);    
% MB  = X(NT+2);  
% MD  = X(2*NT+1);
% D   = Ds+(MD-MDs)*KcD;       % Distillate flow
% B   = Bs+(MB-MBs)*KcB;       % Bottoms flow     
D = uc(4);
B = uc(5);



% Store all inputs and disturbances
u_all(1:2) = [LT;VB];
u_all(3)   = D; 
u_all(4)   = B;
u_all(5)   = F;
u_all(6)   = zF;
u_all(7)   = F_0;


xprime=colamod_cstr(t,X,u_all);