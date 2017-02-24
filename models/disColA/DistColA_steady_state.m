%DISTCOLA_STEADY_STATE Summary of this function goes here
% 
% [OUTPUTARGS] = DISTCOLA_STEADY_STATE(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2016/02/05 23:59:07 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2016

% this code is copied from http://www.nt.ntnu.no/users/skoge/publications/thesis/2011_jacobsen/matlab-gproms-unisim-files/Chapter3-onecolumn/Optimize_col.m
% TO DO:
% - add cstr in the model (read Vlad's model in python!)

%close all;
%clc;
%clear all;
format long;
%Try to optimize distillation column "A" for different disturbances!

% Define nominal inputs!

% Nominal inputs
% LT=2.70629;                          % Reflux
% VB=3.20629;                          % Boilup
% D=0.5;                               % Distillate
% B=0.5;                               % Bottoms
% F=1.0;                               % Feedrate
%LT=0.1;                          % Reflux
%VB=0.1;                          % Boilup

LT  = 2.827;                          % Reflux
VB  = 3.454;                          % Boilup
%D   = 0.627;                          % Distillate
%B   = 0.573;                          % Bottoms
% D   = 0.627000;
% B   = 0.373000;
D   = 0.5;
B   = 0.5;
F   = 1.2;                            % Feedrate
F_0 = 1.2;                            % CSTR feedrate 
zF  = 0.65;                           % Feed composition at CSTR 
qF  = 1.0;                            % Feed liquid fraction
%Uinit = [LT VB D B F zF qF]';
%Uinit = [LT VB D B F zF qF F_0]';
%Uinit = [LT VB D B F zF qF]';
Uinit = [LT VB D B F zF]';
%global NT d
%NT = 41 ;

%Now we split this set in three: LT and VB are our unused inputs which can
%be used for optimization. D and B are used to control the holdups in the
%column. F, zF and qF are considered disturbances (but we leave qF
%unchanged for now).

u0 = Uinit(1:2);
%u0 = Uinit(1:5);   % change control inputs to LT, VB, D, B, F
%u0 = Uinit;
%u0 = Uinit(1:4);

%options = optimset('Display','final','TolFun',1e-6,'TolCon',1e-5,'Algorithm','active-set','MaxIter',75,'MaxFunEvals',1000) ;
options = optimset('Display','iter','TolFun',1e-6,'TolCon',1e-5,'Algorithm','active-set') ;

lb = [0.1 0.1]';
ub = [10 4.008]';
% ub = [inf inf]';
%ub = [4.0 4.0]';
%lb = [0.1 0.1   0.1 0.1 0.1 0.0]';
%ub = [10  4.008 1.0 1.0 2.5 1.0]';

% pV = 0.1:0.02:0.3;
% F = [1.4];
pV = 0.01;
%F  = 1.2;

[u_opt,FVAL,EXITFLAG,OUTPUT,LAMBDA] = fmincon(@col_obj,u0,[],[],[],[],lb,ub,@col_constr,options);

% evaluate optimized control values
global u;
% F  = 1.0;
F  = 1.2;
% D  = 0.627000;
% B  = 0.373000;
% D  = 0.5;
% B  = 0.5;
zF = 0.65;
%u  = [u_opt; D; B; F; zF];
u  = [u_opt; F; zF];
load cstr_init.mat;
x0 = Xinit;
options = optimset('Display','none','TolFun',1e-10);
x = fsolve(@cola_lv_cstr,x0,options);

u_opt_ss = u_opt;
x_opt_ss = x;
save u_opt_cstr.mat u_opt_ss x_opt_ss;

% act = zeros(length(pV),length(F)) ;
% kmax = length(pV)*length(F) ;
% k = 0 ;
% 
% if length(pV) == 1 && length(F) > 1
%     uopt = zeros(2,length(F)) ;
%     copt = zeros(2,length(F)) ;
%     lopt = zeros(3,length(F)) ;
%     jopt = zeros(1,length(F)) ;
%     flags = zeros(1,length(F)) ;
% elseif length(F) == 1 && length(pV) > 1
%     uopt = zeros(2,length(pV)) ;
%     copt = zeros(2,length(pV)) ;
%     lopt = zeros(3,length(pV)) ;
%     jopt = zeros(1,length(pV)) ;
%     flags = zeros(1,length(pV)) ;
% end
% 
% if length(pV) == 1 && length(F) == 1
%     options.Display = 'iter' ;
% end
% 
% 
% 
% for i=1:length(pV)
%     for j=1:length(F)
%     %d = [F(j) zF qF pV(i)]' ;
%     d = [F(j) zF qF pV(i) F_0]';
%     k = k+1;
%     %[u_opt,FVAL,EXITFLAG,OUTPUT,LAMBDA] = fmincon(@col_obj,u0,[],[],[],[],lb,ub,[],options);
%     [u_opt,FVAL,EXITFLAG,OUTPUT,LAMBDA] = fmincon(@col_obj,u0,[],[],[],[],lb,ub,@col_constr,options);
%         if EXITFLAG <= 0 
%         disp(['Optimization number ' num2str(k) ' did not converge in ' num2str(options.MaxIter) ' iterations!']) 
%         act(i,j) = NaN ;
%         elseif EXITFLAG > 0
%         disp(['Optimization number ' num2str(k) ' of ' num2str(kmax) ' completed']) 
%         %Now check which constraint is active
%             if LAMBDA.ineqnonlin(1) == 0 && LAMBDA.ineqnonlin(2) == 0 && LAMBDA.upper(2) == 0 % No active constraints 
%                 act(i,j) = 1 ;
%             elseif LAMBDA.ineqnonlin(1) == 0 && LAMBDA.ineqnonlin(2) == 0 && LAMBDA.upper(2) > 0 % Vmax active
%                 act(i,j) = 2 ;
%             elseif LAMBDA.ineqnonlin(1) == 0 && LAMBDA.ineqnonlin(2) > 0 && LAMBDA.upper(2) > 0 % Vmax and XD active
%                 act(i,j) = 3 ;
%             elseif LAMBDA.ineqnonlin(1) == 0 && LAMBDA.ineqnonlin(2) > 0 && LAMBDA.upper(2) == 0 % XD active
%                 act(i,j) = 4 ;
%             elseif LAMBDA.ineqnonlin(1) > 0 && LAMBDA.ineqnonlin(2) > 0 && LAMBDA.upper(2) == 0 % XD and XB active
%                 act(i,j) = 5 ;
%             elseif LAMBDA.ineqnonlin(1) > 0 && LAMBDA.ineqnonlin(2) > 0 && LAMBDA.upper(2) > 0 % All constraints active! 
%                 act(i,j) = 6 ;
%             elseif LAMBDA.ineqnonlin(1) > 0 && LAMBDA.ineqnonlin(2) == 0 && LAMBDA.upper(2) > 0 % Vmax and XB active 
%                 act(i,j) = 7 ;
%             elseif LAMBDA.ineqnonlin(1) > 0 && LAMBDA.ineqnonlin(2) == 0 && LAMBDA.upper(2) == 0 % XB active
%                 act(i,j) = 8 ; 
%             end    
%         end  
%         
%         if length(pV) == 1 && length(F) > 1
%                 uopt(:,j) = u_opt ;
%                 C = col_constr(u_opt) ;
%                 copt(:,j) = C ;
%                 lopt(:,j) = [LAMBDA.upper(2) ; LAMBDA.ineqnonlin] ;
%                 jopt(j) = FVAL ;
%                 flags(j)=EXITFLAG ;
%         elseif length(F) == 1 && length(pV) > 1
%                 uopt(:,i) = u_opt ;
%                 C = col_constr(u_opt) ;
%                 copt(:,i) = C ;
%                 lopt(:,i) = [LAMBDA.upper(2) ; LAMBDA.ineqnonlin] ;
%                 jopt(i) = FVAL ;
%                 flags(i)=EXITFLAG ;
%          end
%      end
% end
% 
% if length(pV) >1 && length(F) > 1
% act = flipud(act) ;
% disp(act)
% elseif length(pV) == 1  && length(F) > 1
% %figure(1),plot(F,copt(1,:),'r',F,copt(2,:),'b') ;
% %figure(2),plot(F,uopt(1,:),'r',F,uopt(2,:),'b') ;
% %figure(3),plot(F,lopt(1,:),'r',F,lopt(2,:),'b',F,lopt(3,:),'gr') ;
% plot(F,jopt)
% 
% 
% matris = [act ; lopt]';
% disp(matris)
% elseif length(F) == 1 && length(pV) > 1
% figure(1),plot(pV,copt(1,:),'r',pV,copt(2,:),'b') ;
% figure(2),plot(pV,uopt(1,:),'r',pV,uopt(2,:),'b') ;  
% figure(3),plot(pV,lopt(1,:),'r',pV,lopt(2,:),'b',pV,lopt(3,:),'gr') ;
%  figure(4),plot(pV,jopt)
% act = act' ;
% matris = [act ; lopt] ;
% disp(matris)
% 
% elseif length(F) == 1 && length(pV) == 1
% [Jcol, D, B, xD, xB, xf] = col_obj(u_opt) ;
% end
% 
% % load cola_init.mat;
% % x_init = Xinit;
% % save u_opt_ss.mat u_opt x_init;
% save u_opt_ss.mat u_opt xf;
