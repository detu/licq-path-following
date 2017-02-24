%COLLOCATION_EX Summary of this function goes here
% 
% [OUTPUTARGS] = COLLOCATION_EX(INPUTARGS) Explain usage here
% 
% Collocation example from Zavala's thesis
%
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2016/02/12 13:12:22 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2016

% TODO:
% - COMPARE RESULT WITH ZAVALA'S ! SEE GRAPH.M FILE !
% - ADD EQUALITY CONSTRAINTS !
% - ADD BOUND CONSTRAINTS !

import casadi.*

a = [ 0.19681547722366,   0.39442431473909, 0.37640306270047;
      -0.06553542585020,  0.29207341166523, 0.51248582618842;
      0.02377097434822,  -0.04154875212600, 0.11111111111111];
      
% mathematical model parameters
jj     = 100;
cf     = 7.6;          
alpha  = 1.95e-04;
tf     = 300;
k10    = 300;
tc     = 290;
n      = 5;
alpha1 = 1e6;
alpha2 = 2e3;
alpha3 = 1e-03;
nfe    = 40;
ncp    = 3;
r1     = 0.15505102572168;
r2     = 0.64494897427832;
r3     = 1;

% initial and end transition points
c_des  = 0.0944;
t_des  = 0.7766;
u_des  = 340;
c_init = 0.1367;
t_init = 0.7293;
u_init = 390;
time   = 9;
theta  = 20;

% initial guesses of the decision variables
point  = 0;
slopec = (c_des-c_init)/(nfe*ncp);
slopet = (t_des-t_init)/(nfe*ncp);
slopeu = (u_des-u_init)/(nfe*ncp);

%fe = zeros(nfe);  %number of finite elements
%cp = zeros(ncp);  %number of collocation points

% c = zeros(nfe,ncp); 
% t = zeros(nfe,ncp);  
% u = zeros(nfe,ncp); 
% cdot = zeros(nfe,ncp);
% tdot = zeros(nfe,ncp);
c{(nfe),(ncp)}    = [];
t{(nfe),(ncp)}    = [];
u{(nfe),(ncp)}    = [];
cdot{(nfe),(ncp)} = [];
tdot{(nfe),(ncp)} = [];

%param a{cp,cp} ;    # collocation matrix
h = zeros(nfe,1);   %finite element length

% define symbolic variables for decision variables (c(i,j), t(i,j), u(i,j),
% cdot(i,j), and tdot(i,j)
%nv = 5*nfe*ncp;
nv = 3*nfe*ncp;         % c,t,and u
V  = SX.sym('V',nv,1);
vars_init = zeros(nv,1);
vars_lb   = zeros(nv,1);
%vars_ub   = inf*ones(nv,1);
offset = 1;

% assign c(i,j) to symbolic variable
for i=1:nfe
    for j=1:ncp
        c{i,j}              = V(offset,1);
        point               = point + 1;                        % remember to reset (point = 0) after c(i,j) !
        vars_init(offset,1) = slopec*point + c_init;
        offset              = offset + 1;
    end
    h(i) = 1/nfe ;
end

% assign t(i,j) to symbolic variable
point = 0;
for i=1:nfe
    for j=1:ncp
        t{i,j}              = V(offset,1);
        point               = point + 1;
        vars_init(offset,1) = slopet*point + t_init;
        offset              = offset + 1;
    end
end

% assign u(i,j) to symbolic variable
point = 0;
for i=1:nfe
    for j=1:ncp
        u{i,j}              = V(offset,1);
        point               = point + 1;
        vars_init(offset,1) = slopeu*point + u_init;
        offset              = offset + 1;
    end
end

% % assign cdot(i,j) to symbolic variable
% for i=1:nfe
%     for j=1:ncp
%         cdot{i,j}           = V(offset,1);
%         vars_init(offset,1) = 0;
%         offset              = offset + 1;
%     end
% end
% 
% % assign tdot(i,j) to symbolic variable
% for i=1:nfe
%     for j=1:ncp
%         tdot{i,j}           = V(offset,1);
%         vars_init(offset,1) = 0;
%         offset              = offset + 1;
%     end
% end

% for i=1:nfe
%     for j=1:ncp
%         point  = point + 1;
%         c(i,j) = slopec*point+c_init;
%         t(i,j) = slopet*point+t_init;
%         u(i,j) = slopeu*point+u_init;
%     end
%     h(i) = 1/nfe ;
% end

% auxiliary equations
yc = tc/(jj*cf);
yf = tf/(jj*cf);

% states first order derivatives
%var cdot{i in fe, j in cp} = (1-c[i,j])/theta-k10*exp(-n/t[i,j])*c[i,j]                           ;
%var tdot{i in fe, j in cp} = (yf-t[i,j])/theta+k10*exp(-n/t[i,j])*c[i,j]-alpha*u[i,j]*(t[i,j]-yc) ;

for i=1:nfe
    for j=1:ncp
        cdot{i,j} = (1-c{i,j})/theta-k10*exp(-n/t{i,j})*c{i,j};
        %tdot{i,j} = (yf-t{i,j})/theta+k10*exp(-n/t{i,j})*c{i,j}-alpha*u{i,j}*(t{i,j}-yc);
    end
end

for i=1:nfe
    for j=1:ncp
        tdot{i,j} = (yf-t{i,j})/theta+k10*exp(-n/t{i,j})*c{i,j}-alpha*u{i,j}*(t{i,j}-yc);
    end
end

% for i=1:nfe
%     for j=1:ncp
%         cdot(i,j) = (1-c{i,j})/theta-k10*exp(-n/t{i,j})*c{i,j};
%         tdot(i,j) = (yf-t{i,j})/theta+k10*exp(-n/t{i,j})*c{i,j}-alpha*u{i,j}*(t{i,j}-yc);
%     end
% end

% % collocation equations
% fecolc{i in fe diff{1},j in cp}: c[i,j] = c[i-1,ncp]+time*h[i]*sum{k in cp} a[k,j]*cdot[i,k];
% fecolt{i in fe diff{1},j in cp}: t[i,j] = t[i-1,ncp]+time*h[i]*sum{k in cp} a[k,j]*tdot[i,k];
% 
% fecolc0{i in 1..1,j in cp}:       c[i,j] = c_init+time*h[i]*sum{k in cp} a[k,j]*cdot[i,k];
% fecolt0{i in 1..1,j in cp}:       t[i,j] = t_init+time*h[i]*sum{k in cp} a[k,j]*tdot[i,k];

% CHANGE THIS TO EQUALITY CONSTRAINT !
g = {};
for i=1:nfe
    for j=1:ncp
        if i == 1
            
            c_sum = 0;
            t_sum = 0;
            for k=1:ncp
                c_sum = c_sum + a(k,j)*cdot{i,k};
                t_sum = t_sum + a(k,j)*tdot{i,k};
            end
            %c{i,j} = c_init + time*h(i)*c_sum;
            %t{i,j} = t_init + time*h(i)*t_sum;
            %g = {g{:}, c{i,j} - c_init - time*h(i)*c_sum};
            %g = {g{:}, t{i,j} - t_init - time*h(i)*t_sum};
            g = {g{:}, c_init + time*h(i)*c_sum - c{i,j}};
            g = {g{:}, t_init + time*h(i)*t_sum - t{i,j}};
        else
            
            c_sum = 0;
            t_sum = 0;
            for k=1:ncp
                c_sum = c_sum + a(k,j)*cdot{i,k};
                t_sum = t_sum + a(k,j)*tdot{i,k};
            end
            %c{i,j} = c{i-1,ncp} + time*h(i)*c_sum;
            %t{i,j} = t{i-1,ncp} + time*h(i)*t_sum; 
            %g = {g{:}, c{i,j} - c{i-1,ncp} - time*h(i)*c_sum};
            %g = {g{:}, t{i,j} - t{i-1,ncp} - time*h(i)*t_sum};
            g = {g{:}, c{i-1,ncp} + time*h(i)*c_sum - c{i,j}};
            g = {g{:}, t{i-1,ncp} + time*h(i)*t_sum - t{i,j}};
        end
    end
end


% % objective function...
% 
% minimize cost: sum{i in fe} (h[i]*sum{j in cp} ((alpha1*(c[i,j]-c_des)^2+ alpha2*(t[i,j]-t_des)^2+alpha3*(u[i,j]-u_des)^2 )*a[j,ncp]));

% construct objective function 
obj = 0;
for i = 1:nfe
    o_sum = 0;
    for j = 1:ncp
        o_sum = o_sum + ((alpha1*(c{i,j}-c_des)^2+ alpha2*(t{i,j}-t_des)^2+alpha3*(u{i,j}-u_des)^2 )*a(j,ncp));
    end
    obj = obj + h(i)*o_sum;
end

%% SOLVE NLP
nlp = struct('x',V,'f',obj,'g',vertcat(g{:}));

% Set options
opts = struct;
%opts.expand = True;
opts.ipopt.max_iter = 50000;
%opts["linear_solver"] = 'ma27'
opts.ipopt.hessian_approximation = 'exact';
%opts.ipopt.hessian_approximation = 'limited-memory';

solver = nlpsol('solver', 'ipopt', nlp, opts);
arg    = struct;

% Initial condition
arg.x0 = vars_init;

% Bounds on x
arg.lbx = vars_lb;
%arg.ubx = vars_ub;

arg.ubg = 0;
arg.lbg = 0;

% Solve the problem
tic;
res = solver(arg);
toc;


% Print the optimal cost
fprintf('optimal cost: %f \n', full(res.f));

% Retrieve the solution
v_opt = full(res.x);
lamda_opt = full(res.lam_g);