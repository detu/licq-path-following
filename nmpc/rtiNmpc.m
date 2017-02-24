function [Tall, xmeasureAll, uAll, ObjVal, primalPF, params, runtime] = rtiNmpc(optProblem, system, mpciterations, N, T, tmeasure, xmeasure, u0, varargin)
%PFNMPC Summary of this function goes here
% 
% Path-following based Nonlinear Model Predictive Control
%
% [OUTPUTARGS] = PFNMPC(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2016/04/06 20:39:51 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2016

% dimension of state x and control u
nx     = numel(xmeasure);
nu     = size(u0,1);
Tall   = [];
Xall   = zeros(mpciterations,size(xmeasure,1));
Uall   = zeros(mpciterations,size(u0,1));
ObjVal = zeros(mpciterations,1);

xmeasureAll = [];
uAll        = [];
runtime     = [];

u_pf_opt    = [];
x_pf_opt    = [];

x_nlp_opt   = [];
u_nlp_opt   = [];

load noise1pct.mat;
%load noise3pct.mat;
%load noise5pct.mat;
z1 = xmeasure;
mpcStartup = 1;    % noise case
%mpcStartup = 3;   % noise-free
%mpcStartup = 5;
%mpcStartup = 10;

% Start of the NMPC iteration
mpciter = 1;
while(mpciter <= mpciterations)
    
    fprintf('-----------------------------\n');
    fprintf('MPC iteration: %d\n', mpciter);

    %   Obtain new initial value
    [t0, x0] = measureInitialValue ( tmeasure, xmeasure );
    
    % check constraints on boundary
    epsilon  = 0;
   
    %holdupNoise   = x0(43:end) .* noise;
    %holdupNoise   = noise;
    holdupNoise   = noise(:,mpciter);
    concentrationNoise = zeros(42,1);
    measNoise  = [concentrationNoise;holdupNoise];
    x0_measure =  x0 + measNoise;
    
    %x0_measure =  x0;  % noise-free
    %x0_measure = safetyBound(x0_measure,epsilon);

    if mpciter == 1
        % startup NMPC for computing iniital primal and dual
        [primalNLP, dualNLP, lb, ub, objValue, params, elapsedtime] = solveOptimalControlProblem(optProblem, system, N, t0, x0, u0, T, mpciter, u_pf_opt, x_pf_opt, z1);
        [u_nlp_opt, x_nlp_opt] = plotStatesN(primalNLP, lb, ub, N);
    else
        
        p_final = x0_measure; % with  noise
        if mpciter == 2
            p_init  = primalNLP(1:nx);
            xstart  = primalNLP;
            ystart  = dualNLP;
        else
            p_init  = primalPF(1:nx);
            xstart  = primalPF;
            ystart  = dualPF;
        end
        delta_t = 1.0;
        lb_init = lb;
        ub_init = ub;
        [primalPF, dualPF, elapsedqp] = jpredictor_licq_pure_3(@(p)distColACstr_pn(p), p_init, p_final, xstart, ystart, delta_t, lb_init, ub_init, 0, N);
        
        %[u_nlp_opt, x_nlp_opt] = plotStatesN(primalNLP, lb, ub, N);  % distillation column plot
        [u_pf_opt, x_pf_opt] = plotStatesN(primalPF, lb, ub, N);
        z1 = x_pf_opt(1:nx,5); % 4 = (d+1) + 1 (d=number of collocation point)
        
    end
    
    
    Tall   = [Tall t0];
    Xall(mpciter+1,:)   = x0';
    Uall(mpciter+1,:)   = u0(:,1);
    ObjVal(mpciter,1)   = objValue;
    
    % Apply control to process with optimized control from path-following
    % algorithm. it may have noise to measurement...
    x0     = xmeasure;  % from the online step 
    if mpciter == 1
        [tmeasure, xmeasure] = applyControl(system, T, t0, x0, u_nlp_opt);
        uAll                 = [uAll;u_nlp_opt(:,1)];   
        u0                   = shiftHorizon(u_nlp_opt);
        runtime              = [runtime;elapsedtime];
    else
        [tmeasure, xmeasure] = applyControl(system, T, t0, x0, u_pf_opt);
        uAll                 = [uAll;u_pf_opt(:,1)];   % Distilation column A
        u0                   = shiftHorizon(u_pf_opt);
        runtime              = [runtime;elapsedqp];
    end
    
    % PLOT HERE !!!
    xmeasureAll = [xmeasureAll;xmeasure];
    mpciter     = mpciter+1;
end

%xmeasureAll = reshape(xmeasureAll,3,mpciterations);  %bailey-cstr
%xmeasureAll = reshape(xmeasureAll,82,mpciterations);  %distillation column A
xmeasureAll = reshape(xmeasureAll,84,mpciterations);  %CSTR + distillation column A

end

function [t0, x0] = measureInitialValue ( tmeasure, xmeasure )
    t0 = tmeasure;
    x0 = xmeasure;
end

function [tapplied, xapplied] = applyControl(system, T, t0, x0, u0)
%     global sf1 sf2;
%     u0s = u0;
%     u0s(1,:) = u0s(1,:)/sf1;
%     u0s(2,:) = u0s(2,:)/sf2;
    xapplied = dynamic(system, T, t0, x0, u0(:,1));
    %xapplied = dynamic(system, T, t0, x0, u0s(:,1));
    tapplied = t0+T;
end

function u0 = shiftHorizon(u)
    u0 = [u(:,2:size(u,2)) u(:,size(u,2))]; 
end

function [u, lamda, lbw, ubw, objVal, params, elapsednlp] = solveOptimalControlProblem(optProblem, system, N, t0, x0, u0, T, mpciter, u_nlp, x_nlp, x0_measure)

    import casadi.* 
    global ns;
%     global sf1 sf2;
    
    % call ODE15s N-times initial guess in optimization
    x(1,:) = x0';
    for k=1:N
        x(k+1,:) = x0';    % ONE-TIME SIMULATION !
    end
    
%     x(:,1) = x0;
%     load CstrDistXinit.mat;
%     xf     = Xinit(1:84);
%     xtem   = computeInitGuess(x0,xf,N);
%     x      = [x xtem];
%     x      = x';

%     %x(:,1) = x0;
%     if mpciter <=1
%         x(:,1) = x0;
%         load CstrDistXinit.mat;
%         xf     = Xinit(1:84);
%         xtem   = computeInitGuess(x0,xf,N);
%         x      = [x xtem];
%         x      = x';
%         %[J,g,w,lbg,ubg,lbw,ubw] = optProblem(t, x, u);
%         %[J,g,w0,w,lbg,ubg,lbw,ubw] = optProblem(x0, u0);
%         [J,g,w0,w,lbg,ubg,lbw,ubw,params] = optProblem(x, u0, N, x0_measure);
%     else
%         x(:,1) = [x0;zeros(ns,1)];
%         x_nlp(:,1) = [];    % delete the first column
%         x      = [x x_nlp];
%         x      = x';
%         [J,g,w0,w,lbg,ubg,lbw,ubw,params] = optProblem(x, u_nlp, N, x0_measure);
%     end
    
    %[J,g,w,lbg,ubg,lbw,ubw] = optProblem(t, x, u);
    %[J,g,w0,w,lbg,ubg,lbw,ubw] = optProblem(x0, u0);
    [J,g,w0,w,lbg,ubg,lbw,ubw,params] = optProblem(x, u0, N, x0_measure);
   
    prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));
    %solver = nlpsol('solver', 'ipopt', prob); 
    options = struct;
    options.ipopt.tol       = 1e-12;
    options.ipopt.constr_viol_tol    = 1e-10; 
    solver = nlpsol('solver', 'ipopt', prob, options);

    % Solve the NLP
    startnlp = tic;
    sol   = solver('x0', w0, 'lbx', lbw, 'ubx', ubw, 'lbg', lbg, 'ubg', ubg);
    %sol   = solver('x0', w0, 'lbg', lbg, 'ubg', ubg);  % without bound constraint... drive system immediately to steady-state equilibrium! (with combination objective function, IPOPT doesn't converge quickly!)
    elapsednlp = toc(startnlp);
    fprintf('IPOPT solver runtime = %f\n',elapsednlp);

    u      = full(sol.x);
    %lamda  = full(sol.lam_g);
    lamda.lam_g = full(sol.lam_g);
    lamda.lam_x = full(sol.lam_x);
    objVal = full(sol.f);
end

function [x, t_intermediate, x_intermediate] = dynamic(system, T, t0, x0, u0)
    x = system(t0, x0, u0, T);
    x_intermediate = [x0; x];
    t_intermediate = [t0, t0+T];
end

