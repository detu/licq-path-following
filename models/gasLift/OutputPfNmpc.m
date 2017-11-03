function [Tall, xmeasureAll, uAll, ObjVal, primalPF, params, runtime] = OutputPfNmpc(optProblem, system, mpciterations, N, T, tmeasure, xmeasure, u0, paramModel, varargin)
%OUTPUTPFNMPC Summary of this function goes here
% 
% [OUTPUTARGS] = OUTPUTPFNMPC(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2017/10/26 15:07:28 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2017

% dimension of state x and control u
global nx nu;
nx     = numel(xmeasure);
nu     = size(u0,1);
Tall   = [];
Xall   = zeros(mpciterations,size(xmeasure,1));
Uall   = zeros(mpciterations,size(u0,1));

xmeasureAll = [];
uAll        = [];
runtime     = [];

u_pf_opt    = [];
x_pf_opt    = [];

ObjVal      = [];

load noise1pct.mat;
z1 = xmeasure;
% Start of the NMPC iteration
mpciter = 1;
while(mpciter <= mpciterations)
    
    fprintf('-----------------------------\n');
    fprintf('MPC iteration: %d\n', mpciter);

    %   Obtain new initial value
    [t0, x0] = measureInitialValue ( tmeasure, xmeasure );
    
    
    % add measurement error to x0
    %measNoise          = noise(1:nx,mpciter);
    measNoise          = [noise(1:8,mpciter);zeros(30,1)];
    x0_measure         =  x0 + 0*measNoise;
    
    % advanced-step NMPC:
    [primalNLP, dualNLP, lb, ub, ~, params] = solveOptimalControlProblem(optProblem, system, N, t0, x0, u0, T, mpciter, u_pf_opt, x_pf_opt, z1, paramModel);
    
    % re-arrange NLP solutions
    %[~, x_nlp_opt] = plotStatesN(primalNLP, lb, ub, N);
    %[u_nlp_opt, x_nlp_opt] = processNLPResult(primalNLP, paramModel); 
    [u_nlp_opt, x_nlp_opt] = plotStatesGL(primalNLP, lb, ub, N);
    
    p_init  = primalNLP(1:nx);        
    p_final = x0_measure; 
    xstart  = primalNLP;
    ystart  = dualNLP;
    
    % choose number of path-following step
    %delta_t = 0.25;
    %delta_t = 0.5;   
    delta_t = 1.0;
    
    lb_init = lb;
    ub_init = ub;
    
    % NLP sensitivity (predictor-corrector)
    [primalPF, ~, elapsedqp] = jpredictorLicq(@(p,paramModel)gasLiftRiser(p,paramModel), p_init, p_final, xstart, ystart, delta_t, lb_init, ub_init, 0, N, paramModel);
    %[primalPF, ~, elapsedqp] = pf_pc_mfcq(@(p)distColACstr_pn(p), p_init, p_final, xstart, ystart, delta_t, lb_init, ub_init, 0, N);
    
     
    %[u_nlp_opt, x_nlp_opt] = plotStatesN(primalNLP, lb, ub, N);  % distillation column plot
    [u_pf_opt, x_pf_opt] = plotStatesGL(primalPF, lb, ub, N);  
    z1 = x_pf_opt(1:nx,5); % 4 = (d+1) + 1 (d=number of collocation point)

    % store output variables
    Tall   = [Tall t0];
    Xall(mpciter+1,:)   = x0';
    Uall(mpciter+1,:)   = u0(:,1);
    
    % Apply control to process with optimized control from path-following
    % algorithm. 
    x0                   = xmeasure;  % from the online step
    %x0                   = x0_measure; 
    [tmeasure, xmeasure] = applyControl(system, T, t0, x0, u_pf_opt(:,1), paramModel);
    %ObjVal(mpciter) = computeObjectiveFunctionValues(u_pf_opt(:,1),xmeasure); % USING ACTUAL STATE!
    
    
    % collect variables
    xmeasureAll = [xmeasureAll;xmeasure];
    uAll        = [uAll;u_pf_opt(:,1)];
    runtime     = [runtime;elapsedqp];
    
    % prepare restart: u_new is from path-following output !
    u0      = shiftHorizon(u_pf_opt);
    mpciter = mpciter+1;
end

xmeasureAll = reshape(xmeasureAll,38,mpciterations);  

end

function [t0, x0] = measureInitialValue ( tmeasure, xmeasure )
    t0 = tmeasure;
    x0 = xmeasure;
end

function [tapplied, xapplied] = applyControl(system, T, t0, x0, u0, paramModel)
    xapplied = dynamic(system, T, t0, x0, u0(:,1), paramModel);
    tapplied = t0+T;
end

function u0 = shiftHorizon(u)
    u0 = [u(:,2:size(u,2)) u(:,size(u,2))];
end

function [u, lamda, lbw, ubw, objVal, params, elapsednlp] = solveOptimalControlProblem(optProblem, system, N, t0, x0, u0, T, mpciter, u_nlp, x_nlp, x0_measure, paramModel)

    import casadi.*
%     % call ODE15s N-times initial guess in optimization
%     x(1,:) = x0';
%     for k=1:N
%         x(k+1,:) = x0';    % ONE-TIME SIMULATION !
%     end
    x = x0;

    [J,g,w0,w,p,lbg,ubg,lbw,ubw,params] = optProblem(x, u0, N, x0_measure, paramModel);
    prob = struct('p', p, 'f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));
    options = struct;
    %options.ipopt.tol       = 1e-12;
    %options.acceptable_compl_inf_tol    = 1e-6; 
    solver = nlpsol('solver', 'ipopt', prob, options);

    
    % Solve the NLP
    startnlp   = tic;
    %sol        = solver('x0', w0, 'lbx', lbw, 'ubx', ubw, 'lbg', lbg, 'ubg', ubg);
    sol         = solver('x0', w0, 'p', paramModel.GOR, 'lbx', lbw, 'ubx', ubw, 'lbg', lbg, 'ubg', ubg);
    elapsednlp = toc(startnlp);
    fprintf('IPOPT solver runtime = %f\n',elapsednlp);
    success = strcmp(solver.stats.return_status,'Infeasible_Problem_Detected');
    if (success)
        keyboard;
    end

    %% EXTRACT OPTIMAL SOLUTIONS !
    u           = full(sol.x);
    lamda.lam_g = full(sol.lam_g);
    lamda.lam_x = full(sol.lam_x);
    objVal      = full(sol.f);
    
    
end

function [x, t_intermediate, x_intermediate] = dynamic(system, T, t0, x0, u0, paramModel)
    x = system(t0, x0, u0, T, paramModel);
    x_intermediate = [x0; x];
    t_intermediate = [t0, t0+T];
end

function [u_nlp, x_nlp] = processNLPResult(w_opt, par)
    nu    = par.n_w;
    nz    = 30;
    nx    = 8;
    d     = 3;
    sv    = par.sv;
    
    n_w_i         = nx+nz+nu+ (nu + (nx+nz+sv)*d+nx)*(par.N) ;
    u_opt1        = [w_opt((nx+nz+nu+1):d*(nx+nz+sv)+nx+nu:n_w_i);NaN];
    u_opt2        = [w_opt((nx+nz+nu+2):d*(nx+nz+sv)+nx+nu:n_w_i);NaN];
    x_nlp.w_po1   = w_opt([(nx+nz-11),nu+(nx+nz-11)+(nx+nz+nu):d*(nx+nz+sv)+nx+nu:n_w_i]);
    x_nlp.w_po2   = w_opt([(nx+nz-10),nu+(nx+nz-10)+(nx+nz+nu):d*(nx+nz+sv)+nx+nu:n_w_i]);
    x_nlp.w_oEx   = w_opt([(nx+nz-1),nu+(nx+nz-1)+(nx+nz+nu):d*(nx+nz+sv)+nx+nu:n_w_i]);
    x_nlp.w_pg1   = w_opt([(nx+nz-13),nu+(nx+nz-13)+(nx+nz+nu):d*(nx+nz+sv)+nx+nu:n_w_i]);
    x_nlp.w_pg2   = w_opt([(nx+nz-12),nu+(nx+nz-12)+(nx+nz+nu):d*(nx+nz+sv)+nx+nu:n_w_i]);
    x_nlp.w_gEx   = w_opt([(nx+nz-0),nu+(nx+nz-0)+(nx+nz+nu):d*(nx+nz+sv)+nx+nu:n_w_i]);
    x_nlp.p_bh1   = w_opt([(nx+nz-23),nu+(nx+nz-23)+(nx+nz+nu):d*(nx+nz+sv)+nx+nu:n_w_i]);
    x_nlp.p_bh2   = w_opt([(nx+nz-22),nu+(nx+nz-22)+(nx+nz+nu):d*(nx+nz+sv)+nx+nu:n_w_i]);
    
    % implement the first sample on the simulator
    u_in_1  = u_opt1(1) ;
    u_in_2  = u_opt2(1) ;
    u_nlp   = [u_in_1;u_in_2;par.GOR];

end
