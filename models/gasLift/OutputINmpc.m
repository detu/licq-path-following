function [Tall, xmeasureAll, uAll, ObjVal, primalNLP, params, runtime] = OutputINmpc(optProblem, system, mpciterations, N, T, tmeasure, xmeasure, u0, paramModel, varargin)
%OUTPUTINMPC Summary of this function goes here
% 
% [OUTPUTARGS] = OUTPUTINMPC(INPUTARGS) Explain usage here
% 
% Output ideal NMPC consist of NMPC and Observer (EKF)
%
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2017/10/25 12:37:55 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2017

nx      = numel(xmeasure);
nu      = size(u0,1);
Tall    = [];
Xall    = zeros(mpciterations,size(xmeasure,1));
Uall    = zeros(mpciterations,size(u0,1));

xmeasureAll = [];
uAll        = [];
runtime     = [];
x_nlp_opt   = [];
u_nlp_opt   = [];
ObjVal      = [];
mpciter = 1;

load noise1pct.mat;

%global nx;

while(mpciter <= mpciterations)
    
    fprintf('-----------------------------\n');
    fprintf('MPC iteration: %d\n', mpciter);
    
    % obtain new initial value
    [t0, x0] = measureInitialValue ( tmeasure, xmeasure );
    
    % add measurement error to x0
    measNoise          = noise(1:nx,mpciter);
    %measNoise          = [noise(1:8,mpciter);zeros(30,1)];
    x0_measure         =  x0 + 0*measNoise;     

    % ideal NMPC:
    [primalNLP, ~, lb, ub, ~, params, elapsedtime] = solveOptimalControlProblem(optProblem, system, N, t0, x0, u0, T, mpciter, u_nlp_opt, x_nlp_opt, x0_measure, paramModel);
    
    % re-arrange NLP results
    %[u_nlp_opt, x_nlp_opt] = plotStatesN(primalNLP, lb, ub, N);  
    %[u_nlp_opt, x_nlp_opt] = processNLPResult(primalNLP, paramModel);  
    [u_nlp_opt, x_nlp_opt] = plotStatesGL(primalNLP, lb, ub, N);
    
    % save open loop solution for error computation
    %iNmpcData(mpciter).z = x_nlp_opt;

    %z1 = x_nlp_opt(1:nx,5);  % 5 = (d+1) + 1 (d=number of collocation point)
    
    
    % record information
    Tall                = [Tall t0];
    Xall(mpciter+1,:)   = x0';
    %Uall(mpciter+1,:)   = u0(:,1);
    Uall(mpciter+1,:)   = u0(1:2);
    
    
    % Apply control to process with optimized control from path-following
    % algorithm. 
    x0 = xmeasure;  % from the online step 
    %x0                   = x0_measure; 
    [tmeasure, xmeasure] = applyControl(system, T, t0, x0, u_nlp_opt, paramModel);
    
    %ObjVal(mpciter) = computeObjectiveFunctionValues(u_nlp_opt(:,1),xmeasure); % USING ACTUAL STATE!
    
    
    % store output variables
    xmeasureAll = [xmeasureAll;xmeasure];
    %uAll        = [uAll;u_nlp_opt(:,1)];   
    uAll        = [uAll;u_nlp_opt];
    runtime     = [runtime;elapsedtime];
    
    % prepare restart
    u0 = shiftHorizon(u_nlp_opt);
    
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
    u      = full(sol.x);
    lamda  = full(sol.lam_g);
    objVal = full(sol.f);
    
    
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

