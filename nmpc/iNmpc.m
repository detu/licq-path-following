function [Tall, xmeasureAll, uAll, ObjVal, primalNLP, params, runtime] = iNmpc(optProblem, system, mpciterations, N, T, tmeasure, xmeasure, u0, varargin)
%INMPC Summary of this function goes here
% 
% [OUTPUTARGS] = INMPC(INPUTARGS) Explain usage here
% 
% Ideal NMPC (single step)
%
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2016/04/25 19:24:15 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2016


Tall    = [];
Xall    = zeros(mpciterations,size(xmeasure,1));
Uall    = zeros(mpciterations,size(u0,1));

xmeasureAll = [];
uAll        = [];
runtime     = [];
x_nlp_opt   = [];
u_nlp_opt   = [];
mpciter = 1;

load noise1pct.mat;

global nx;

while(mpciter <= mpciterations)
    
    fprintf('-----------------------------\n');
    fprintf('MPC iteration: %d\n', mpciter);
    
    % obtain new initial value
    [t0, x0] = measureInitialValue ( tmeasure, xmeasure );
    
    % add measurement error to x0
    holdupNoise        = noise(:,mpciter);
    concentrationNoise = zeros(42,1);
    measNoise          = [concentrationNoise;holdupNoise];
    x0_measure         =  x0 + measNoise;
    

    % ideal NMPC:
    [primalNLP, ~, lb, ub, ~, params, elapsedtime] = solveOptimalControlProblem(optProblem, system, N, t0, x0, u0, T, mpciter, u_nlp_opt, x_nlp_opt, x0_measure);
    
    % re-arrange NLP results
    [u_nlp_opt, x_nlp_opt] = plotStatesN(primalNLP, lb, ub, N);  
    
    % save open loop solution for error computation
    %iNmpcData(mpciter).z = x_nlp_opt;

    z1 = x_nlp_opt(1:nx,5);  % 5 = (d+1) + 1 (d=number of collocation point)
    
    
    % record information
    Tall                = [Tall t0];
    Xall(mpciter+1,:)   = x0';
    Uall(mpciter+1,:)   = u0(:,1);
    
    
    % Apply control to process with optimized control from path-following
    % algorithm. 
    x0 = xmeasure;  % from the online step 
    [tmeasure, xmeasure] = applyControl(system, T, t0, x0, u_nlp_opt);
    
    ObjVal(mpciter) = computeObjectiveFunctionValues(u_nlp_opt(:,1),xmeasure); % USING ACTUAL STATE!
    
    
    % store output variables
    xmeasureAll = [xmeasureAll;xmeasure];
    uAll        = [uAll;u_nlp_opt(:,1)];   
    runtime     = [runtime;elapsedtime];
    
    % prepare restart
    u0 = shiftHorizon(u_nlp_opt);
    
    mpciter = mpciter+1;
end
xmeasureAll = reshape(xmeasureAll,84,mpciterations); 
%save iNmpcData.mat iNmpcData;
end

function [t0, x0] = measureInitialValue ( tmeasure, xmeasure )
    t0 = tmeasure;
    x0 = xmeasure;
end

function [tapplied, xapplied] = applyControl(system, T, t0, x0, u0)
    xapplied = dynamic(system, T, t0, x0, u0(:,1));
    tapplied = t0+T;
end

function u0 = shiftHorizon(u)
    u0 = [u(:,2:size(u,2)) u(:,size(u,2))];
end

function [u, lamda, lbw, ubw, objVal, params, elapsednlp] = solveOptimalControlProblem(optProblem, system, N, t0, x0, u0, T, mpciter, u_nlp, x_nlp, x0_measure)

    import casadi.*
    % call ODE15s N-times initial guess in optimization
    x(1,:) = x0';
    for k=1:N
        x(k+1,:) = x0';    % ONE-TIME SIMULATION !
    end

    [J,g,w0,w,lbg,ubg,lbw,ubw,params] = optProblem(x, u0, N, x0_measure);
    prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));
    options = struct;
    %options.ipopt.tol       = 1e-12;
    %options.acceptable_compl_inf_tol    = 1e-6; 
    solver = nlpsol('solver', 'ipopt', prob, options);
    
    % Solve the NLP
    startnlp   = tic;
    sol        = solver('x0', w0, 'lbx', lbw, 'ubx', ubw, 'lbg', lbg, 'ubg', ubg);
    elapsednlp = toc(startnlp);
    fprintf('IPOPT solver runtime = %f\n',elapsednlp);

    u      = full(sol.x);
    lamda  = full(sol.lam_g);
    objVal = full(sol.f);
end

function [x, t_intermediate, x_intermediate] = dynamic(system, T, t0, x0, u0)
    x = system(t0, x0, u0, T);
    x_intermediate = [x0; x];
    t_intermediate = [t0, t0+T];
end

