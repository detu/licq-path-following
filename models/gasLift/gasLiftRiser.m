function prob = gasLiftRiser(p,paramModel)
%GASLIFTRISER Summary of this function goes here
% 
% [OUTPUTARGS] = GASLIFTRISER(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2017/10/28 02:15:25 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2017

import casadi.*

prob = struct('neq',{0},'niq',{0},'cin',{0},'ceq',{0},'dp_in',{0},'dp_eq',{0},'hess',{0},'lxp',{0},'x',0,'name',0);
prob.neq  = 2000;         % HARD-CODE !    % number of equality constraint
prob.niq  = 0;            % number of inequality constraint
prob.name = 'Gas Lift Riser';
prob.x    = zeros(2,1);

prob.obj  = (@(x,y,p,N,paramModel)(objective(x,y,p,N,paramModel)));

end

function [f,g,H,Lxp,cst,J,cp,Jeq,dpe,Hobj,bIndex] = objective(x,y,p,N,par)
    import casadi.*

    nPrimal = numel(x);
    nDual   = numel(y.lam_g);
    nParam  = numel(p);
    nu = 2;
    nz = 30;
    nx = 8;
    bIndex = [];

    % nomimal model
    sf = nomModel(par);

    % Direct Collocation
    [B,C,D,d] = collocationSetup();

    %Build NLP
    V    = {};
    obj  = 0;
    cons = {};

    % initial conditions for each scenario
    X0  = MX.sym('X0',nx);
    Z0  = MX.sym('Z0',nz);
    V0  = [X0;Z0];
    V   = {V{:}, X0,Z0};
    cons    = {cons{:}, [X0;Z0] - x(1:nx+nz,1)};
    cons_x0 = [X0;Z0] - x(1:nx+nz,1);
    bIndex  = [bIndex;ones(nx+nz,1)];

    % Formulate NLP
    Xk  = X0;
    Xkj = {};
    Zkj = {};

    load QmaxGL.mat;
    load xuSS.mat;
    
    for k = 0:par.N-1

        Uk  = MX.sym(['U_' num2str(k)],nu);
        V   = {V{:},Uk};
        
        Jcontrol   = (QmaxGL(nx+nz+1:nx+nz+nu,1).*(Uk - uSS))' * (Uk - uSS);
        
        Xkj = {};
        Zkj = {};

        for j = 1:d
            Xkj{j} = MX.sym(['X_' num2str(k) '_' num2str(j)],nx);
            Zkj{j} = MX.sym(['Z_' num2str(k) '_' num2str(j)],nz);
            V      = {V{:},Xkj{j},Zkj{j}};
        end

        % Loop over collocation points
        Xk_end  = D(1)*Xk;

        for j = 1:d
            % Expression for the state derivative at the collocation point
            xp  = C(1,j+1)*Xk;  % helper state
            for r = 1:d
                xp = xp + C(r+1,j+1)*Xkj{r};
            end
            [fj,zj,qj] =  sf(Xkj{j},Zkj{j},vertcat(Uk,par.GOR));

            cons    = {cons{:},par.tf*fj-xp,zj};  % dynamics and algebraic constraints
            bIndex  = [bIndex;ones(nx+nz,1)];

            % Add contribution to the end states
            Xk_end  = Xk_end + D(j+1)*Xkj{j};
            
        end

        % New NLP variable for state at end of interval
        Xk      = MX.sym(['X_' num2str(k+1) ], nx);
        V       = {V{:},Xk};
        
        % Shooting Gap constraint
        cons   = {cons{:},Xk_end-Xk};
        bIndex = [bIndex;ones(nx,1)];
        
%         % Gas capacity constraints
%         cons   = {cons{:},sum(Zkj{j}(17:18))};
%         bIndex  = [bIndex;0];
        
%         % Gas Lift constraint
%         cons   = {cons{:},sum(Uk)};
%         bIndex  = [bIndex;0];
        
        Jstate =(QmaxGL(1:nx+nz,1).*([Xk;Zkj{j}] - xSS))' * ([Xk;Zkj{j}] - xSS);
        Jecon  = -Zkj{j}(29) + sum(Uk);
        obj    = obj + Jcontrol + Jstate + Jecon;
        
    end
    

    
    V = vertcat(V{:});
    % Concatenate constraints
    cons  = vertcat(cons{:});
    
    % objective function and constraint functions
    f = Function('f', {V}, {obj}, char('V'), char('objective'));
    c = Function('c', {V}, {cons}, char('V'), char('constraint'));
    %cx0 = Function('cx0', {X0}, {cons_x0}, char('X0'), char('constraint'));
    cx0 = Function('cx0', {V0}, {cons_x0}, char('V0'), char('constraint'));
    
    % construct Lagrangian
    lag_expr = obj + y.lam_g'*cons;
    
    g    = f.gradient();
    lagr = Function('lagr', {V}, {lag_expr}, char('V'), char('lag_expr'));
    H    = Function(lagr.hessian('V','lag_expr'));
    Hobj = f.hessian('V','objective');
    J    = c.jacobian('V','constraint');
    Jp   = cx0.jacobian('V0','constraint');
    
    f   = f(x); 
    g   = g(x);
    H   = H(x);
    Lxp = H(1:nPrimal,1:nParam);
    J   = J(x);
    %cp  = J(1:nDual,1:nParam);
    Jtemp  = zeros(nDual,nParam);
    cp  = Jp(x(1:nParam));
    Jtemp(1:nParam,1:nParam) = full(cp);
    cp  = sparse(Jtemp);
    cst = c(x);
    
    % Evaluation of objective function used for Greshgorin bound
    Hobj = Hobj(x);
    Hobj = sparse(Hobj);
    
    f   = full(f);
    g   = sparse(g);
    H   = sparse(H);
    Lxp = sparse(Lxp);
    J   = sparse(J);
    cp  = sparse(cp);
    cst = full(cst);
    
    % Equality constraint
    Jeq = J;
    dpe = cp;

end
