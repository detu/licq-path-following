function [prob] = prob4_casadi(p)
%PROB4_CASADI Summary of this function goes here
% 
% Modification of prob4 with CasADi automatic differentation tool
% [OUTPUTARGS] = PROB4_CASADI(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: detu $	$Date: 2015/11/15 17:28:48 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2015

% add CasADi 
% addpath('C:\Users\detu\Documents\casadi-matlabR2014b-v2.4.1')
import casadi.*

prob = struct('neq',{0},'niq',{0},'cin',{0},'ceq',{0},'dp_in',{0},'dp_eq',{0},'hess',{0},'lxp',{0},'x',0,'name',0);
% prob.neq  = 0;            % number of equality constraint
% prob.niq  = 2;            % number of inequality constraint
prob.neq  = 2;            % number of equality constraint
prob.niq  = 0;            % number of inequality constraint
prob.name = 'Problem 4';
prob.x    = zeros(2,1);

prob.obj   = (@(x,y,p)(objective(x,y,p)));

end


function [f,g,H,Lxp,cst,J,cp,Jeq,dpe] = objective(x,y,p)
    import casadi.*
   
    % construct Lagrangian
    v = SX.sym('v',2);  % x variable resemble
    t = SX.sym('t',2);  % p variable resemble 
    f        = t(1)*v(1)^3+v(2)^2;
    c_expr   = [exp(-v(1))-v(2);t(2)-v(1)];
    lag_expr = f + y'*c_expr;
    
    % objective function related derivatives
    %opf = struct;
    %opf.input_scheme = char('v','t');
    %opf.output_scheme = char('f');
    %f = SXFunction('f',{v,t},{f},opf);
    %g = SXFunction(f.gradient);
    f = Function('f',{v,t},{f},char('v', 't'), char('f'));
    g = f.gradient();
    
    % Lagrangian and constraint related derivatives
    %op       = struct;
    %op.input_scheme  = char('v','t');
    %op.output_scheme = char('lag_expr');
    %lagr       = SXFunction('lagr',{v,t},{lag_expr},op);
    %H          = SXFunction(lagr.hessian('v','lag_expr'));
    %jac_lagr_x = SXFunction(lagr.jacobian('v','lag_expr'));
    %Lxp        = jac_lagr_x.jacobian(1,0);
    lagr       = Function('lagr',{v,t},{lag_expr},char('v','t'),char('lag_expr'));
    H          = Function(lagr.hessian('v','lag_expr'));
    jac_lagr_x = Function(lagr.jacobian('v','lag_expr'));
    Lxp        = jac_lagr_x.jacobian(1,0);
    
    % Inequality constraint
    %opc = struct;
    %opc.input_scheme  = char('v','t');
    %opc.output_scheme = char('c_expr');
    %cst = SXFunction('cst',{v,t},{c_expr},opc);
    %J   = cst.jacobian('v','c_expr');
    %cp  = cst.jacobian('t','c_expr');
    cst = Function('cst', {v,t}, {c_expr}, char('v','t'), char('c_expr'));
    J   = cst.jacobian('v','c_expr');
    cp  = cst.jacobian('t','c_expr');
    
    % evaluate and obtain their values
%     f.setInput(x,'v');
%     f.setInput(p,'t');
%     g.setInput(x,'v');
%     g.setInput(p,'t');
%     H.setInput(x,'v');
%     H.setInput(p,'t');
%     Lxp.setInput(x,'v');
%     Lxp.setInput(p,'t');
%     J.setInput(x,'v');
%     J.setInput(p,'t');
%     cp.setInput(x,'v');
%     cp.setInput(p,'t');
%     cst.setInput(x,'v');
%     cst.setInput(p,'t');
    
    f   = f(x,p);
    g   = g(x,p);
    H   = H(x,p);
    Lxp = Lxp(x,p);
    J   = J(x,p);
    cp  = cp(x,p);
    cst = cst(x,p);

%     f.evaluate();
%     g.evaluate();
%     H.evaluate();   
%     J.evaluate();
%     cp.evaluate();
%     Lxp.evaluate();
%     cst.evaluate();
%     
%     f   = full(f.getOutput());
%     g   = sparse(g.getOutput());
%     H   = sparse(H.getOutput());
%     Lxp = sparse(Lxp.getOutput());
%     J   = sparse(J.getOutput());
%     cp  = sparse(cp.getOutput());
%     cst = full(cst.getOutput());
%     % Equality constraint
%     Jeq = [];
%     dpe = [];

    f   = full(f);
    g   = sparse(g);    
    H   = sparse(H);
    Lxp = sparse(Lxp);
    J   = sparse(J);
    cp  = sparse(cp);
    cst = full(cst);

    % Equality constraint
    %Jeq = [];
    %dpe = [];
    Jeq = J;
    dpe = cp;
end
