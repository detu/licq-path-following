function [prob] = kyparisis(p)
%KYPARISIS Summary of this function goes here
% 
% [prob] = KYPARISIS(p) Explain usage here
% a case example from Kyparisis's paper
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2015/10/05 14:27:54 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2015

prob = struct('neq',{0},'niq',{0},'cin',{0},'ceq',{0},'dp_in',{0},'dp_eq',{0},'hess',{0},'lxp',{0},'x',0,'name',0);
prob.neq  = 0;            % number of equality constraint
prob.niq  = 3;            % number of inequality constraint
prob.name = 'Kyparisis';
prob.x    = zeros(2,1);

prob.cin   = (@(x,p)(constraint(x,p)));  % inequality constraint and prob.ceq for equality constraint
prob.ceq   = (@(x,p)(cons_equal(x,p)));
prob.dp_in = (@(x,p)(deriv_param_in(x,p)));
prob.dp_eq = (@(x,p)(deriv_param_eq(x,p)));
prob.hess  = (@(x,y,p)(hessian(x,y,p)));
prob.lxp   = (@(x,y,p)(lag_der_xp(x,y,p)));

end

function [c,J] = constraint(x,p)
    c      = [-x(1) ; -x(2) ; -x(1)-x(2)+p(1)+p(2)];
    J      = zeros(3,2);
    J(1,1) = -1;
    J(2,2) = -1;
    J(3,1) = -1;
    J(3,2) = -1;
end

function [Jeq] = cons_equal(x,p)
    Jeq = [];
end

function [cp] = deriv_param_in(x,p)
    cp      = zeros(3,2);
    cp(3,1) = 1;
    cp(3,2) = 1;
end

function [cp] = deriv_param_eq(x,p)
    cp = [];
end

function [H] = hessian(x,y,p)
    H      = zeros(2,2);
    H(1,1) = 1;
    H(2,2) = 1;
end

function [Lxp] = lag_der_xp(x,y,p)
    Lxp      = zeros(2);
    Lxp(1,1) = -1;
    Lxp(2,2) = -1;
end