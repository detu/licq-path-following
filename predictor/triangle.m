function [prob] = triangle(p)
%TRIANGLE Summary of this function goes here
% 
% A case example from book Guddat et. al
% Ex. 3.4.2 page 87
% Testing of EQUALITY constraints 
% [OUTPUTARGS] = TRIANGLE(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2015/10/19 14:45:18 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2015

prob = struct('neq',{0},'niq',{0},'cin',{0},'ceq',{0},'dp_in',{0},'dp_eq',{0},'hess',{0},'lxp',{0},'x',0,'name',0);
prob.neq  = 3;            % number of equality constraint
prob.niq  = 0;            % number of inequality constraint
prob.name = 'Triangle';
prob.x    = zeros(5,1);

prob.cin   = (@(x,p)(constraint(x,p)));  % inequality constraint and prob.ceq for equality constraint
prob.ceq   = (@(x,p)(cons_equal(x,p)));
prob.dp_in = (@(x,p)(deriv_param_in(x,p)));
prob.dp_eq = (@(x,p)(deriv_param_eq(x,p)));
prob.hess  = (@(x,y,p)(hessian(x,y,p)));
prob.lxp   = (@(x,y,p)(lag_der_xp(x,y,p)));
prob.obj   = (@(x,p)(objective(x,p)));

end

function [f,g] = objective(x,p)
    f = -p*x(1)*x(2)*x(4) + (1-p)*( ( x(1) - 1.5 )^2 + x(2)^2 + ( x(3) - 1.5 )^2 + ( x(4) - 1 )^2 + x(5)^2 );
    g = zeros(5,1);
    g(1) = -p*x(2)*x(4) + (1-p)*2*(x(1)-1.5);
    g(2) = -p*x(1)*x(4) + (1-p)*2*x(2);
    g(3) = (1-p)*2*(x(3)-1.5);
    g(4) = -p*x(1)*x(4) + (1-p)*2*(x(4)-1);
    g(5) = (1-p)+2*x(5);
end

function [c,J] = constraint(x,p)
    c = [];
    J = [];
end

function [Jeq] = cons_equal(x,p)
    Jeq = zeros(3,5);
    Jeq(1,1)   = 1;
    Jeq(1,2)   = 1;
    Jeq(1,3)   = 1;
    Jeq(2,1)   = (2*x(1)) - (2*x(2)*x(5));
    Jeq(2,2)   = (2*x(2)) - (2*x(1)*x(5));
    Jeq(2,3)   = -2*x(3);
    Jeq(2,5)   = -2*x(1)*x(2);
    Jeq(3,4)   = 2*x(4);
    Jeq(3,5)   = 2*x(5);
end

function [cp] = deriv_param_in(x,p)
    cp = [];
end

function [cp] = deriv_param_eq(x,p)
    cp = zeros(3,1);
end

function [H] = hessian(x,y,p)
    H      = zeros(5,5);
    H(1,1) = 2*(1-p) + 2*y(2);
    H(1,2) = -p*x(4) - (2*y(2)*x(5));
    H(1,4) = -p*x(2);
    H(1,5) = -2*y(2)*x(2);
    H(2,1) = -p*x(4) - (2*y(2)*x(5));
    H(2,2) = 2*(1-p) + 2*y(2);
    H(2,4) = -p*x(1);
    H(2,5) = -2*y(2)*x(1);
    H(3,3) = 2*(1-p) - (2*y(2));
    H(4,1) = -p*x(2);
    H(4,2) = -p*x(1);
    H(4,4) = 2*(1-p) + 2*y(3);
    H(5,1) = -2*y(2)*x(2);
    H(5,2) = -2*y(2)*x(1);
    H(5,5) = 2*(1-p) + 2*y(3);
end

function [Lxp] = lag_der_xp(x,y,p)
    Lxp      = zeros(5,1);
    Lxp(1,1) = (-x(2)*x(4)) - (2*(x(1)-1.5));
    Lxp(2,1) = (-x(1)*x(4)) - (2*x(2));
    Lxp(3,1) = -2*(x(3) - 1.5);
    Lxp(4,1) = (-x(1)*x(2)) - (2*(x(4) - 1));
    Lxp(5,1) = -2*x(5);
end