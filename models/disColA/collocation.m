function F = collocation(dae, tf, nsteps)
% Joel Andersson, joel@casadi.org, 2016
import casadi.*

% Degree of interpolating polynomial
d = 3;

% Get collocation points
tau_root = collocationPoints(d, 'legendre');

% Coefficients of the collocation equation
C = zeros(d+1,d+1);

% Coefficients of the continuity equation
D = zeros(d+1, 1);

% Coefficients of the quadrature function
B = zeros(d+1, 1);

% Construct polynomial basis
for j=1:d+1
  % Construct Lagrange polynomials to get the polynomial basis at the collocation point
  coeff = 1;
  for r=1:d+1
    if r ~= j
      coeff = conv(coeff, [1, -tau_root{r}]);
      coeff = coeff / (tau_root{j}-tau_root{r});
    end
  end
  % Evaluate the polynomial at the final time to get the coefficients of the continuity equation
  D(j) = polyval(coeff, 1.0);

  % Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
  pder = polyder(coeff);
  for r=1:d+1
    C(j,r) = polyval(pder, tau_root{r});
  end

  % Evaluate the integral of the polynomial to get the coefficients of the quadrature function
  pint = polyint(coeff);
  B(j) = polyval(pint, 1.0);  
end

% Continuous time dynamics
f = Function('f', {dae.x, dae.z}, {dae.ode, dae.alg});

% Variables for one finite element
X0 = MX.sym('X0', dae.x.sparsity());
Z0 = MX.sym('Z0', dae.z.sparsity());
X = {};
Z = {};
V = {};
for k=1:d
    X{k} = MX.sym(['X' num2str(k)], dae.x.sparsity());
    V{end+1} = X{k};
    Z{k} = MX.sym(['Z' num2str(k)], dae.z.sparsity());
    V{end+1} = Z{k};
end

% Collocation equations and quadrature
g = {};
xf = D(1)*X0;
zf = D(1)*Z0;
h = tf/nsteps;
for j=1:d
  % Expression for the state derivative at the collocation point
  xp = C(1,j+1)*X0;
  for r=1:d
    xp = xp + C(r+1,j+1)*X{r};
  end
      
  % Append collocation equations
  [fj, aj] = easycall(f, X{j},Z{j});
  g{end+1} = h*fj - xp;

  % Append algebraic equations
  g{end+1} = aj;
  
  % Add contribution to the end state
  xf = xf + D(j+1)*X{j};
  zf = zf + D(j+1)*Z{j};
end

% Rootfinder function, implicitly defines X and Q as a function of X0 and U
rfp = Function('rfp', {vertcat(V{:}), X0, Z0},...
                      {vertcat(g{:}), xf, zf});                  
rf = rootfinder('rf', 'newton', rfp);

% Discrete time dynamics
XF = X0;
ZF = Z0;
for j=1:nsteps
    [~, XF, ZF] = easycall(rf, repmat([XF;ZF], d, 1), XF, ZF);
end
F = Function('F', {X0, Z0}, {XF, ZF}, ...
    char('x0', 'z0'), char('xf', 'zf'));
