function [xnew, lam1, mu, lam2, fold, f1, f2, cc] = ESIIO(...
    m, n, loop, xval, xmin, xmax, f0val, df0dx, fval, dfdx, ...
    move, p, lam1, mu, lam2, fold, f1, f2, cc)
%
%    This is the file ESIIO.m
%
%    Written in May 2026 by
%    Zeng Meng <mengz@hfut.edu.cn>
%    School of Civil Engineering, Hefei University of Technology.
%    Hefei 230009, PR China.
%
%*** INPUT:
%
%      m  = Number of constraint functions.
%      n  = Number of design variables.
%   loop  = Number of iterations.
%   xval  = Design variables at k-th iterative step.
%   xmin  = Lower bounds of design variables.
%   xmax  = Upper bounds of design variables.
%  f0val  = Objective function value at current iteration.
%  df0dx  = Objective function sensitivity.
%   fval  = Constraint functions value at current iteration.
%   dfdx  = Constraint function sensitivity.
%   move  = The moving step size.
%      p  = Power of p-norm.
%   lam1  = Lagrange multipliers from the previous iteration.
%     mu  = Quadratic penalty factor from the previous iteration.
%   lam2  = Chaos control factor from the previous iteration.
%   fold  = Constraint function value from previous iteration.
%     f1  = Objective function value from one iterations ago.
%     f2  = Objective function value from two iterations ago.
%     cc  = Number of consecutive boundary crossings.
%     
%*** OUTPUT:
%
%   xnew  = Updated design variables.
%   lam1  = Updated Lagrange multipliers.
%     mu  = Updated quadratic penalty factor.
%   lam2  = Updated chaos control factor.
%   fold  = Constraint function value from current iteration.
%     f1  = Objective function value from current iteration.
%     f2  = Objective function value from previous iteration.
%     cc  = Updated number of consecutive boundary crossings.
%

% Constraint interval transformation method
x_mu = (xmax(:) + xmin(:)) / 2;
sigma = (xmax(:) - xmin(:)) / 2;
uval = (xval(:) - x_mu) ./ sigma;
umax = min(1, uval + 2 * move); umin = max(-1, uval - 2 * move);
% Explicit stable interval iteration
dLdx = df0dx(:) + (1/m) * dfdx' * (lam1(:) + mu(:).*fval(:));
df0du = dLdx .* sigma;
lam = norm(df0du, p/(p-1));
unew = -sign(df0du).*(abs(df0du)/lam).^(1/(p-1));
xnew = sigma .* max(umin, min(umax, unew)) + x_mu;
lam1(:) = lam1(:) + mu(:) .* max(fval(:), -lam1(:)./mu(:));
mu(:) = min(5000, mu(:) .* 1.02);
% Directional stability transformation method
if loop > 10
    cc = any(fval(:) .* fold(:) < 0) * (cc + 1);
    if cc >= 5
        lam2 = max(0.1, 0.95 * lam2);
        cc = min(cc, 4);
    elseif ((f1-f2)*(f0val-f1)>0 && all(abs(fval(:))<=1e-3))
        lam2 = min(1.0, 1.05 * lam2);
    end
    xnew = xval + lam2 * (xnew - xval);
    xnew = xnew/norm(xnew, inf);
else
    xnew = xval + lam2 * (xnew - xval);
end
f2 = f1;
f1 = f0val;
fold = fval;
