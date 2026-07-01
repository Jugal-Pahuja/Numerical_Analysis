%% ========================================================================
%  LAB 9 — COMPOSITE NUMERICAL QUADRATURE
%  (Composite Trapezoidal Rule and Composite Simpson's Rule)
%  MTH 308: Numerical Analysis and Scientific Computing-I, IIT Kanpur
%  ------------------------------------------------------------------------
%  BIG PICTURE
%  Both rules approximate a definite integral by replacing the function
%  with a simple interpolant over each of n equal subintervals, then
%  summing the (exactly known) integral of each piece:
%    Trapezoidal : approximates f as PIECEWISE LINEAR -> error O(h^2)
%    Simpson's   : approximates f as PIECEWISE QUADRATIC -> error O(h^4),
%                  so it should noticeably out-perform Trapezoidal for
%                  the same number of subintervals n (provided f is
%                  smooth enough).
%  ========================================================================

%% Problem 1 — Approximate the integral of 1/(x ln x) from e to e+2
% Exact value: ln(ln(e+2))  [since d/dx[ln(ln x)] = 1/(x ln x)]
% We use n = 6 subintervals for both rules and compare accuracy.

n = 6;
x = linspace(exp(1), exp(1)+2, n+1);    % n+1 equally-spaced nodes
f = @(x) 1./(x.*log(x));

T = CompositeTrapezoidalRule(x, f(x));
S = CompositeSimpsonsRule(x, f(x));

I_exact = log(log(exp(1)+2));

fprintf('-------- Problem 1: integral of 1/(x ln x) from e to e+2, n=%d --------\n\n', n);
fprintf('Composite Trapezoidal Rule: %.6f\n', T)
fprintf('Composite Simpson''s Rule: %10.6f\n', S)
fprintf('Exact Value: %10.6f\n\n', I_exact)

fprintf('Absolute Error (Trapezoidal): %.6e\n', abs(T - I_exact))
fprintf('Absolute Error (Simpson):     %.6e\n', abs(S - I_exact))
% Expect Simpson's error to be several orders of magnitude smaller than
% Trapezoidal's, consistent with their O(h^4) vs O(h^2) error orders.


%% =======================================================================
%  FUNCTIONS
%  =======================================================================

function approx = CompositeTrapezoidalRule(x, f)
% COMPOSITETRAPEZOIDALRULE  Approximates the integral of f over [x(1),
% x(end)] using the composite trapezoidal rule:
%   integral ~ h/2 * [ f(x0) + 2*sum(f(x1..x_{n-1})) + f(xn) ]
%
%   INPUT
%     x : vector of n+1 evenly-spaced nodes
%     f : vector of function values at those nodes
    h = (x(end)-x(1))/(length(x)-1);
    approx = h/2*(f(1)+2*sum(f(2:end-1))+f(end));
end


function approx = CompositeSimpsonsRule(x, f)
% COMPOSITESIMPSONSRULE  Approximates the integral of f over [x(1),
% x(end)] using the composite Simpson's rule (n must be even):
%   integral ~ h/3 * [ f(x0) + 4*sum(odd-indexed f) + 2*sum(even-indexed
%   interior f) + f(xn) ]
%
%   INPUT
%     x : vector of n+1 evenly-spaced nodes (n even)
%     f : vector of function values at those nodes
    h = (x(end)-x(1))/(length(x)-1);
    approx = 1/3*h*(f(1)+4*sum(f(2:2:end-1))+2*sum(f(3:2:end-1))+f(end));
end
