%% ========================================================================
%  LAB 10 — GAUSS-LEGENDRE QUADRATURE vs. CLENSHAW-CURTIS QUADRATURE
%  MTH 308: Numerical Analysis and Scientific Computing-I, IIT Kanpur
%  ------------------------------------------------------------------------
%  BIG PICTURE
%  Unlike composite Trapezoidal/Simpson rules (equally-spaced nodes,
%  fixed weights), both methods here choose nodes and weights to make an
%  N-point rule as accurate as possible:
%    Gauss-Legendre  : nodes = roots of the Legendre polynomial P_N;
%                      EXACT for polynomials up to degree 2N-1 ("degree
%                      of exactness" 2N-1).
%    Clenshaw-Curtis : nodes = Chebyshev points (cos(k*pi/N)); EXACT for
%                      polynomials up to degree N (a lower degree of
%                      exactness than Gauss-Legendre at the same N, but
%                      computed via an FFT-friendly, numerically stable
%                      formula and often nearly as accurate in practice
%                      for smooth functions).
%
%  We test both rules, for N = 1..50, on four integrals over [-1,1] with
%  very different smoothness properties, to see how integrand smoothness
%  affects each method's convergence rate.
%  ========================================================================

Nmax = 50;

% --- exact values of the four test integrals on [-1,1] ---
Itrue = [ ...
    1/2;                                        % I1 = integral of |x|^3
    2*(exp(-1) + sqrt(pi)*(erf(1)-1));          % I2 = integral of exp(-x^-2)
    pi/2;                                       % I3 = integral of 1/(1+x^2)
    2/11                                        % I4 = integral of x^10
];

labels = {'|x|^3','exp(-x^{-2})','1/(1+x^2)','x^{10}'};

% --- storage for the absolute error of each method, each function, ---
% --- and each value of N ---
E_cc = zeros(4, Nmax);   % Clenshaw-Curtis errors
E_ga = zeros(4, Nmax);   % Gauss-Legendre errors


%% Problem 1 — Compute both quadrature rules for N = 1..50, all 4 functions
for N = 1:Nmax

    [x_cc, w_cc] = clencurt(N);
    [x_ga, w_ga] = gauss(N);

    F_cc = test_functions(x_cc);
    F_ga = test_functions(x_ga);

    for k = 1:4
        E_cc(k,N) = abs(w_cc * F_cc(k,:)' - Itrue(k));
        E_ga(k,N) = abs(w_ga * F_ga(k,:)' - Itrue(k));
    end
end


%% Problem 2 — Plot error vs. N for both methods (semilog scale)
plot_errors(E_cc, labels, Nmax, 'Clenshaw–Curtis Quadrature Errors');
plot_errors(E_ga, labels, Nmax, 'Gauss–Legendre Quadrature Errors');


%% Problem 3 — Convergence behavior and the effect of integrand smoothness
% OBSERVATIONS (inspect the plots above to confirm):
%   - 1/(1+x^2) and exp(-x^-2) are SMOOTH (infinitely differentiable) on
%     [-1,1], so both quadrature rules converge geometrically/spectrally
%     fast: the error drops by many orders of magnitude as N grows, until
%     it bottoms out near machine precision (~1e-16).
%   - x^10 is a polynomial, so once N is large enough to exactly
%     represent a degree-10 polynomial, the error drops to ~machine
%     precision and STAYS there (see Problem 4 below for exactly how
%     large N must be for each method).
%   - |x|^3 is only C^2 (its third derivative is discontinuous at x=0),
%     so its convergence is much SLOWER and only algebraic (a fixed
%     negative power of N) rather than spectral — the error curve looks
%     close to a straight line on a log-log plot rather than a sharp drop
%     on the semilog plot used here. This is the expected behavior:
%     quadrature accuracy for non-smooth integrands is limited by how
%     many continuous derivatives the function actually has.


%% Problem 4 — Degree of exactness test, using f(x) = x^10
% Gauss-Legendre with N nodes is theoretically EXACT for polynomials up
% to degree 2N-1. Since x^10 has degree 10, we expect the error to drop
% to (near) machine precision once 2N-1 >= 10, i.e. N >= 5.5, so N = 6.
fprintf('\nDegree of Exactness Test for Gauss–Legendre (f(x) = x^10):\n')
for N = 1:6
    [x,w] = gauss(N);
    approx = w * (x.^10);
    err = abs(approx - 2/11);
    fprintf('N = %d, Error for x^{10} = %.2e\n', N, err);
end
% Watch the error column: it should remain non-negligible for N=1..5
% (since 2N-1 < 10) and then collapse to ~1e-15/1e-16 at N=6 (since
% 2*6-1 = 11 >= 10), confirming the degree-of-exactness formula.


%% =======================================================================
%  FUNCTIONS
%  =======================================================================

function F = test_functions(x)
% TEST_FUNCTIONS  Evaluates all four test integrands at every node in x.
%   Row 1: |x|^3            (limited smoothness -> slow convergence)
%   Row 2: exp(-x^-2)       (smooth but flat/singular-looking near 0)
%   Row 3: 1/(1+x^2)        (smooth, classic Runge-type test function)
%   Row 4: x^10              (polynomial -> finite degree of exactness)
    F = zeros(4,length(x));

    F(1,:) = abs(x).^3;

    xtmp = x;
    xtmp(xtmp==0) = eps;     % avoid 0^-2 = Inf; exp(-large) underflows to 0 anyway
    F(2,:) = exp(-xtmp.^(-2));

    F(3,:) = 1./(1 + x.^2);

    F(4,:) = x.^10;
end


function plot_errors(E, labels, Nmax, title_str)
% PLOT_ERRORS  Plots absolute error vs. N (semilog scale) for each of
% the four test functions, in a 2x2 grid of subplots.
    figure('Position',[100 100 900 900])
    sgtitle(title_str,'FontSize',18,'FontWeight','bold')

    for k = 1:4
        subplot(2,2,k)
        % "+1e-100" avoids a log(0) warning if the error ever hits zero
        % to within floating-point precision
        semilogy(1:Nmax, E(k,:) + 1e-100,'o-','LineWidth',1.2)
        grid on

        xlim([1 Nmax])
        ylim([1e-18 1e+1])

        xlabel('N','FontSize',14)
        ylabel('Error','FontSize',14)
        title(labels{k},'FontSize',15,'FontWeight','bold')
    end
end


function [x,w] = clencurt(N)
% CLENCURT  Clenshaw-Curtis quadrature nodes and weights on [-1,1], using
% N+1 Chebyshev points x_j = cos(j*pi/N), j = 0..N.
% (Standard closed-form construction; see Trefethen, "Spectral Methods
% in MATLAB", for the derivation.)
    theta = pi*(0:N)'/N;
    x = cos(theta);

    w = zeros(1,N+1);
    ii = 2:N;
    v = ones(N-1,1);

    if mod(N,2)==0
        w(1) = 1/(N^2 - 1);
        w(N+1) = w(1);
        for k = 1:(N/2 - 1)
            v = v - 2*cos(2*k*theta(ii)) / (4*k^2 - 1);
        end
        v = v - cos(N*theta(ii))/(N^2 - 1);
    else
        w(1) = 1/N^2;
        w(N+1) = w(1);
        for k = 1:((N-1)/2)
            v = v - 2*cos(2*k*theta(ii)) / (4*k^2 - 1);
        end
    end

    w(ii) = 2*v/N;
end


function [x,w] = gauss(N)
% GAUSS  Gauss-Legendre quadrature nodes and weights on [-1,1], via the
% Golub-Welsch algorithm: the nodes are eigenvalues of a symmetric
% tridiagonal "Jacobi matrix" built from the three-term recurrence of
% the Legendre polynomials, and the weights come from the eigenvectors.
    beta = .5 ./ sqrt(1 - (2*(1:N-1)).^(-2));
    T = diag(beta,1) + diag(beta,-1);

    [V,D] = eig(T);
    x = diag(D);
    [x,ind] = sort(x);

    w = 2 * V(1,ind).^2;
end
