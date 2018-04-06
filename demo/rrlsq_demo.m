%% Relaxed regularized least squares RRLSQ demo
%
% In this file we demonstrate the use of the |rrlsq| MATLAB(R) routine
% on a few examples in regularized least squares fitting.
%
% The RRLSQ framework is an approach to problems of the form 
%
% $$ \min_x \frac{1}{2} \|Ax-b\|_{\ell_2}^2 + \lambda \rho(D x) $$ 
%
% where $\rho(Dx)$ represents some regularizer. For example, in LASSO
% style regression, we have that $D$ is the identity and 
% $\rho(x)=\|x\|_{\ell_1}$.
%
% The relaxation comes in by introducing a variable $w$ which
% does a good job of minimizing the regularizer, while $x$ is allowed
% to vary a little around that. Consider instead the minimization
% problem
%
% $$ \min_{x,w} \frac{1}{2} \|Ax-b\|_{\ell_2}^2 + \lambda \rho(w)
%   + \frac{\kappa}{2} \|Dx-w\|_{\ell_2}^2 $$
%
% When $D$ is the identity, the variable $w$ tends to be a good
% approximation of the solution of the original problem. For sparse
% signal problems in particular, $w$ performs well at recovering
% the 'true' support. When $D$ is not the identity, the splitting
% into $x$ and $w$ allows for a simple prox-type operation to be
% applied to $w$ (for prox friendly $\rho$) as opposed to dealing with
% the more complicated penalty $\rho(Dx)$.

% initialize
clear; clf; close all;
iseed = 8675309;
rng(iseed);

%% Problem 1: $\ell_1$ vs $\ell_0$ penalties
%
% In this problem we test the performance of using RRLSQ to 
% recover a sparse signal. Here $\rho$ is either the $\ell_1$ 
% or $\ell_0$ penalty. We also plot the results of a least squares
% fit and the built in |lasso| function (if available). This is a 
% relatively easy problem and all of the regularizers perform well, 
% beating the standard least squares approach for obvious reasons.

% matrix dimensions
m = 200;
n = 1000;
k = 10; % number of non-zeros in true solution
sigma = 1e-1; % additive noise

A = randn(m,n);

y = zeros(n,1);
ind = randperm(n,k);
y(ind) = sign(randn(k,1));

b = A*y+sigma*randn(m,1);

% set up parameters of fit
lam1 = 0.01; % good for l_1 regularizer
lam0 = 0.002; % good for l_0 regularizer

% apply solver
[x0, w0] = rrlsq(A, b, 'mode', '0', 'lam',lam0,'ptf',0);
[x1, w1] = rrlsq(A, b, 'lam',lam1,'ptf',0);

% built-ins
xl2 = A\b;
if exist('lasso')
    xl1 = lasso(A,b,'Lambda',lam1);
end

% plot solution
% both regularizers perform well on this problem, though the $\ell_1$
% regularizer introduces a little more bias
figure(); hold on;
plot(y, '-*b'); plot(x0, '-xr'); plot(w0, '-og'); plot(x1, '-xc');
plot(w1, '-om'); scatter(1:length(xl2),xl2,'ok', ...
    'MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.25);
if exist('lasso')
    plot(xl1,'-ok'); 
    legend('true signal', 'x0', 'w0', 'x1', 'w1','backslash','lasso');
else
    legend('true signal', 'x0', 'w0', 'x1', 'w1','backslash');    
end

