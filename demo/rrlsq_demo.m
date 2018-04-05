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
% $\rho(x)=\|x\|_{\ell_1}$

% initialize
clear; clf; close all;
iseed = 8675309;
rng(iseed);

%% Problem 1: l_1 vs l_0 penalties
%
% In this problem we test the performance of using RRLSQ to 
% recover a sparse signal. Here $\rho$ is either the $\ell_1$ 
% or $\ell_0$ penalty. We also plot the results of a least squares
% fit and the built in |lasso| function. This is a relatively easy
% problem and all of the regularizers perform well, beating the 
% standard least squares approach for obvious reasons.

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
[x1, w1] = rrlsq(A, b, 'mode', '1', 'lam',lam1,'ptf',0); 

% built-ins
xl2 = A\b;
xl1 = lasso(A,b,'Lambda',lam1);

% plot solution
% both regularizers perform well on this problem, though the $\ell_1$
% regularizer introduces a little more bias
figure(); hold on;
plot(y, '-*b'); plot(x0, '-xr'); plot(w0, '-og'); plot(x1, '-xc');
plot(w1, '-om'); plot(xl1,'-ok'); scatter(1:length(xl2),xl2,'ok', ...
    'MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.25);
legend('true signal', 'x0', 'w0', 'x1', 'w1','lasso','backslash');

%% Problem 2: edge detection with l_0, leveraging sparse routines
% 
% Because of the flexibility in choosing D, it is possible to
% create a simple-minded edge detection system using the $\ell_0$ penalty.
% In this example, we load an image, corrupt it with noise, define
% an edge-detection penalty, and use our solver with the option
% |ifusenormal| (this uses the normal equations and a Cholesky 
% factorization instead of a QR decomposition for the least squares 
% solves within |rrlsq|, which has significant advantages for a fit 
% with sparse A and D matrices). 
%
% We set $D$ to map to centered approximations of the $x$ and $y$ 
% derivatives of the image and $A$ is simply the identity. The non-zero
% entries in $w$ then correspond to 'edges'.

% load the popular 'cameraman' image
img = imread('cameraman.tiff');

% set up as a regularized least squares problem
% doing everything sparse for efficiency
b = double(img);
[m,~] = size(b); % image is square
b = b(:);
% corrupt with noise
sigma = 10;
b_w_noise = b + sigma*randn(size(b));
% w is the x and y derivatives of the image stacked on each other
% we simply leave out the pixels on the border here...
e = ones(m-2,1); diff = spdiags([-[e;0;0],[0;0;e]],[-1,1],m,m);
dmat = [kron(diff,speye(m)); kron(speye(m),diff)];
% x should approximate original image
amat = speye(m*m); 

lam0 = 500;

[x,w] = rrlsq(amat,b_w_noise,'D',dmat,'mode','0','lam',lam0,'ifusenormal',1);

ww = abs(w(1:end/2))+abs(w(end/2+1:end)); % get sum of x and y derivatives

% plot edges over noisy image

figure
h=pcolor((reshape(b_w_noise,m,m)));
set(h,'EdgeColor','none');
colormap('gray');
hold on
spy(reshape(ww,m,m));
