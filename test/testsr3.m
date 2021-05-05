%% unit tests for sr3 routine

% initialize

clear; clf; close all;
iseed = 8675309;
rng(iseed);
set(groot, 'defaultLineMarkerSize',10)
set(groot, 'defaultLineLineWidth',2)


%% synthetic problem

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
lam0 = 0.004; % good for l_0 regularizer

l0w = lam0*0.5; % weights for mixed regularizer
l1w = lam1*0.5;

% apply solver
[x0, w0] = sr3(A, b, 'mode', '0', 'lam',lam0,'ptf',0);
[x1, w1] = sr3(A, b, 'lam',lam1,'ptf',0);
[xmix,wmix] = sr3(A, b, 'l0w',l0w,'l1w',l1w,'mode','mixed');

% built-ins
xl2 = A\b;
if exist('lasso','builtin')
    xl1 = lasso(A,b,'Lambda',lam1);
end

% plot solution
% both regularizers perform well on this problem, though the $\ell_1$
% regularizer introduces a little more bias
figure(); hold on;
plot(y, '-*b'); plot(x0, '-xr'); plot(w0, '-or'); plot(x1, '-xc');
plot(w1, '-oc'); scatter(1:length(xl2),xl2,'ok', ...
    'MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.25);
plot(xmix,'-xg'); plot(wmix, '-og')
if exist('lasso','builtin')
    plot(xl1,'-xk'); 
    legend('true signal', 'x0', 'w0', 'x1', 'w1','backslash','xmix',...
        'wmix','lasso');
else
    legend('true signal', 'x0', 'w0', 'x1', 'w1','backslash','xmix',...
        'wmix');    
end
