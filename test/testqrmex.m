
m = 3000;
n = 200;
nrhs = 10;

opts.UT = true;

A = randn(m,n) + 1i*randn(m,n);
xtrue = randn(n,nrhs) + 1i*randn(n,nrhs);

B = A*xtrue;

tol = 20*eps(1);
fprintf('Testing LAPACK based QR with tolerance = %e\n',tol)


%% test 1: complex least squares

[AOUT,jpvt,tau] = xgeqp3_m(A);
C = xormqr_m('L','T',AOUT,tau,B);
y = linsolve(AOUT,C,opts);
x = zeros(size(y));
x(jpvt(1:length(y)),:) = y;

assert(norm(x-xtrue,'fro')/norm(B,'fro') < cond(A)*tol,...
    'failed forward error test - complex mat')
assert(norm(A*x-B)/norm(B) < tol,...
    'failed backward error test - complex mat')

%% test 2: real least squares

A = real(A);
xtrue = real(xtrue);
B = A*xtrue;

[AOUT,jpvt,tau] = xgeqp3_m(A);
C = xormqr_m('L','T',AOUT,tau,B);
y = linsolve(AOUT,C,opts);
x = zeros(size(y));
x(jpvt(1:length(y)),:) = y;

assert(norm(x-xtrue,'fro')/norm(B,'fro') < cond(A)*tol,...
    'failed forward error test - real mat')
assert(norm(A*x-B)/norm(B) < tol,...
    'failed backward error test - real mat')
