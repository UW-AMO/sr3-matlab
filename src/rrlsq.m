function [x, w] = rrlsq(A,b,varargin)
%RRLSQ Relaxed pursuit method for regularized least squares problems
% of the form:
%   0.5*norm(A*x-b,2)^2 + lam*rho(w) + 0.5*kap*norm(D*x-w,2)^2
% over x and w. The output w represents a regularized solution of 
% the least squares problem described by A and b. 
%
% Required input (positional):
%
%   A   double precision real or complex matrix (dimension, say, MxN)
%   b   double precision real or complex vector (length M)
%
% Parameter input:
%
%   'x0'        initial guess, decision variable (default zeros(N,1))
%   'w0'        initial guess, regularized decision variable (default
%               zeros(N,1))
%   'D'         regularization pre-multiplication matrix as in formula
%               (default eye(N))
%   'lam'       hyper-parameter, control strength of rho (default 1.0)
%   'kap'       hyper-parameter, control strength of the quadratic penalty
%               (default 1.0)
%   'itm'       maximum number of iterations (default 100)
%   'tol'       terminate if change in w (in l2 norm) is less than tol
%               (default 1e-6)
%   'ptf'       print every ptf iterations (don't print if 0). (default 0)
%   'mode'      '2': rho = 0.5*squared 2 norm, i.e. 0.5*sum(abs(x).^2)
%               '1': rho = 1 norm, i.e. sum(abs(x))
%               '0': rho = 0 norm, i.e. nnz(x)
%               'mixed': rho = sum of 0, 1, and squared 2 norms with 
%                weights l0w, l1w, and l2w
%               'other': rho and rhoprox must be provided
%               (default '1')
%   'l0w'       weight of l0 norm for 'mixed' mode (default 0.0)
%   'l1w'       weight of l1 norm for 'mixed' mode (default 1.0)
%   'l2w'       weight of l2 norm for 'mixed' mode (default 0.0)
%   'rho'       function evaluating regularizer rho (default
%   'rhoprox'   proximal function which, for any alpha, evaluates 
%               rhoprox(x,alpha) = argmin_y alpha*rho(y)+0.5*norm(x-y,2)^2
%
% output:
%   x, w the computed minimizers of the objective
%
% Example:
%
%   >> m = 100; n = 2000; k = 10;
%   >> A = randn(m,n);
%   >> y = zeros(n,1); y(randperm(n,k)) = sign(randn(k,1));
%   >> lam = A.'*b;
%   >> [x,w] = rrlsq(A,b,'lam',lam);
%
% See also RRLSQ_PARAMS, LASSO, LINSOLVE

% Copyright 2018 Travis Askham and Peng Zheng
% Available under the terms of the MIT License

l0rho = @(x) nnz(x);
l0rhoprox = @(x,alpha) wthresh(x,'h',sqrt(2*alpha));
l1rho = @(x) sum(abs(x));
l1rhoprox = @(x,alpha) wthresh(x,'s',alpha);
l2rho = @(x) 0.5*sum(abs(x).^2);
l2rhoprox = @(x,alpha) x/(1.0+alpha);

%% parse inputs

[m,n] = size(A);

p = rrlsq_parse_input(A,b,m,n,varargin{:});

x = p.Results.x0;
w = p.Results.w0;
D = p.Results.D;
lam = p.Results.lam;
kap = p.Results.kap;
itm = p.Results.itm;
tol = p.Results.tol;
ptf = p.Results.ptf;
mode = p.Results.mode;
l0w = p.Results.l0w;
l1w = p.Results.l1w;
l2w = p.Results.l2w;
ifusenormal = p.Results.ifusenormal;

[md,~] = size(D);
if md ~= n
    w = zeros(md,1);
end

if strcmp(mode,'0')
    rho = l0rho;
    rhoprox = l0rhoprox;
elseif strcmp(mode,'1')
    rho = l1rho;
    rhoprox = l1rhoprox;
elseif strcmp(mode,'2')
    rho = l2rho;
    rhoprox = l2rhoprox;
elseif strcmp(mode,'mixed')
    
elseif strcmp(mode,'other')
    rho = p.Results.rho;
    rhoprox = p.Results.rhoprox;
else
    error('incorrect value for mode')
end
            
%% pre-process data

rootkap = sqrt(kap);
alpha = lam/kap;
if ifusenormal
   atareg = (A.'*A) + kap*(D.'*D);
   if issparse(atareg)
    [atacholfac,p,s] = chol(atareg,'upper','vector');
   else
    [atacholfac,p] = chol(atareg,'upper');
    s = 1:n;
   end
   if p ~= 0 
       error('error using normal equations');
   end
   atb = A.'*b;
else
    [AOUT,jpvt,tau] = xgeqp3_m([full(A);rootkap*full(D)]); % qr decomposition
    opts.UT = true;
end

%% start iteration

wm  = w;
err = 2.0*tol;
noi = 0;

normb = norm(b,2);

while err >= tol
    % xstep
    if ifusenormal
        u = atb(s) + kap*w(s);
        x(s) = atacholfac\(atacholfac.'\u);
    else
        u = xormqr_m('L','T',AOUT,tau,[b;rootkap*w]); % apply q* from qr 
        x(jpvt(1:n)) = linsolve(AOUT,u,opts); % solve rx = u
    end
    
    % store D*x
    y = D*x; 
    
    % wstep
    w = rhoprox(y,alpha);
    
    % update convergence information
    obj = 0.5*sum((A*x-b).^2) + lam*rho(w) + 0.5*kap*sum((y-w).^2);
    err = sqrt(sum((w - wm).^2))/normb;
    wm  = w;
    
    % print information
    noi = noi + 1;
    if mod(noi, ptf) == 0
        fprintf('iter %4d, obj %1.2e, err %1.2e\n', noi, obj, err);
    end
    if noi >= itm
        break;
    end
end

end

function p = rrlsq_parse_input(A,b,m,n,varargin)
%RRLSQ_PARSE_INPUT parse the input to RRLSQ
% Sets default values and checks types (within reason)
% See also RRLSQ for details

    l1rho = @(x) sum(abs(x));
    l1rhoprox = @(x,alpha) wthresh(x,'s',alpha);

    defaultx0 = zeros(n,1);
    defaultw0 = zeros(n,1);
    defaultD = speye(n);
    defaultlam = 1.0;
    defaultkap = 1.0;
    defaultitm = 100;
    defaulttol = 1e-6;
    defaultptf = 0;
    defaultmode = '1';
    defaultl0w = 0.0;
    defaultl1w = 1.0;
    defaultl2w = 0.0;
    defaultrho = l1rho;
    defaultrhoprox = l1rhoprox;
    defaultifusenormal = 0;
    
    p = inputParser;
    isdouble = @(x) isa(x,'double');
    isdoublep = @(x) isa(x,'double') && x > 0;
    isdoublepp = @(x) isa(x,'double') && x >= 0;
    isdoublem = @(x) isa(x,'double') && length(x)==m;
    isdoublen = @(x) isa(x,'double') && length(x)==n;
    isnumericp = @(x) isnumeric(x) && x > 0;
    isnumericpp = @(x) isnumeric(x) && x >= 0;    
    isfunhandle = @(x) isa(x,'function_handle');
    
    addRequired(p,'A',isdouble);
    addRequired(p,'b',isdoublem);
    addParameter(p,'x0',defaultx0,isdoublen);
    addParameter(p,'w0',defaultw0,isdoublen);
    addParameter(p,'D',defaultD,isdouble);
    addParameter(p,'lam',defaultlam,isdoublep);
    addParameter(p,'kap',defaultkap,isdoublep);
    addParameter(p,'itm',defaultitm,isnumericp);
    addParameter(p,'tol',defaulttol,isdoublep);
    addParameter(p,'ptf',defaultptf,isnumericpp);
    addParameter(p,'mode',defaultmode,@ischar);
    addParameter(p,'l0w',defaultl0w,isdoublepp);
    addParameter(p,'l1w',defaultl1w,isdoublepp);
    addParameter(p,'l2w',defaultl2w,isdoublepp);
    addParameter(p,'rho',defaultrho,isfunhandle);
    addParameter(p,'rhoprox',defaultrhoprox,isfunhandle);
    addParameter(p,'ifusenormal',defaultifusenormal,@isnumeric);

    parse(p,A,b,varargin{:});
end