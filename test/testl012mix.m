
%% test 1: correctness of l012mixrhoprox function

ntest = 100;
nx = 101;
nw = 10001; % must be odd to test the point w = 0

fl012mix = @(w,x,alpha0,alpha1,alpha2) alpha0*(abs(w)>eps(1)) ...
    + alpha1*abs(w) + alpha2*0.5*abs(w).^2 + 0.5*abs(w-x).^2;

ifprox = 1;
iffail= 0;

for i = 1:ntest
    % set up random alpha0, alpha1, alpha2 weights
    alpha0 = rand(); alpha1 = rand(); alpha2 = rand();
        
    xs = linspace(-2,2,nx);
    ws = linspace(-2,2,nw);
    
    hw = ws(2)-ws(1);
    
    for i = 1:length(xs)
        x = xs(i);
        
        fs = fl012mix(ws,x,alpha0,alpha1,alpha2);
        [fmin,imin] = min(fs);
        wmin = ws(imin);
        wprox = l012mixrhoprox(x,alpha0,alpha1,alpha2,ifprox);
        fwprox = fl012mix(wprox,x,alpha0,alpha1,alpha2); 
        
        % compare prox with brute force min
        
        if ~(abs(wmin-wprox) < 2*hw)
            iffail = iffail+1;
            fprintf('fail\n');
            fprintf('f(wmin)  = %e\n',fmin);
            fprintf('f(wprox) = %e\n',fwprox);
        end
        
    end
end

assert(iffail==0, ...
    sprintf('failed prox correctness test % i times out of %i',...
        iffail,nx*ntest));