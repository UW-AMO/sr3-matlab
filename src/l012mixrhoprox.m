
function z = l012mixrhoprox(x,alpha0,alpha1,alpha2,ifprox)
%%L012MIXRHOPROX utility function for a mixed l0, l1, and l2 penalty
% this function returns either the value of the penalty or 
% the solution of the prox problem. 
%
% if ifprox:
%
% argmin_z alpha0* nnz(z) + alpha1*\|z\|_1 + alpha2*0.5\|z\|_2^2 
%                   + 0.5*\|x-z\|_2^2
%
% where p is one of 0, 1, 2. 
%
% if ~ifprox:
%
% alpha0*nnz(x) + alpha1*\|x\|_1 + 0.5*alpha2*\|x\|_2^2 is
% returned
%
% input:
%   x - vector, as above
%   alpha0, alpha1, alpha2 - weights, as above
%   ifprox - flag, as above
%

if ifprox
    z = sign(x).*(abs(x)-alpha1).*(abs(x)>alpha1)/(1+alpha2);
    fz = (abs(z)~=0)*alpha0+abs(z)*alpha1+0.5*alpha2*abs(z).^2+0.5*abs(z-x).^2;
    z = z.*(fz < 0.5*abs(x).^2);
    %z = wthresh(x,'s',alpha1)/(1+alpha2);
else
    z = alpha0*nnz(x)+alpha1*sum(abs(x))+alpha2*0.5*sum(abs(x).^2);
end
