function x = loggaminv(p,a,b,u)
%LOGGAMINV Inverse of the Log-Gamma cumulative distribution function
%   X = LOGGAMINV(P,A,B,U) returns the inverse cumulative distribution 
%   function of the Log-Gamma distribution with shape parameters A and B
%   and location parameter U. 
%
%   Type: Continuous, semi-bounded
%   Restrictions:
%     0<=P<=1
%     A,B>0
%     U>=0
%
%   The size of X is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%

%   Mike Sheppard
%   Last Modified 25-Jun-2011


if nargin < 4
   error('loggaminv:TooFewInputs','Requires four input arguments.');
end

[errorcode, p,a,b,u] = distchck(4,p,a,b,u);

if errorcode > 0
    error('loggaminv:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize X to zero.
if isa(p,'single') || isa(a,'single') || isa(b,'single') || isa(u,'single')
   x = zeros(size(p),'single');
else
   x = zeros(size(p));
end

k=(a>0 & b>0 & u>=0 & p>=0 & p<=1);
if any(k)
    term=gammaincinv(p(k),a(k));
    x(k)=exp(b(k).*term)+u(k)-1;
end

%Edge cases
k0=(a>0 & b>0 & u>=0 & p==0); x(k0)=u(k0);
k1=(a>0 & b>0 & u>=0 & p==1); x(k1)=Inf;


end