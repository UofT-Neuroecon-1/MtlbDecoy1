function x = lapinv(p,u,b)
%LAPINV Inverse of the Laplace cumulative distribution function
%   X = LAPINV(P,U,B) returns the inverse of the Laplace cumulative
%   distribution function with mean U and scale parameter B, at the values
%   in X.
%
%   Type: Continuous, unbounded
%   Restrictions:
%     0<=P<=1
%     B>0
%
%   The size of X is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%

%   Mike Sheppard
%   Last Modified 24-Jun-2011


if nargin < 3
   error('lapinv:TooFewInputs','Requires three input arguments.');
end

[errorcode, p,u,b] = distchck(3,p,u,b);

if errorcode > 0
    error('lapinv:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize y to zero.
if isa(p,'single') || isa(u,'single') || isa(b,'single')
   x = zeros(size(p),'single');
else
   x = zeros(size(p));
end

k=(b>0 & p>0 & p<1);
if any(k)
    term2=log(1-2.*abs(p(k)-0.5));
    term1=b(k).*sign(p(k)-0.5);
    x(k)=u(k)-term1.*term2;
end

%Edge cases
x(b>0 & p==0)=-Inf;
x(b>0 & p==1)=Inf;


end