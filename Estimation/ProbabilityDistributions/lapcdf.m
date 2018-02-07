function p = lapcdf(x,u,b)
%LAPCDF Laplace cumulative distribution function
%   P = LAPCDF(X,U,B) returns the Laplace cumulative distribution function
%   with mean U and scale parameter B, at the values in X.
%
%   Type: Continuous, unbounded
%   Restrictions:
%     B>0
%
%   The size of P is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%

%   Mike Sheppard
%   Last Modified 24-Jun-2011


if nargin < 3
   error('lapcdf:TooFewInputs','Requires three input arguments.');
end

[errorcode, x,u,b] = distchck(3,x,u,b);

if errorcode > 0
    error('lapcdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize y to zero.
if isa(x,'single') || isa(u,'single') || isa(b,'single')
   p = zeros(size(x),'single');
else
   p = zeros(size(x));
end

k=(b>0);

if any(k)
    term2=1-exp(-abs(x(k)-u(k))./b(k));
    term1=sign(x(k)-u(k));
    p(k)=0.5.*(1+term1.*term2);
end

% Return NaN for out of range parameters.
p(b<=0) = NaN;

%Round-off
p(p<0)=0; p(p>1)=1;


end