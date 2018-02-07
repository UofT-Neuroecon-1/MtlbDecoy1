function p = loglogcdf(x,a,b)
%LOGLOGCDF Log-logistic cumulative distribution function
%   P = LOGLOGCDF(X,A,B) returns the Log-logistic cumulative distribution
%   function with scale parameter A and shape parameter B, at the values
%   in X.
%
%   Type: Continuous, semi-bounded
%   Restrictions:
%     A>0
%     B>0
%
%   The size of P is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%

%   Mike Sheppard
%   Last Modified 23-Jun-2011


if nargin < 3
   error('loglogcdf:TooFewInputs','Requires three input arguments.');
end

[errorcode, x,a,b] = distchck(3,x,a,b);

if errorcode > 0
    error('loglogcdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize y to zero.
if isa(x,'single') || isa(a,'single') || isa(b,'single')
   p = zeros(size(x),'single');
else
   p = zeros(size(x));
end

k=(a>0 & b>0 & x>=0);

if any(k)
    xa=x(k)./a(k);
    p(k)=1./(1+xa.^(-b(k)));
end

% Return NaN for out of range parameters.
p(a<=0 | b<=0 | x<0) = NaN;

%Round-off
p(p<0)=0; p(p>1)=1;


end