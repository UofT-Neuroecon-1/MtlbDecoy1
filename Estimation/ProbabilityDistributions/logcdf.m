function p = logcdf(x,mu,sigma)
%LOGCDF Logistic cumulative distribution function
%   P = LOGCDF(X,MU,SIGMA) returns the Logistic cumulative distribution
%   function with mean U and scale parameter SIGMA.
%
%   Type: Continuous, unbounded
%   Restrictions:
%     SIGMA>0
%
%   The size of P is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%

%   Mike Sheppard
%   Last Modified 24-Jun-2011


if nargin < 3
   error('logcdf:TooFewInputs','Requires three input arguments.');
end

[errorcode, x,mu,sigma] = distchck(3,x,mu,sigma);

if errorcode > 0
    error('logcdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

k=(sigma>0);
if any(k)
    p(k) = (1/2)+(1/2).*tanh((x(k)-mu(k))./(2.*sigma(k)));
end

%Else not valid
p(~k)=NaN;

%Round off
p(p<0)=0; p(p>1)=1;

end
