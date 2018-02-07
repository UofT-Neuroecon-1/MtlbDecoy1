function y = logpdf(x,mu,sigma)
%LOGPDF Logistic probability density function
%   Y = LOGPDF(X,MU,SIGMA) returns the Logistic probability density
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
   error('logpdf:TooFewInputs','Requires three input arguments.');
end

[errorcode, x,mu,sigma] = distchck(3,x,mu,sigma);

if errorcode > 0
    error('logpdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize y to zero.
if isa(x,'single') || isa(mu,'single') || isa(single,'single')
   y = zeros(size(x),'single');
else
   y = zeros(size(x));
end


k=(sigma>0);
if any(k)
    y(k) = ((sech((x(k)-mu(k))./(2.*sigma(k)))).^2)./(4.*sigma(k));
end

%Else not valid
y(~k)=NaN;

%Round off
y(y<0)=0;

end
