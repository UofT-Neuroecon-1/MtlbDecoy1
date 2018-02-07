function y = logipdf(x,mu,sigma)
%LOGIPDF Logistic probability density function
%   Y = LOGIPDF(X,MU,SIGMA) returns the probability density function of the
%   Logistic Distribution with mean MU and scale parameter SIGMA, evaluated
%   at the values in X.
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Default values for MU and SIGMA are 0 and 1, respectively.
%
%   Distribution: Continuous, unbounded, (-Inf, Inf)
%   Restrictions:
%      SIGMA > 0
%
%   See also LOGICDF, LOGIINV, LOGISTAT, LOGIFIT,
%            LOGILIKE, LOGIRND, LOGISF, LOGIHAZ
%

%   Mike Sheppard
%   Last Modified 10-Dec-2011


if nargin < 1
    error('logipdf:TooFewInputs',...
          'Requires at least one input argument.'); 
end
if nargin < 2, mu = 0; end
if nargin < 3, sigma = 1; end


[errorcode, x,mu,sigma] = distchck(3,x,mu,sigma);

if errorcode > 0
    error('logipdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize y to zero.
if isa(x,'single') || isa(mu,'single') || isa(sigma,'single')
   y = zeros(size(x),'single');
else
   y = zeros(size(x));
end

% Return NaN for out of range parameters.
sigma(sigma <= 0) = NaN;

%From addlogi.m in statistics toolbox
%-----
z = (x - mu) ./ sigma;
k = (z>350); if any(k), z(k) = -z(k); end % prevent Inf/Inf
y = exp(z) ./ ((1 + exp(z)).^2 .* sigma);
%-----

%Round-off
y(y<0)=0;


end