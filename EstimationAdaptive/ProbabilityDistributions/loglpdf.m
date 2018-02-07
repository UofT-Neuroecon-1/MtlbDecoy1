function y = loglpdf(x,mu,sigma)
%LOGLPDF Log-Logistic probability density function
%   Y = LOGLPDF(X,MU,SIGMA) returns the probability density function of the
%   Log-Logistic Distribution with log-mean MU and log-scale SIGMA, evaluated
%   at the values in X.
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Default values for MU and SIGMA are 0 and 1, respectively.
%
%   Distribution: Continuous, semi-bounded, (0, Inf)
%   Restrictions:
%      SIGMA>0
%
%   See also LOGLCDF, LOGLINV, LOGLSTAT, LOGLFIT,
%            LOGLLIKE, LOGLRND, LOGLSF, LOGLHAZ
%

%   Mike Sheppard
%   Last Modified 10-Dec-2011


if nargin < 1
    error('loglpdf:TooFewInputs',...
          'Requires at least one input argument.'); 
end
if nargin < 2, mu = 0; end
if nargin < 3, sigma = 1; end


[errorcode, x,mu,sigma] = distchck(3,x,mu,sigma);

if errorcode > 0
    error('loglpdf:InputSizeMismatch',...
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

%Modified from addlogi.m in statistics toolbox
%-----
nonpos = (x < 0);
x(nonpos) = realmin;
z = (log(x) - mu) ./ sigma;
c = ones(size(z));
k = (z>350); % prevent Inf/Inf
if any(k)
    z(k) = -z(k);
    c(k) = -1;
end
y = exp(z.*(1-c.*sigma) - mu) ./ ((1 + exp(z)).^2 .* sigma);
y(nonpos) = 0;
% the first and third of these would happen automatically for x==0, but
% generate LogOfZero warnings.  the second would be NaN.
y(x==0 & sigma<1) = 0;
y(x==0 & sigma==1) = 1;
y(x==0 & sigma>1) = Inf;
%-----

%Round-off
y(y<0)=0;

end