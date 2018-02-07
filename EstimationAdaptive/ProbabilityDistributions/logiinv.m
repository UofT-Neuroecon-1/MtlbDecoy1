function x = logiinv(p,mu,sigma)
%LOGIINV Inverse of the logistic cumulative distribution function
%   X = LOGIINV(P,MU,SIGMA) returns the inverse cumulative distribution 
%   function of the Logistic distribution with mean MU and scale 
%   parameter SIGMA, evaluated at the values in P.
%
%   The size of X is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Default values for MU and SIGMA are 0 and 1, respectively.
%
%   Type: Continuous, unbounded
%   Restrictions:
%     SIGMA>0
%
%   See also LOGIPDF, LOGICDF, LOGISTAT, LOGIFIT, LOGILIKE, LOGIRND
%

%   Mike Sheppard
%   Last Modified 10-Dec-2011


if nargin < 1
    error('logiinv:TooFewInputs',...
          'Requires at least one input argument.'); 
end
if nargin < 2, mu = 0; end
if nargin < 3, sigma = 1; end


[errorcode, p,mu,sigma] = distchck(3,p,mu,sigma);

if errorcode > 0
    error('logiinv:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize y to zero.
if isa(p,'single') || isa(mu,'single') || isa(sigma,'single')
   x = zeros(size(p),'single');
else
   x = zeros(size(p));
end

% Return NaN for out of range parameters.
sigma(sigma <= 0) = NaN;


%From addlogi.m in statistics toolbox
%-----
x = logit(p).*sigma + mu;

function logitp = logit(p)
%LOGIT Logistic transformation, handling edge and out of range.
logitp = zeros(size(p));
logitp(p==0) = -Inf;
logitp(p==1) = Inf;
ok = (0<p & p<1);
logitp(ok) = log(p(ok)./(1-p(ok)));
end

%-----



end