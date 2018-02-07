function y = nakapdf(x,mu,omega)
%NAKAPDF Nakagami probability density function
%   Y = NAKAPDF(X,MU,OMEGA) returns the probability density function of the
%   Nakagami Distribution with shape parameter MU and scale parameter OMEGA, 
%   at the values in X.
%
%   The size of Y is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.     
%
%   Distribution: Continuous, semi-bounded, (0, Inf)
%   Restrictions:
%      MU, OMEGA > 0
%
%   See also NAKACDF, NAKAINV, NAKASTAT, NAKAFIT,
%            NAKALIKE, NAKARND, NAKASF, NAKAHAZ
%

%   Mike Sheppard
%   Last Modified 10-Dec-2011



if nargin ~= 3
    error('nakapdf:TooFewInputs',...
          'Requires three input arguments.');
end

[errorcode x mu omega] = distchck(3,x,mu,omega);

if errorcode > 0
    error('nakapdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize Y to zero.
if isa(x,'single') || isa(mu,'single') || isa(omega,'single')
   y = zeros(size(x),'single');
else
   y = zeros(size(x));
end

% Return NaN for out of range parameters.
mu(mu <= 0) = NaN;
omega(omega <= 0) = NaN;

%From addnaka.m in statistics toolbox
%-----
x(x<0) = 0;
% equivalent to y = 2.*x .* gampdf(x.^2, mu, omega./mu), but the version here
% puts all the x terms into gampdf, so Inf*0, etc. is handled there.
y = 2.*sqrt(omega./mu).*exp(gammaln(mu+.5) - gammaln(mu)) .* gampdf(x.^2, mu+.5, omega./mu);
%-----

%Round-off
y(y<0)=0;

end