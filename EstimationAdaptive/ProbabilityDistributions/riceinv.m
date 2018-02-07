function x = riceinv(p,s,sigma)
%RICEINV Inverse of the Rician cumulative distribution function (cdf)
%   X = RICEINV(P,S,SIGMA) returns the inverse cumulative distribution 
%   function of the Rician distribution with noncentrality parameter S 
%   and scale parameter SIGMA, at the values in P.
%
%   The size of P is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.     
%
%   Type: Continuous, semi-bounded
%   Restrictions:
%     S >= 0
%     SIGMA > 0
%
%   See also RICEPDF, RICECDF, RICESTAT, RICEFIT, RICELIKE, RICERND
%

%   Mike Sheppard
%   Last Modified 10-Dec-2011


if nargin < 3
    error('riceinv:TooFewInputs',...
          'Requires three input arguments.');
end

[errorcode p s sigma] = distchck(3,p,s,sigma);

if errorcode > 0
    error('riceinv:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize X to zero.
if isa(p,'single') || isa(s,'single') || isa(sigma,'single')
   x = zeros(size(p),'single');
else
   x = zeros(size(p));
end

% Return NaN for out of range parameters.
s(s < 0) = NaN;
sigma(sigma <= 0) = NaN;

%From addrice.m in statistics toolbox
%-----
x = sigma .* sqrt(ncx2inv(p, 2, (s./sigma).^2));
%-----

end