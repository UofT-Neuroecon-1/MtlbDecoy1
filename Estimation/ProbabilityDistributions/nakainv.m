function x = nakainv(p,mu,omega)
%NAKAINV Inverse of the Nakagami cumulative distribution function (CDF)
%   X = NAKAINV(P,MU,OMEGA) returns the inverse cumulative distribution 
%   function of the Nakagami distribution with shape parameter MU 
%   and scale parameter OMEGA, at the values in P.
%
%   The size of X is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.     
%
%   Type: Continuous, semi-bounded
%   Restrictions:
%     MU, OMEGA > 0
%
%   See also NAKAPDF, NAKACDF, NAKASTAT, NAKAFIT, NAKALIKE, NAKARND
%

%   Mike Sheppard
%   Last Modified 10-Dec-2011



if nargin < 3
    error('nakainv:TooFewInputs',...
          'Requires three input arguments.');
end

[errorcode p mu omega] = distchck(3,p,mu,omega);

if errorcode > 0
    error('nakainv:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize x to zero.
if isa(p,'single') || isa(mu,'single') || isa(omega,'single')
   x = zeros(size(p),'single');
else
   x = zeros(size(p));
end


% Return NaN for out of range parameters.
mu(mu <= 0) = NaN;
omega(omega <= 0) = NaN;

%From addnaka.m in statistics toolbox
%-----
x(x<0) = 0;
x = sqrt(gaminv(p,mu,omega./mu));
%-----


end