function x = tlsinv(p,mu,sigma,nu)
%TLSINV Inverse of the T Location-Scale cumulative distribution function
%   X = TLSINV(P,MU,SIGMA,NU) returns the inverse of the cumulative 
%   distribution function of the T Location-Scale Distribution with 
%   location parameter MU, scale parameter SIGMA, and NU degrees of
%   freedom, at the values in P.
%
%   The size of X is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.     
%
%   Type: Continuous, unbounded
%   Restrictions:
%     SIGMA, NU > 0
%
%   See also TLSPDF, TLSCDF, TLSSTAT, TLSFIT, TLSLIKE, TLSRND
%

%   Mike Sheppard
%   Last Modified 10-Dec-2011



if nargin < 4
    error('tlsinv:TooFewInputs',...
          'Requires four input arguments.');
end

[errorcode p mu sigma nu] = distchck(4,p,mu,sigma,nu);

if errorcode > 0
    error('tlsinv:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize X to zero.
if isa(p,'single') || isa(mu,'single') || isa(sigma,'single') || isa(nu,'single')
   x = zeros(size(p),'single');
else
   x = zeros(size(p));
end

% Return NaN for out of range parameters.
sigma(sigma <= 0) = NaN;
nu(nu <= 0) = NaN;

%From addtls.m in statistics toolbox
%-----
x = tinv(p,nu).*sigma + mu;
%-----

end