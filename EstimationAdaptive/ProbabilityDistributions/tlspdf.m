function y = tlspdf(x,mu,sigma,nu)
%TLSPDF T Location-Scale probability density function
%   Y = TLSPDF(X,MU,SIGMA,NU) returns the probability density function 
%   of the T Location-Scale Distribution with location parameter MU, scale
%   parameter SIGMA, and NU degrees of freedom, at the values in X.
%
%   The size of Y is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.     
%
%   Distribution: Continuous, unbounded, (-Inf, Inf)
%   Restrictions:
%      SIGMA, NU > 0
%
%   See also TLSCDF, TLSINV, TLSSTAT, TLSFIT,
%            TLSLIKE, TLSRND, TLSSF, TLSHAZ
%

%   Mike Sheppard
%   Last Modified 10-Dec-2011



if nargin ~= 4
    error('tlspdf:TooFewInputs',...
          'Requires four input arguments.');
end

[errorcode x mu sigma nu] = distchck(4,x,mu,sigma,nu);

if errorcode > 0
    error('tlspdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize Y to zero.
if isa(x,'single') || isa(mu,'single') || isa(sigma,'single') || isa(nu,'single')
   y = zeros(size(x),'single');
else
   y = zeros(size(x));
end

% Return NaN for out of range parameters.
sigma(sigma <= 0) = NaN;
nu(nu <= 0) = NaN;

%From addtls.m in statistics toolbox
%-----
y = tpdf((x - mu)./sigma,nu)./sigma;
%-----

%Round-off
y(y<0)=0;

end