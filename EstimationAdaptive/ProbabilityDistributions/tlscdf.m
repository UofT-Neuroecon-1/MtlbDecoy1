function p = tlscdf(x,mu,sigma,nu)
%TLSCDF T Location-Scale cumulative distribution function
%   P = TLSCDF(X,MU,SIGMA,NU) returns the cumulative distribution function 
%   of the T Location-Scale Distribution with location parameter MU, scale
%   parameter SIGMA, and NU degrees of freedom, at the values in X.
%
%   The size of P is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.     
%
%   Distribution: Continuous, unbounded, (-Inf, Inf)
%   Restrictions:
%      SIGMA, NU > 0
%
%   See also TLSPDF, TLSINV, TLSSTAT, TLSFIT,
%            TLSLIKE, TLSRND, TLSSF, TLSHAZ
%

%   Mike Sheppard
%   Last Modified 10-Dec-2011



if nargin ~= 4
    error('tlscdf:TooFewInputs',...
          'Requires four input arguments.');
end

[errorcode x mu sigma nu] = distchck(4,x,mu,sigma,nu);

if errorcode > 0
    error('tlscdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize P to zero.
if isa(x,'single') || isa(mu,'single') || isa(sigma,'single') || isa(nu,'single')
   p = zeros(size(x),'single');
else
   p = zeros(size(x));
end

% Return NaN for out of range parameters.
sigma(sigma <= 0) = NaN;
nu(nu <= 0) = NaN;

%From addtls.m in statistics toolbox
%-----
p = tcdf((x - mu)./sigma,nu);
%-----


%Catch round off
p(p<0)=0; p(p>1)=1;

end