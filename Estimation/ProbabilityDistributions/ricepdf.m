function y = ricepdf(x,s,sigma)
%RICEPDF Rician probability density function
%   Y = RICEPDF(X,S,SIGMA) returns the probability density function of the
%   Rician Distribution with noncentrality parameter S and scale parameter
%   SIGMA, at the values in X.
%
%   The size of Y is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.     
%
%   Distribution: Continuous, semi-bounded, (0, Inf)
%   Restrictions:
%      S >= 0
%      SIGMA > 0
%
%   See also RICECDF, RICEINV, RICESTAT, RICEFIT,
%            RICELIKE, RICERND, RICESF, RICEHAZ
%

%   Mike Sheppard
%   Last Modified 10-Dec-2011



if nargin ~= 3
    error('ricepdf:TooFewInputs',...
          'Requires three input arguments.');
end

[errorcode x s sigma] = distchck(3,x,s,sigma);

if errorcode > 0
    error('ricepdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize Y to zero.
if isa(x,'single') || isa(s,'single') || isa(sigma,'single')
   y = zeros(size(x),'single');
else
   y = zeros(size(x));
end

% Return NaN for out of range parameters.
s(s < 0) = NaN;
sigma(sigma <= 0) = NaN;

%From addrice.m in statistics toolbox
%-----
x(x<0) = 0;
sigsq = sigma.^2;
rsq = (x.^2 + s.^2)./(2.*sigsq);
z = x./sigsq;
expon = rsq - z.*s;
y = z .* exp(-expon) .* besseli(0,z.*s,1);
y(expon > (log(realmax(class(x)))-1)) = 0; % fix up 0*Inf
%-----


%Round-off
y(y<0)=0;

end