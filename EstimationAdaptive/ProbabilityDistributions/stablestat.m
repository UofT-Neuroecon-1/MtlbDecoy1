function [m,v]=stablestat(type,a,b,mu,sigma)
%STABLESTAT Mean and variance for the Stable Distribution
%   [M,V]=STABLESTAT(TYPE,A,B,MU,SIGMA) returns the mean and variance of
%   the Stable distribution of type TYPE with index of stability A,
%   skewness parameter B, location parameter MU, and scale parameter SIGMA.
%
%   Mean is only defined if 1<A<=2
%   Variance is only defined if A=2

%   Mike Sheppard
%   Last Modified 5-Jun-2011


if nargin < 5
    error('stablestat:TooFewInputs',...
          'Requires five input arguments.');
end

[errorcode,a,b,mu,sigma] = distchck(5,type,a,b,mu,sigma);

if errorcode > 0
    error('stablestat:InputSizeMismatch',...
          'Non-scalar arguments must match in size.');
end

% Initialize x to NaN
if isa(type,'single') || isa(a,'single') || isa(b,'single') || isa(mu,'single') || isa(sigma,'single')
   m = zeros(size(a),'single');
else
   m = zeros(size(a));
end
v=m;

k0=((type==0)&(a>1)&(a<=2));
k1=((type==1)&(a>1)&(a<=2));
k2=(a==2);

if any(k0)
    m(k0)=mu(k0)-b(k0).*sigma(k0).*tan(alpha(k0)*pi/2);
end
if any(k1)
    m(k1)=mu(k1);
end
if any(k2)
    v(k2)=2.*(sigma(k2)).^2;
end


end