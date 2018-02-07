function [m,v]=tukeystat(lambda,mu,sigma)
%TUKEYSTAT Mean and variance of the Tukey Lambda distribution
%   [M,V]=TUKEYSTAT(LAMBDA,MU,SIGMA) returns the mean and variance of the
%   Tukey Lambda distribution with shape parameter LAMBDA, location 
%   parameter MU, and scale parameter SIGMA.
%
%   TUKEYSTAT(LAMBDA) uses the default values for MU=0, SIGMA=1
%   TUKEYSTAT(LAMBDA, MU) uses the default value SIGMA=1;
%
%   The size of the output is the common size of the input arguments. 
%   A scalar input functions as a constant matrix of the same size as 
%   the other inputs.
%

%   Mike Sheppard
%   Last Modified 30-May-2011


if nargin < 1
    error('vonmstat:TooFewInputs',...
          'Requires at least one input argument.');
end

if nargin==1
    sigma=1;
end

if nargin==2
    mu=0;
    sigma=1;
end



[errorcode,lambda,mu, sigma] = distchck(3,lambda,mu,sigma);

if errorcode > 0
    error('tukeystat:InputSizeMismatch',...
          'Non-scalar arguments must match in size.');
end

%for ease of programming
L=lambda; u=mu; s=sigma;

% Initialize x to NaN
if isa(L,'single') || isa(u,'single') || isa(s,'single')
   m = zeros(size(L),'single');
else
   m = zeros(size(L));
end
v=m;

k1=(L>-1);
if any(k1)
    m(k1)=u(k1);
end

k2=(L>(-1/2));
if any(k2)
    uk=u(k2); sk=s(k2); Lk=L(k2);
    %Use gammaln for accuracy
    term1=2.*(sk.^2)./(Lk.^2);
    term2=1./(1+2.*Lk);
    term3=exp(2*gammaln(1+Lk)-gammaln(2.*(Lk+1)));
    v(k2)=term1.*(term2-term3);
end


end