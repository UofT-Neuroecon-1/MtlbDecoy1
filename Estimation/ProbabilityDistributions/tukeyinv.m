function x = tukeyinv(p,lambda,mu,sigma)
%TUKEYINV Inverse of the Tukey Lambda distribution
%   X = TUKEYINV(P,LAMBDA, MU, SIGMA) returns the inverse of the Tukey
%   Lambda distribution with shape parameter LAMBDA, location parameter MU,
%   and scale parameter SIGMA.
%
%   TUKEYINV(P, LAMBDA) uses the default values for MU=0, SIGMA=1
%   TUKEYINV(P, LAMBDA, MU) uses the default value SIGMA=1;
%
%   The size of X is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%

%   Mike Sheppard
%   Last Modified 30-May-2011


if nargin < 2
    error('tukeyinv:TooFewInputs',...
          'Requires at least two input arguments.');
end

if nargin==3
    sigma=1;
end

if nargin==2
    mu=0;
    sigma=1;
end

[errorcode, p,lambda,mu, sigma] = distchck(4,p,lambda,mu,sigma);

if errorcode > 0
    error('tukeyinv:InputSizeMismatch',...
          'Non-scalar arguments must match in size.');
end

%for ease of programming
L=lambda; u=mu; s=sigma;

% Initialize x to NaN
if isa(p,'single') || isa(L,'single') || isa(u,'single') || isa(s,'single')
   x = zeros(size(p),'single');
else
   x = zeros(size(p));
end


%Six cases
k1=(L~=0)&(p>0)&(p<1);
k2=(L==0)&(p>0)&(p<1);
k3=(L<=0)&(p==0);
k4=(L>0)&(p==0);
k5=(L<=0)&(p==1);
k6=(L>0)&(p==1);
%otherwise indeterminate

if any(k1)
    uk=u(k1); sk=s(k1); pk=p(k1); Lk=L(k1);
    x(k1)=uk+(sk.*((pk.^Lk)-((1-pk).^Lk))./Lk);
elseif any(k2)
    uk=u(k2); sk=s(k2); pk=p(k2); Lk=L(k2);
    x(k2)=uk+sk.*(log(pk)-log(1-pk));
elseif any(k3)
    x(k3)=-Inf;
elseif any(k4)
    uk=u(k4); sk=s(k4); pk=p(k4); Lk=L(k4);
    x(k4)=uk-(sk./Lk);
elseif any(k5)
    x(k5)=Inf;
elseif any(k6)
    x(k4)=uk+(sk./Lk);
else
    %Do nothing, as indeterminate and already predefined x to be NaN
end


end