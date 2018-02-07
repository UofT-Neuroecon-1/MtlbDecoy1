function x = gtukeyinv(p,lambda1,lambda2,mu,sigma1,sigma2)
%GTUKEYINV Inverse of the Generalized Tukey Lambda distribution
%   X = GTUKEYINV(P,LAMBDA1, LAMBDA2, MU, SIGMA1, SIGMA2) returns the 
%   inverse of the Generalized Tukey Lambda distribution with location
%   parameter MU, scale parameters SIGMA1 and SIGMA2, and shape parameters
%   LAMBDA1, LAMBDA2
%
%   The size of X is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%

%   Mike Sheppard
%   Last Modified 30-May-2011


if nargin < 6
    error('gtukeyinv:TooFewInputs',...
          'Requires six input arguments.');
end


[errorcode, p, lambda1, lambda2, mu, sigma1, sigma2] = distchck(6,p,lambda1,lambda2,mu,sigma1,sigma2);

if errorcode > 0
    error('gtukeyinv:InputSizeMismatch',...
          'Non-scalar arguments must match in size.');
end

%for ease of programming
L1=lambda1; L2=lambda2; u=mu; s1=sigma1; s2=sigma2;

% Initialize x to NaN
if isa(p,'single') || isa(L1,'single') || isa(L2,'single') || isa(u,'single') || isa(s1,'single') || isa(s2,'single')
   x = zeros(size(p),'single');
else
   x = zeros(size(p));
end


%Eight cases
k1=(L1~=0)&(L2~=0)&(p>0)&(p<1);
k2=(L1==0)&(L2~=0)&(p>0)&(p<1);
k3=(L1~=0)&(L2==0)&(p>0)&(p<1);
k4=(L1==0)&(L2==0)&(p>0)&(p<1);
k5=(L1<=0)&(p==0);
k6=(L1>0)&(p==0);
k7=(L2<=0)&(p==1);
k8=(L2>0)&(p==1);
%otherwise indeterminate

if any(k1)
    pk=p(k1); L1k=L1(k1); L2k=L2(k1); uk=u(k1); s1k=s1(k1); s2k=s2(k1);
    term1=(-1+(pk.^L1k)).*s1k./L1k;
    term2=(-1+((1-pk).^L2k)).*s2k./L2k;
    x(k1)=term1-term2+uk;
elseif any(k2)
    pk=p(k2); L2k=L2(k2); uk=u(k2); s1k=s1(k2); s2k=s2(k2);
    term1=s1k.*log(pk);
    term2=(-1+((1-pk).^L2k)).*s2k./L2k;
    x(k2)=term1-term2+uk;
elseif any(k3)
    pk=p(k3); L1k=L1(k3); uk=u(k3); s1k=s1(k3); s2k=s2(k3);
    term1=(-1+(pk.^L1k)).*s1k./L1k;
    term2=s2k.*log(1-pk);
    x(k3)=term1-term2+uk;
elseif any(k4)
    pk=p(k4); uk=u(k4); s1k=s1(k4); s2k=s2(k4);
    term1=s1k.*log(pk);
    term2=s2k.*log(1-pk);
    x(k4)=term1-term2+uk;
elseif any(k5)
    x(k5)=-Inf;
elseif any(k6)
    L1k=L1(k6); uk=u(k6); s1k=s1(k6);
    x(k6)=uk-(s1k./L1k);
elseif any(k7)
    x(k7)=Inf;
elseif any(k8)
    L2k=L2(k1); uk=u(k1); s2k=s2(k1);
    x(k8)=uk+(s2k./L2k);
else
    %Do nothing, as indeterminate and already predefined x to be NaN
end



end