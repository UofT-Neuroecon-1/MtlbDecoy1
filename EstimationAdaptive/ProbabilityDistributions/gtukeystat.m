function [m,v]=gtukeystat(lambda1,lambda2,mu,sigma1,sigma2)
%GTUKEYSTAT Mean and variance of the Generalized Tukey Lambda distribution
%   [M,V]=GTUKEYSTAT(LAMBDA1, LAMBDA2, MU, SIGMA1, SIGMA2) returns the 
%   mean and variance of the Generalized Tukey Lambda distribution with location
%   parameter MU, scale parameters SIGMA1 and SIGMA2, and shape parameters
%   LAMBDA1, LAMBDA2
%
%   The size of the output is the common size of the input arguments. 
%   A scalar input functions as a constant matrix of the same size as 
%   the other inputs.
%

%   Mike Sheppard
%   Last Modified 30-May-2011




if nargin < 5
    error('gtukeystat:TooFewInputs',...
          'Requires five input arguments.');
end


[errorcode, lambda1, lambda2, mu, sigma1, sigma2] = distchck(5,lambda1,lambda2,mu,sigma1,sigma2);

if errorcode > 0
    error('gtukeystat:InputSizeMismatch',...
          'Non-scalar arguments must match in size.');
end

%for ease of programming
L1=lambda1; L2=lambda2; u=mu; s1=sigma1; s2=sigma2;

% Initialize x to NaN
if isa(L1,'single') || isa(L2,'single') || isa(u,'single') || isa(s1,'single') || isa(s2,'single')
   m = zeros(size(L1),'single');
else
   m = zeros(size(L1));
end
v=m;


k1=(L1>-1)&(L2>-1);
if any(k1)
    uk=u(k1); L1k=L1(k1); L2k=L2(k1); s1k=s1(k1); s2k=s2(k1);
    m(k1)=uk+(s2k./(1+L2k))-(s1k./(1+L1k));
end

k2=(L1>(-1/2))&(L2>(-1/2));
if any(k2)
    L1k=L1(k2); L2k=L2(k2); s1k=s1(k2); s2k=s2(k2);
    term1=(s1k.^2)./((1+L1k).^2.*(1+2.*L1k));
    term2=(2.*s1k.*s2k)./(L1k.*(1+L1k).*(L2k+(L2k.^2)));
    term3=(s2k.^2)./((L2k.^2).*((1+L2k).^2));
    term4=((s2k.^2).*beta(1,1+2.*L2k))./(L2k.^2);
    term5=2.*s1k.*s2k.*exp(gammaln(L1k)+gammaln(L2k)-gammaln(2+L1k+L2k));
    v(k2)=term1+term2-term3+term4-term5;
end




end