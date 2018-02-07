function [m,v]=wkbstat(a,b,g,d,mu)
%WKBSTAT Mean and variance for the Wakeby distribution
%   [M,V] = WKBSTAT(A,B,G,D,M) returns the mean and variance of the Wakey
%   distribution with shape parameters B and D, scale parameters A
%   and G, and location parameter MU.
%

%   Mike Sheppard
%   Last Modified 27-Mar-2011




if nargin < 5
    error('wkbstat:TooFewInputs',...
          'Requires at least five input arguments.');
end

[errorcode, a, b, g, d, mu] = distchck(6,a,b,g,d,mu);

if errorcode > 0
    error('wkbstat:InputSizeMismatch',...
          'Non-scalar arguments must match in size.');
end

% Initialize x to NaN
if isa(a,'single') || isa(b,'single') || isa(g,'single') || isa(d,'single') || isa(mu,'single')
   m = zeros(size(a),'single');
else
   m = zeros(size(a));
end
v=m;


%Both mean and variance are valid
k=(a>0 & b>0 & g>0 & d<1/2);
if any(k)
    ak=a(k);
    bk=b(k);
    gk=g(k);
    dk=d(k);
    muk=mu(k);
    m(k)=muk+(gk./(1-dk))+(ak./(1+bk));
    v1=(ak^2)./((1+bk).^2*(1+2*bk));
    v2=(2*ak*gk)./((1+bk).*(1+bk-dk).*(-1+dk));
    v3=(gk.^2)/(((-1+dk).^2).*(-1+2*dk));
    v(k)=v1-v2-v3;
end

%Only the mean is valid
k1=(a>0 & b>0 & g>0 & d>=1/2 & d<1);
if any(k1)
    ak=a(k1);
    bk=b(k1);
    gk=g(k1);
    dk=d(k1);
    muk=mu(k1);
    m(k1)=muk+(gk./(1-dk))+(ak./(1+bk));
end


end
