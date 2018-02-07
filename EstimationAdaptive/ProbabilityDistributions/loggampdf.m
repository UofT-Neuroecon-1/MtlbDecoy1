function y = loggampdf(x,a,b,u)
%LOGGAMPDF Log-Gamma probability density function
%   Y = LOGGAMPDF(X,A,B,U) returns the probability density function of the
%   Log-Gamma Distribution with shape parameters A and B and location
%   parameter U. 
%
%   Type: Continuous, semi-bounded
%   Restrictions:
%     X >= U
%     A, B > 0
%     U >= 0
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%

%   Mike Sheppard
%   Last Modified 25-Jun-2011


if nargin < 4
   error('loggampdf:TooFewInputs','Requires four input arguments.');
end

[errorcode, x,a,b,u] = distchck(4,x,a,b,u);

if errorcode > 0
    error('loggampdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize y to zero.
if isa(x,'single') || isa(a,'single') || isa(b,'single') || isa(u,'single')
   y = zeros(size(x),'single');
else
   y = zeros(size(x));
end


k=(a>0 & b>0 & u>=0 & x>=u);
if any(k)
    term1=(-a(k).*log(b(k)));
    term2=(-(1+b(k))./b(k)).*log(1+x(k)-u(k));
    term3=gammaln(a(k));
    y(k)=(log(1+x(k)-u(k)).^(-1+a(k))).*exp(term1+term2-term3);
end

%Else not valid or zero
k1=(a>0 & b>0 & u>=0 & x<u);
y(k1)=0; y(~(k|k1))=NaN;

%Round off
y(y<0)=0;

end