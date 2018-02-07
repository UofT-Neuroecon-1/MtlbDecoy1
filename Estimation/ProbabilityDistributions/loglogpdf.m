function y = loglogpdf(x,a,b)
%LOGLOGPDF Log-logistic probability density function
%   Y = LOGLOGPDF(X,A,B) returns the Log-logistic probability density
%   function with scale parameter A and shape parameter B, at the values in
%   X.
%
%   Type: Continuous, semi-bounded
%   Restrictions:
%     A>0
%     B>0
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%

%   Mike Sheppard
%   Last Modified 23-Jun-2011


if nargin < 3
   error('loglogpdf:TooFewInputs','Requires three input arguments.');
end

[errorcode, x,a,b] = distchck(3,x,a,b);

if errorcode > 0
    error('loglogpdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize y to zero.
if isa(x,'single') || isa(a,'single') || isa(b,'single')
   y = zeros(size(x),'single');
else
   y = zeros(size(x));
end

k=(a>0 & b>0 & x>=0);

if any(k)
    xa=x(k)./a(k);
    ba=b(k)./a(k);
    num=ba.*(xa.^(b(k)-1));
    den=(1+(xa.^b(k))).^2;
    y(k)=num./den;
end

% Return NaN for out of range parameters.
y(a<=0 | b<=0 | x<0) = NaN;


end