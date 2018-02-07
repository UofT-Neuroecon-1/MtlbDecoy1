function [m,v] = loglogstat(a,b)
%LOGLOGSTAT Mean and variance for the Log-logistic Distribution
%   [M,V] = LOGLOGSTAT(A,B) returns the mean and variance for the
%   Log-logistic Distribution with scale parameter A and shape parameter B
%
%   Type: Continuous, semi-bounded
%   Restrictions:
%     A>0
%     B>0
%
%   The size of the outputs is the common size of the input arguments. 
%   A scalar input functions as a constant matrix of the same size as
%   the other inputs.
%

%   Mike Sheppard
%   Last Modified 23-Jun-2011


if nargin < 2
    error('loglogstat:TooFewInputs',...
          'Requires two input argument.');
end


[errorcode, a,b] = distchck(2,a,b);

if errorcode > 0
    error('loglogstat:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end


% Initialize y to zero.
if isa(a,'single') || isa(b,'single')
    m=NaN(size(a),'single');
else
    m=NaN(size(a));
end
v=m;

c=pi./b;

k1=(b>1);
if any(k1)
    m(k1)=a(k1).*c(k1).*csc(c(k1));  %Defined for B>1
end

k2=(b>2);
if any(k2)
    term1=2*c(k2).*csc(2.*c(k2));
    term2=(c(k2).*csc(c(k2))).^2;
    v(k2)=(a(k2).^2).*(term1-term2);  %Defined for B>2
end

end
