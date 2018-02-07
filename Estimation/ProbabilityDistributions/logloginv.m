function x = logloginv(p,a,b)
%LOGLOGINV Inverse of the Log-logistic cumulative distribution function
%   X = LOGLOGINV(P,A,B) returns the inverse of the log-logistic cumulative
%   distribution function with scale parameter A and shape parameter B, at
%   the values in P.
%
%   Type: Continuous, semi-bounded
%   Restrictions:
%     A>0
%     B>0
%     0<=P<=1
%
%   The size of X is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%

%   Mike Sheppard
%   Last Modified 23-Jun-2011


if nargin < 3
   error('logloginv:TooFewInputs','Requires three input arguments.');
end

[errorcode, p,a,b] = distchck(3,p,a,b);

if errorcode > 0
    error('logloginv:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize y to zero.
if isa(p,'single') || isa(a,'single') || isa(b,'single')
   x = zeros(size(p),'single');
else
   x = zeros(size(p));
end

k=(a>0 & b>0 & p>=0 & p<=1);

if any(k)
    x(k)=a(k).*((p(k)./(1-p(k))).^(1./b(k)));
end


end