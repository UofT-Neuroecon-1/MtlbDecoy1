function [m,v] = logstat(mu,sigma)
%LOGLOGSTAT Mean and variance for the Logistic Distribution
%   [M,V] = LOGSTAT(MU, SIGMA) returns the mean and variance for the
%   Logistic Distribution with mean U and scale parameter SIGMA.
%
%   Type: Continuous, unbounded
%   Restrictions:
%     SIGMA>0
%
%   The size of the outputs is the common size of the input arguments. 
%   A scalar input functions as a constant matrix of the same size as
%   the other inputs.
%

%   Mike Sheppard
%   Last Modified 25-Jun-2011


if nargin < 2
    error('logstat:TooFewInputs',...
          'Requires two input argument.');
end


[errorcode, mu,sigma] = distchck(2,mu,sigma);

if errorcode > 0
    error('logstat:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end


% Initialize y to zero.
if isa(mu,'single') || isa(sigma,'single')
    m=NaN(size(mu),'single');
else
    m=NaN(size(mu));
end
v=m;

k=(sigma>0);
if any(k)
    m(k)=mu(k);
    v(k)=(pi^2/3).*(sigma(k)).^2;
end

end
