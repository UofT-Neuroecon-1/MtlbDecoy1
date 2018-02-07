function [m,v] = lapstat(u,b)
%LAPSTAT Mean and variance for the Laplace Distribution
%   [M,V] = LAPSTAT(U,B) returns the mean and variance for the Laplace
%   distribution with mean U and scale parameter B.
%
%   Type: Continuous, unbounded
%   Restrictions:
%     B>0
%
%   The size of the outputs is the common size of the input arguments. 
%   A scalar input functions as a constant matrix of the same size as
%   the other inputs.
%

%   Mike Sheppard
%   Last Modified 24-Jun-2011


if nargin < 2
    error('lapstat:TooFewInputs',...
          'Requires two input argument.');
end


[errorcode, u,b] = distchck(2,u,b);

if errorcode > 0
    error('lapstat:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end


% Initialize y to zero.
if isa(u,'single') || isa(b,'single')
    m=NaN(size(u),'single');
else
    m=NaN(size(u));
end
v=m;


k=(b>0);
if any(k)
    m(k)=u(k);
    v(k)=2.*(b(k).^2);
end


end
