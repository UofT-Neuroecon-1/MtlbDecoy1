function y = rcospdf(x,u,s)
%RCOSPDF Raised Cosine probability density function
%   Y = RCOSPDF(X,U,S) returns the probability density function of the
%   Raised Cosine Distribution with mean U and scale parameter S, at the
%   values in X.
%
%   Type: Continuous, Bounded
%   Restrictions:
%        S > 0
%        U-S <= X <= U+S
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%

%   Mike Sheppard
%   Last Modified 13-Jun-2011


if nargin < 3
    error('rcospdf:TooFewInputs',...
          'Requires three input arguments.');
end


[errorcode,x,u,s] = distchck(3,x,u,s);

if errorcode > 0
    error('rcospdf:InputSizeMismatch',...
          'Non-scalar arguments must match in size.');
end


% Initialize Y to zero.
if isa(x,'single') || isa(u,'single') || isa(s,'single')
    y=zeros(size(x),'single');
else
    y=zeros(size(x));
end

y((x<u-s | x>u+s) & s<0)=NaN;


k=(x>=u-s & x<=u+s & s>0);
if any(k)
    y(k)=(1+cos(pi*(x(k)-u(k))./s(k)))./(2*s(k));
end

%round off
y(y<0)=0;


end