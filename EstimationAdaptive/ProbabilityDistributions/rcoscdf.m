function p = rcoscdf(x,u,s)
%RCOSCDF Raised Cosine cumulative distribution function
%   P = RCOSCDF(X,U,S) returns the raised cosine cumulative distribution
%   function with mean U and scale parameter S.
%
%   Type: Continuous, Bounded
%   Restrictions:
%   S>0
%   U-S<=X<=U+S
%
%   The size of P is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%

%   Mike Sheppard
%   Last Modified 13-Jun-2011


if nargin < 3
    error('rcoscdf:TooFewInputs',...
          'Requires three input arguments.');
end


[errorcode,x,u,s] = distchck(3,x,u,s);

if errorcode > 0
    error('rcoscdf:InputSizeMismatch',...
          'Non-scalar arguments must match in size.');
end


% Initialize Y to zero.
if isa(x,'single') || isa(u,'single') || isa(s,'single')
    p=zeros(size(x),'single');
else
    p=zeros(size(x));
end

p((x<u-s | x>u+s) & s<0)=NaN;


k=(x>=u-s & x<=u+s & s>0);
if any(k)
    sc_x=(x(k)-u(k))./s(k);
    p(k)=(1+sc_x+(sin(sc_x*pi)./pi))./2;
end

%Bounds
p(x<u-s & s>0)=0;
p(x>u+s & s>0)=1;


end