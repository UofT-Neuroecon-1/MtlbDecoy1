function p = asincdf(x,a,b)
%ASINCDF Arcsine cumulative distribution function
%   P = ASINCDF(X,A,B) returns the cumulative distribution function of the
%   Arcsine Distribution on the interval [A,B] at the values in X.
%
%   The size of P is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Default values for A and B are 0 and 1 respectively.
%
%   Distribution: Continuous, bounded, [A,B]
%   Restrictions:
%         A < B
%
%   See also ASINPDF, ASININV, ASINSTAT, ASINFIT,
%            ASINLIKE, ASINRND, ASINSF, ASINHAZ
%

%   Mike Sheppard
%   Last Modified 14-Dec-2011


if (nargin < 1)
    error('asincdf:TooFewInputs',...
        'Requires at least one input argument.');
end

if nargin==1
    a=0; b=1;
elseif nargin==2
    b=1;
end

[errorcode x a b] = distchck(3,x,a,b);

if errorcode > 0
    error('asincdf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

% Initialize P to zero.
if isa(x,'single') || isa(a,'single') || isa(b,'single')
    p = zeros(size(x),'single');
else
    p = zeros(size(x));
end


% Return NaN if any arguments are outside of their respective limits.
okparam = (-Inf<a & a<Inf) & (-Inf<b & b<Inf) & (a<b);
okvar=(a <= x & x <= b);
ok = (okparam & okvar);
p(~okparam)=NaN;
p(okparam & x<a)=0;
p(okparam & x>b)=1;

if any(ok)
    x=x(ok); a=a(ok); b=b(ok);
    p(ok) = (2/pi) * asin( ((x-a) ./ (b-a)).^(1/2) );
end

%Catch round off
p(p<0)=0; p(p>1)=1;

end