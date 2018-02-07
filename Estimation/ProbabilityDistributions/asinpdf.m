function y = asinpdf(x,a,b)
%ASINPDF Arcsine probability density function
%   Y = ASINPDF(X,A,B) returns the probability density function of the
%   Arcsine Distribution on the interval [A,B] at the values in X.
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Default values for A and B are 0 and 1 respectively.
%
%   Distribution: Continuous, bounded, [A,B]
%   Restrictions:
%         A < B
%
%   See also ASINCDF, ASININV, ASINSTAT, ASINFIT, 
%            ASINLIKE, ASINRND, ASINSF, ASINHAZ
%

%   Mike Sheppard
%   Last Modified 13-May-2012


if (nargin ~=1)&&(nargin ~=3)
    error('asinpdf:TooFewInputs',...
        'Requires either one or three input arguments.');
end

if nargin == 1
    a = 0;
    b = 1;
end

[errorcode x a b] = distchck(3,x,a,b);

if errorcode > 0
    error('asinpdf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

% Initialize Y to zero.
if isa(x,'single') || isa(a,'single') || isa(b,'single')
    y = zeros(size(x),'single');
else
    y = zeros(size(x));
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (-Inf<a & a<Inf) & (-Inf<b & b<Inf) & (a<b);
okvar=(a <= x) & (x <= b);
ok = (okparam & okvar);
y(okparam & ~okvar)=0;
y(~okparam)=NaN;

if any(ok)
    x=x(ok); a=a(ok); b=b(ok);
    y(ok) = 1 ./ (pi*((x-a).*(b-x)).^(1/2));
end

%Catch round off
y(y<0)=0;

end