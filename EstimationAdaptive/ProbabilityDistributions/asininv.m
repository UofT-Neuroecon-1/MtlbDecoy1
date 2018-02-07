function x = asininv(p,a,b)
%ASININV Inverse of the arcsine cumulative distribution function
%   X = ASININV(P,A,B) returns the inverse cumulative distribution function
%   of the Arcsine Distribution on the interval [A,B] at the values in P.
%
%   The size of X is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Default values for A and B are 0 and 1 respectively.
%
%   Distribution: Continuous, bounded, [A,B]
%   Restrictions:
%      A < B
%
%   See also ASINPDF, ASINCDF, ASINSTAT, ASINFIT,
%            ASINLIKE, ASINRND, ASINSF, ASINHAZ
%

%   Mike Sheppard
%   Last Modified 16-May-2012


if (nargin <1)
    error('asininv:TooFewInputs',...
        'Requires at least one input argument.');
end

if nargin ==1
    a=0; b=1;
elseif nargin==2
    b=1;
end

[errorcode p a b] = distchck(3,p,a,b);

if errorcode > 0
    error('asininv:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

% Initialize X to zero.
if isa(p,'single') || isa(a,'single') || isa(b,'single')
    x = zeros(size(p),'single');
else
    x = zeros(size(p));
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (-Inf<a & a<Inf) & (-Inf<b & b<Inf) & (a<b);
okvar = (0<=p & p<=1);
ok=(okparam & okvar);
x(~ok)=NaN;


%Compute transformed inverse
if any(ok)
    p=p(ok); a=a(ok); b=b(ok);
    x(ok) = a +  (b-a).*(sin(pi.*p./2)).^2;
end


end
