function x = gasininv(p,alpha,a,b)
%GASININV Inverse of the generalized arcsine cumulative distribution function
%   X = GASININV(P,ALPHA, A,B) returns the inverse cdf of the
%   Arcsine Distribution with shape parameter ALPHA on the interval [A,B]
%   at the values in P.
%
%   The size of X is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Default values for A and B are 0 and 1 respectively.
%
%   Distribution: Continuous, bounded, [A,B]
%   Restrictions:
%      A < B
%      0 < ALPHA < 1
%
%   See also GASINPDF, GASINCDF, GASINSTAT, GASINFIT,
%            GASINLIKE, GASINRND, GASINSF, GASINHAZ
%

%   Mike Sheppard
%   Last Modified 14-Dec-2011


if (nargin<2)
    error('gasininv:TooFewInputs',...
        'Requires at least two input arguments.');
end

if nargin == 2
    a = 0;
    b = 1;
elseif nargin==3
    b=1;
end

[errorcode p alpha a b] = distchck(4,p,alpha,a,b);

if errorcode > 0
    error('gasininv:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


% Return NaN if any arguments are outside of their respective limits.
okparam = (-Inf<a & a<Inf) & (-Inf<b & b<Inf) & (a<b) & (0<alpha & alpha<1);
okvar = (0<p & p<1);
ok=(okparam & okvar);
x(~okparam)=NaN;
x(okparam & p==0)=a(okparam & p==0);
x(okparam & p==1)=b(okparam & p==1);

if any(ok)
    p=p(ok); alpha=alpha(ok); a=a(ok); b=b(ok);
    %Use BETAINV to catch additional errors
    x(ok) = a + (b-a).*betainv(p,1-alpha,alpha);
end

end