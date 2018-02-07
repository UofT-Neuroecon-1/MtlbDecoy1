function x = cauchyinv(p,a,b)
%CAUCHYINV Inverse of the Cauchy cumulative distribution function (cdf)
%   X = CAUCHYINV(P,A,B) returns the inverse cumulative distribution
%   function of the Cauchy Distribution with location parameter A and
%   scale parameter B, at the values in P.
%
%   The size of X is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Default values for A and B are 0 and 1, respectively.
%
%   Distribution: Continuous, unbounded, (-Inf,Inf)
%   Restrictions:
%      B > 0
%
%   NOTE: The Cauchy Distribution is also known as the Lorentz Distribution
%
%   See also CAUCHYPDF, CAUCHYCDF, CAUCHYSTAT, CAUCHYFIT,
%            CAUCHYLIKE, CAUCHYRND, CAUCHYSF, CAUCHYHAZ
%

%   Mike Sheppard
%   Last Modified 29-Mar-2011


if nargin<1
    error('cauchyinv:TooFewInputs','Requires at least one input argument.');
end
if nargin<2, a=0; end
if nargin<3, b=1; end

[errorcode p a b] = distchck(3,p,a,b);

if errorcode > 0
    error('cauchyinv:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

% Initialize X to zero.
if isa(p,'single') || isa(a,'single') || isa(b,'single')
    x=zeros(size(p),'single');
else
    x=zeros(size(p));
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (-Inf<a & a<Inf) & (0<b & b<Inf);
okvar = (0 < p) & (p < 1);
ok=(okparam & okvar);
x(~ok)=NaN;
x(okparam & p==0)=-Inf;
x(okparam & p==1)=Inf;

if any(ok)
    p = p(ok); a = a(ok); b = b(ok);
    x(ok) = a + b.*tan(pi.*(p-(1/2)));
end

end