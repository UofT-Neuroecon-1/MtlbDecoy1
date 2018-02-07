function y = bbinopdf(x,n,a,b)
%BBINOPDF Beta Binomial probability density function
%   Y = BBINOPDF(X,N,A,B) returns the probability density function of the
%   Beta Binomial Distribution with parameters A and B, with N trials,
%   at the values in X.
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   The Beta Binomial Distribution is defined as the distribution of
%   observing X successes in N binomial trials where the probability of
%   success is the Beta Distribution with parameters A and B.
%
%   Distribution: Discrete, bounded, {0,...,N}
%   Restrictions:
%         A, B > 0
%         N >= 1       (N integer)
%
%   See also BBINOCDF, BBINOINV, BBINOSTAT, BBINOFIT,
%            BBINOLIKE, BBINORND, BBINOSF, BBINOHAZ
%

%   Mike Sheppard
%   Last Modified 13-May-2012


if nargin ~= 4
    error('bbinopdf:TooFewInputs',...
        'Requires four input arguments.');
end

[errorcode x n a b] = distchck(4,x,n,a,b);

if errorcode > 0
    error('bbinopdf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


% Initialize Y to zero.
if isa(x,'single') || isa(n,'single') || isa(a,'single') || isa(b,'single')
    y=zeros(size(x),'single');
else
    y=zeros(size(x));
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (0<a & a<Inf) & (0<b & b<Inf) & (n>=1 & n==round(n));
okvar = (0 <= x) & (x <= n) & (x==round(x));
ok = (okparam & okvar);
y(okparam & ~okvar)=0;
y(~okparam)=NaN;

if any(ok)
    x=x(ok); n=n(ok); a=a(ok); b=b(ok);
    num=beta(x+a,n-x+b);
    den=beta(a,b);
    %nchoosek is not vectorized, use gammaln
    prodtemp=exp(gammaln(n+1)-gammaln(x+1)-gammaln(n-x+1));
    y(ok)=prodtemp.*num./den;
end

%Catch round off
y(y<0)=0;

end
