function y = bnbinpdf(x,n,a,b)
%BNBINPDF Beta Negative Binomial probability density function
%   Y = BNBINPDF(X,N,A,B) returns the probability density function of the
%   Beta Negative Binomial Distribution with parameters A and B,
%   with N successes, at the values in X.
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Distribution: Discrete, semi-bounded, {0,...,Inf}
%   Restrictions:
%        A , B > 0
%            N >= 0 (integer)
%
%   See also BNBINCDF, BNBININV, BNBINSTAT, BNBINFIT,
%            BNBINLIKE, BNBINRND, BNBINSF, BNBINHAZ
%

%   Mike Sheppard
%   Last Modified 16-Dec-2011


if nargin ~= 4
    error('bnbinpdf:TooFewInputs',...
        'Requires four input arguments.');
end

[errorcode x n a b] = distchck(4,x,n,a,b);

if errorcode > 0
    error('bnbinpdf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

% Initialize Y to zero.
if isa(x,'single') || isa(n,'single') || isa(a,'single') || isa(b,'single')
    y=zeros(size(x),'single');
else
    y=zeros(size(x));
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (0<a & a<Inf) & (0<b & b<Inf) & (n>=0 & n==round(n));
okvar = (0 <= x) & (x < Inf) & (x==round(x));
ok = (okparam & okvar);
y(okparam & ~okvar)=0;
y(~okparam)=NaN;

if any(ok)
    x=x(ok); n=n(ok); a=a(ok); b=b(ok);
    %Use logarithms to avoid overflows
    num=logpoch(n,x)+logpoch(a,n)+logpoch(b,x);
    den=gammaln(x+1)+logpoch(a+b,n)+logpoch(n+a+b,x);
    y(ok)=exp(num-den);
end

%Catch round off
y(y<0)=0;

end


function lpch=logpoch(a,b)
%Returns the logarithm of the rising Pochhammer symbol x^(n)
%x^(n) = (x+n-1)! / (x-1)! = gamma(x+n) / gamma(x)
%Therefore, log(x^(n))=gammaln(x+n)-gamma(x);
lpch=gammaln(a+b)-gammaln(a);
end
