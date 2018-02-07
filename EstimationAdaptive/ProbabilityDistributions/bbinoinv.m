function x = bbinoinv(p,n,a,b)
%BBINOINV Inverse of the Beta Binomial cumulative distribution function
%   X = BBINOINV(P,N,A,B) returns the inverse cumulative distribution
%   function of the Beta Binomial distribution with parameters A and B,
%   at the values in P.
%
%   The size of X is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Since the Beta Binomial distribution is discrete, BBINOINV
%   returns the least integer X such that the Beta-Binomial cdf
%   evaluated at X, equals or exceeds P.
%
%   Distribution: Discrete, bounded, {0,...,N}
%   Restrictions:
%      A, B > 0
%      N >= 1       (N integer)
%
%   See also BBINOPDF, BBINOCDF, BBINOSTAT, BBINOFIT,
%            BBINOLIKE, BBINORND, BBINOSF, BBINOHAZ
%

%   Mike Sheppard
%   Last Modified 23-Mar-2011


if nargin ~= 4
    error('bbinoinv:TooFewInputs',...
        'Requires four input arguments.');
end


[errorcode p n a b] = distchck(4,p,n,a,b);

if errorcode > 0
    error('bbinoinv:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

%Initialize X to 0.
if isa(p,'single') || isa(n,'single') || isa(a,'single') || isa(b,'single')
    x=zeros(size(p),'single');
else
    x=zeros(size(p));
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (0<a & a<Inf) & (0<b & b<Inf) & (n>=1 & n==round(n));
okvar = (0 < p) & (p < 1);
ok=(okparam & okvar);
x(~ok)=NaN;
x(okparam & p==0)=0;
x(okparam & p==1)=n(okparam & p==1);

if any(ok)
    k=find(ok);
    cumdist=x;
    if isempty(k), return; end
    cumdist(k)=bbinopdf(0,n(k),a(k),b(k));
    count=0;
    k=k(cumdist(k)<p(k));
    while ~isempty(k)
        x(k)=x(k)+1;
        count=count+1;
        cumdist(k)=cumdist(k)+bbinopdf(count,n(k),a(k),b(k));
        k=k(cumdist(k)<p(k));
    end
end

end