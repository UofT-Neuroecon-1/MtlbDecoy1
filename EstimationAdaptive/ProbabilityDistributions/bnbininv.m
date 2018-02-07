function x = bnbininv(p,n,a,b)
%BNBININV Inverse of the Beta Negative Binomial cumulative distribution function
%   X = BNBININV(P,N,A,B) returns the inverse cumulative distribution
%   function of the Beta Negative Binomial distribution with
%   parameters A and B, with N successes, at the values in P.
%
%   The size of X is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Since the Beta Negative Binomial distribution is discrete, BNBININV
%   returns the least integer X such that the Beta Negative Binomial cdf
%   evaluated at X, equals or exceeds P.
%
%   Distribution: Discrete, semi-bounded, {0,...,Inf}
%   Restrictions:
%      A , B > 0
%      N >= 0 (integer)
%
%   See also BNBINPDF, BNBINCDF, BNBINSTAT, BNBINFIT,
%            BNBINLIKE, BNBINRND, BNBINSF, BNBINHAZ
%

%   Mike Sheppard
%   Last Modified 16-Dec-2011



if nargin ~= 4
    error('bnbininv:TooFewInputs',...
        'Requires four input arguments.');
end


[errorcode p n a b] = distchck(4,p,n,a,b);

if errorcode > 0
    error('bnbininv:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

%Initialize X to 0.
if isa(p,'single') || isa(n,'single') || isa(a,'single') || isa(b,'single')
    x=zeros(size(p),'single');
else
    x=zeros(size(p));
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (0<a & a<Inf) & (0<b & b<Inf) & (n>=0 & n==round(n));
okvar = (0 < p) & (p < 1);
ok=(okparam & okvar);
x(~ok)=NaN;
x(okparam & p==0)=0;
x(okparam & p==1)=Inf;

if any(ok)
    k=find(ok);
    cumdist=x;
    if isempty(k), return; end
    cumdist(k)=bnbinpdf(0,n(k),a(k),b(k));
    count=0;
    k=k(cumdist(k)<p(k));
    while ~isempty(k)
        x(k)=x(k)+1;
        count=count+1;
        cumdist(k)=cumdist(k)+bnbinpdf(count,n(k),a(k),b(k));
        k=k(cumdist(k)<p(k));
    end
end


end