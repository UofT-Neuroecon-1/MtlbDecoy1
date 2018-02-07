function x = benfordinv(p,b)
%BENFORDINV Inverse of the Benford cumulative distribution function
%   X = BENFORDINV(P,B) returns the inverse cumulative distribution
%   function of the Benford distribution with parameter B at the
%   values in P.
%
%   Since the Benford distribution is discrete, BENFORDINV returns the
%   least integer X, such that the Benford cdf evaluated at X, equals
%   or exceeds P.
%
%   Default value for B is 10. (Base 10)
%
%   Distribution: Discrete, bounded, {1,...,B-1}
%   Restrictions:
%      B >= 2    (B integer)  [Base]
%
%   See also BENFORDPDF, BENFORDCDF, BENFORDSTAT, BENFORDFIT,
%            BENFORDLIKE, BENFORDRND, BENFORDSF, BENFORDHAZ
%

%   Mike Sheppard
%   Last Modified 15-Dec-2011


if nargin < 1
    error('benfordinv:TooFewInputs',...
        'Requires at least one input argument.');
end

if nargin==1, b=10; end


try
    p=p+zeros(size(b));
    b=b+zeros(size(p));
catch err
    error('benfordinv:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

%Initialize X to 0.
if isa(p,'single') || isa(b,'single')
    x=zeros(size(p),'single');
else
    x=zeros(size(p));
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (b>=2 & b<Inf & b==round(b));
okvar = (0 < p & p < 1);
ok=(okparam & okvar);
x(~ok)=NaN;
x(okparam & p==0)=1;
x(okparam & p==1)=b(okparam & p==1)-1;

if any(ok)
    k=find(ok);
    cumdist = x;
    cumdist(k) = benfordpdf(0,b(k));
    count = 0;
    k = k(cumdist(k) < p(k));
    while ~isempty(k)
        x(k) = x(k) + 1;
        count = count + 1;
        cumdist(k) = cumdist(k) + benfordpdf(count,b(k));
        k(cumdist(k) > p(k)) = [];
        k(x(k) >= b(k)) = [];
    end
end

end
