function x = bnkwblinv(p,a,b)
%BNKWBLINV Inverse of the Benktander-Weibull cumulative distribution function
%   X = BNKWBLINV(P,A,B) returns the inverse cumulative distribution
%   function of the Benktander-Weibull distribution with parameters 
%   A and B at the values in P.
%
%   The size of X is the common sizes of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other input.
%
%   Distribution: Continuous, semi-bounded, [1,Inf)
%   Restrictions:
%      A > 0
%      0 < B <= 1
%
%   Note: The Benktander-Weibull Distribution is also known as the
%   Benktander Distribution of Type II
%
%   See also BNKWBLPDF, BNKWBLCDF, BNKWBLSTAT, BNKWBLFIT,
%            BNKWBLLIKE, BNKWBLRND, BNKWBLSF, BNKWBLHAZ
%

%   Mike Sheppard
%   Last Modified: 16-Dec-2011


if nargin ~= 3
    error('bnkwblinv:TooFewInputs',...
        'Requires three input arguments.');
end

[errorcode p a b] = distchck(3,p,a,b);

if errorcode > 0
    error('bnkwblinv:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

% Initialize x to zero.
if isa(p,'single') || isa(a,'single') || isa(b,'single')
    x=zeros(size(p),'single');
else
    x=zeros(size(p));
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (0<a & a<Inf) & (0<b & b<=1);
okvar = (0 < p & p < 1);
ok=(okparam & okvar);
x(~ok)=NaN;
x(okparam & p==0)=1;
x(okparam & p==1)=Inf;

%Special case, simplifies
k1=(ok & b==1);
if any(k1);
    x(k1)=1-(log(1-p(k1))./a(k1));
end


%More general case, uses Lambert's W function
k=(ok & b~=1);
if any(k)
    pk = p(k); ak = a(k); bk = b(k);
    term=(-ak.*exp(-ak./(-1+bk)).*((1-pk).^(bk./(-1+bk))))./(-1+bk);
    %lambert_term=lambertw(term); %available in symbolic toolbox
    lambert_term=lambertw2(term); %solve below using iteration
    x(k)=(-1*(-1+bk).*lambert_term./ak).^(1./bk);
end

end


function w=lambertw2(x)
%Source:
%http://www.whim.org/nebula/math/lambertw.html
%Rework to close to -1/e part

wnew=x;
%Initialize
wnew(x<10)=0;
wnew(x>10)=log(x(x>10))-log(log(x(x>10)));
crit=sqrt(eps);
w=Inf(size(wnew));
count_limit=100;
count=0;
flag=1;

while( (count<count_limit) && flag )
    %Test flag
    k=isfinite(w) & isfinite(wnew) & isfinite(x);
    if any(k)
        flag=any(abs(w(k)'-wnew(k)')>(crit.*x(k)'));
    end
    count=count+1;
    w=wnew;
    wnew=(x.*exp(-w)+w.^2)./(w+1);
end

end
