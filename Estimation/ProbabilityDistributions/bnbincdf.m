function p = bnbincdf(x,n,a,b)
%BNBINCDF Beta Negative Binomial cumulative distribution function
%   P = BNBINCDF(X,N,A,B) returns the cumulative distribution function
%   of the Beta Negative Binomial Distribution with parameters A and B, 
%   with N successes, at the values in X.
%
%   The size of P is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Distribution: Discrete, semi-bounded, {0,...,Inf}
%   Restrictions:
%        A , B > 0
%            N >= 0 (integer)
%
%   See also BNBINPDF, BNBININV, BNBINSTAT, BNBINFIT, 
%            BNBINLIKE, BNBINRND, BNBINSF, BNBINHAZ
%

%   Mike Sheppard
%   Last Modified 16-Dec-2011


if nargin ~= 4
    error('bnbincdf:TooFewInputs',...
        'Requires four input argument.');
end


[errorcode x n a b] = distchck(4,x,n,a,b);

if errorcode > 0
    error('bnbincdf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

% Initialize Y to zero.
if isa(x,'single') || isa(n,'single') || isa(a,'single') || isa(b,'single')
    p = zeros(size(x),'single');
else
    p = zeros(size(x));
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (0<a & a<Inf) & (0<b & b<Inf) & (n>=0 & n==round(n));
x=floor(x); %Floor for CDF for discrete
okvar = (0 <= x) & (x < Inf) & (x==round(x));
ok=(okparam & okvar);
p(~okparam)=NaN;
p(okparam & x<0)=0;
p(okparam & x==Inf)=1;


scalarnab=isscalar(b);
%In theory there is a closed form expression for the CDF, however it
%involves the generalized hypergeometric function. This method sums the PDF
%to the required amount.
if any(ok)
    x=x(ok); n=n(ok); a=a(ok); b=b(ok);
    val=max(x(:));
    I=(0:val)';
    if scalarnab
        tmp=cumsum(bnbinpdf(I,n(1),a(1),b(1)));
        p(ok)=tmp(floor(x+1));  %CDF is right continuous
    else
        n1=sum(ok); n2=numel(I);
        compare=repmat(I,1,n1);
        index=repmat(x,n2,1);
        nbig=repmat(n,n2,1);
        abig=repmat(a,n2,1);
        bbig=repmat(b,n2,1);
        y0=bnbinpdf(compare,nbig,abig,bbig);
        y0(compare>index)=0;
        p(ok)=sum(y0,1);
    end
end

%Catch round off
p(p<0)=0; p(p>1)=1;

end