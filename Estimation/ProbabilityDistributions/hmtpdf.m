function y = hmtpdf(x,n)
%HMTPDF Heads-Minus-Tails probability density function
%   Y = HMTPDF(X,N) returns the probability density function of the
%   Heads-Minus-Tails Distribution of having an absolute difference
%   of heads and tails of 2X given a fair coin is tossed 2N number
%   of times.
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Distribution: Discrete, bounded, {0,...,N}
%   Restrictions:
%      N >= 1   (integer)
%
%   See also HMTCDF, HMTINV, HMTSTAT, HMTFIT, 
%            HMTLIKE, HMTRND, HMTSF, HMTHAZ
%

%   Mike Sheppard
%   Last Modified 24-Dec-2011


if nargin ~= 2
    error('htmpdf:TooFewInputs',...
        'Requires two input arguments.');
end


try
    %Expand size if necessary
    x=x+zeros(size(n));
    n=n+zeros(size(x));
catch err
    error('htmpdf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


% Initialize Y to zero.
if isa(x,'single') || isa(n,'single')
    y = zeros(size(x),'single');
else
    y = zeros(size(x));
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (n>=1 & n<Inf & n==round(n));
okvar = (0 <= x & x <= n) & (x==round(x));
ok=(okparam & okvar);
y(~okparam)=NaN;
y(okparam & ~okvar)=0;

if any(ok)
    x=x(ok); n=n(ok);
    %nchoosek is not vectorized, use gammaln
    term1=log(2)+(-2*n).*log(2);
    term1(x==0)=term1(x==0)-log(2); %exception if x==0
    term2=gammaln(2*n+1)-gammaln(n+x+1)-gammaln(n-x+1);
    y(ok) = exp(term1+term2);
end

%Catch round off
y(y<0)=0;


end