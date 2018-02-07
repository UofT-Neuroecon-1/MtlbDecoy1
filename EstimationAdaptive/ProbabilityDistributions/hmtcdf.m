function p = hmtcdf(x,n)
%HMTCDF Heads-Minus-Tails cumulative distribution function
%   P = HMTCDF(X,N) returns the cumulative distribution function of the 
%   Heads-Minus-Tails Distribution of having an absolute difference of 2X,
%   or less, given a fair coin is tossed 2N number of times.
%
%   The size of P is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Distribution: Discrete, bounded, {0,...,N}
%   Restrictions:
%      N >= 1   (integer)
%
%   See also HMTPDF, HMTINV, HMTSTAT, HMTFIT,
%            HMTLIKE, HMTRND, HMTSF, HMTHAZ
%

%   Mike Sheppard
%   Last Modified 24-Dec-2011


if nargin ~= 2
    error('hmtcdf:TooFewInputs',...
          'Requires two input argument.'); 
end


try
    %Expand size if necessary
    x=x+zeros(size(n));
    n=n+zeros(size(x));
catch err
    error('hmtcdf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end



% Initialize Y to zero.
if isa(x,'single') || isa(n,'single')
   p = zeros(size(x),'single');
else
   p = zeros(size(x));
end


% Return NaN if any arguments are outside of their respective limits.
okparam = (n>=1 & n<Inf & n==round(n));
okvar = (0 <= x & x <= n) & (x==round(x));
ok=(okparam & okvar);
p(~okparam)=NaN;
p(okparam & x<0)=0;
p(okparam & x>n)=1;


scalarn=isscalar(n);

if any(ok)
    x=x(ok); n=n(ok); 
    val=max(x(:));
    I=(0:val)';
    if scalarn
        tmp=cumsum(hmtpdf(I,n(1)));
        p(ok)=tmp(x+1);
    else
        n1=sum(ok); n2=numel(I);
        compare=repmat(I,1,n1);
        index=repmat(x,n2,1);
        nbig=repmat(n,n2,1);
        y0=hmtpdf(compare,nbig);
        y0(compare>index)=0;
        p(ok)=sum(y0,1);
    end
end

%Catch round off
p(p<0)=0; p(p>1)=1;

end