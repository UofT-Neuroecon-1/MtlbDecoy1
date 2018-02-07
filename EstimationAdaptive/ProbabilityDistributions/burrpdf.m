function y = burrpdf(x,c,k)
%BURRPDF Burr probability density function
%   Y = BURRPDF(X,C,K) returns the probability density function of the
%   Burr Distribution with shape parameters C and K, at the values in X.
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Distribution: Continuous, semi-bounded, (0,Inf)
%   Restrictions:
%        C, K > 0
%
%   Note: The Burr Distribution is also known as the Burr Type XII
%   Distribution, and is a special case of the Sing-Maddala Distribution.
%   The notation for the distribution used here is F(X<x) = 1 - (1+x^C)^-K
%
%   See also BURRCDF, BURRINV, BURRSTAT, BURRFIT, 
%            BURRLIKE, BURRRND, BURRSF, BURRHAZ
%

%   Mike Sheppard
%   Last Modified 17-Dec-2011


if nargin ~= 3
    error('burrpdf:TooFewInputs',...
        'Requires three input arguments.');
end


[errorcode x c k] = distchck(3,x,c,k);

if errorcode > 0
    error('burrpdf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

% Initialize Y to zero.
if isa(x,'single') || isa(c,'single') || isa(k,'single')
    y=zeros(size(x),'single');
else
    y=zeros(size(x));
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (0<c & c<Inf) & (0<k & k<Inf);
okvar = (0<x & x<Inf);
ok=(okparam & okvar);
y(okparam & ~okvar)=0;
y(~okparam)=NaN;


if any(ok),
    x = x(ok); c = c(ok); k = k(ok);
    y(ok)=c.*k.*(x.^(c-1))./((1+(x.^c)).^(k+1));
end

%Catch round off
y(y<0)=0;

end