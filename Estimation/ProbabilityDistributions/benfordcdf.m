function p = benfordcdf(x,b)
%BENFORDCDF Benford cumulative distribution function
%   P = BENFORDCDF(X,B) returns the cumulative distribution function of the
%   Benford distribution with parameter of base B, at the values of X.
%
%   The size of P is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Default value for B is 10. (Base 10)
%
%   A set of numbers is said to satisfy Benford's Law if the leading digit
%   X={1,...,B-1} occurs with probability BENFORDPDF(X,B)
%
%   Distribution: Discrete, bounded, {1,...,B-1}
%   Restrictions:
%        B >= 2    (B integer)  [Base]
%
%   See also BENFORDPDF, BENFORDINV, BENFORDSTAT, BENFORDFIT,
%            BENFORDLIKE, BENFORDRND, BENFORDSF, BENFORDHAZ
%

%   Mike Sheppard
%   Last Modified 15-Dec-2011


if nargin < 1
    error('benfordcdf:TooFewInputs',...
        'Requires at least one input argument.');
end

if nargin==1
    b=10;
end

try
    x=x+zeros(size(b));
    b=b+zeros(size(x));
catch err
    error('benfordcdf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end



% Return NaN if any arguments are outside of their respective limits.
okparam = (b>=2 & b<Inf & b==round(b));
%okvar = (1 <= x & x < b) & (x==round(x));
%Formula
p=log(1+floor(x))./log(b);
%Catch edge cases
p(~okparam)=NaN;
p(okparam & x<1)=0;
p(okparam & x>=b)=1;
%Catch round off
p(p<0)=0; p(p>1)=1;


end