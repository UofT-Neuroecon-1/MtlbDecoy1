function y = benfordpdf(x,b)
%BENFORDPDF Benford probability density function
%   Y = BENFORDPDF(X,B) returns the probability density function of the
%   Benford Distribution with base B at the values of X.
%
%   The size of Y is the common size of the input arguments. A scalar input
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
%   See also BENFORDCDF, BENFORDINV, BENFORDSTAT, BENFORDFIT, 
%            BENFORDLIKE, BENFORDRND, BENFORDSF, BENFORDHAZ
%

%   Mike Sheppard
%   Last Modified 15-Dec-2011


if nargin < 1
    error('benfordpdf:TooFewInputs',...
          'Requires at least one input argument.'); 
end

if nargin==1
  b=10;
end


try
    %Match dimensions
    x=x+zeros(size(b));
    b=b+zeros(size(x));
catch err
    error('benfordpdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end


okparam = (b>=2 & b<Inf & b==round(b));
okvar = (1 <= x & x < b) & (x==round(x));
y=log(1+(1./x))./log(b);
y(~okparam)=NaN;
y(okparam & ~okvar)=0;


%Catch round off
y(y<0)=0;

end