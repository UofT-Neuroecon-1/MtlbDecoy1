function h = benfordhaz(x,b)
%BENFORDHAZ Benford hazard function
%   H = BENFORDHAZ(X,B) returns the Benford hazard function
%   with parameter of base B, at the values of X.
%
%   The size of H is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Default value for B is 10. (Base 10)
%
%   A set of numbers is said to satisfy Benford's Law if the leading digit
%   X={1,...,B-1} occurs with probability BENFORDPDF(X,B)
%
%   Type: Discrete, bounded, {1,...,B-1}
%   Restrictions:
%        B >= 2    (B integer)  [Base]
%
%   See also BENFORDPDF, BENFORDCDF, BENFORDINV, BENFORDSTAT, BENFORDFIT,
%            BENFORDLIKE, BENFORDRND, BENFORDSF
%

%   Mike Sheppard
%   Last Modified 21-Dec-2011


if nargin < 1
    error('benfordhaz:TooFewInputs',...
        'Requires at least one input argument.');
end

if nargin==1
    b=10;
end


try
    yt = benfordpdf(x,b);
    st = benfordsf(x,b);
    h = yt ./ (yt+st);  % +yt term in denominator for discrete r.v.
catch
    error('benfordhaz:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


end