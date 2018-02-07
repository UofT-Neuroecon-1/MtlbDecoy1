function s = gbenfordsf(x,b,n)
%GBENFORDSF Generalized Benford survival function
%   S = GBENFORDSF(X,B,N) returns the survival function of the Generalized
%   Benford distribution with parameters of base B, at the N'th digit,
%   at the values of X.
%    
%   The size of S is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   A set of numbers is said to satisfy Generalized Benford's Law if the 
%   N'th (N>=2) digit X={0,...,B-1} occurs with probability 
%   GBENFORDPDF(X,B,N)
%
%   Type: Discrete, bounded, {0,...,B-1}
%   Restrictions:
%        N >= 2         (N integer)  [Digit Position]
%        B >= 2         (B integer)  [Base]
%
%   See also GBENFORDPDF, GBENFORDCDF, GBENFORDINV, GBENFORDSTAT, 
%            GBENFORDFIT, GBENFORDLIKE, GBENFORDRND, GBENFORDHAZ
%

%   Mike Sheppard
%   Last Modified: 21-Dec-2011


if nargin < 3
    error('gbenfordsf:TooFewInputs',...
          'Requires at least three input argument.'); 
end


try
    s = 1 - gbenfordcdf(x,b,n);
catch
    error('gbenfordsf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end