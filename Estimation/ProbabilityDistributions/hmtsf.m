function s = hmtsf(x,n)
%HMTSF Heads-Minus-Tails survival function
%   S = HMTSF(X,N) returns the survival function of the 
%   Heads-Minus-Tails Distribution of having an absolute difference of 
%   more than 2X, given a fair coin is tossed 2N number of times.
%
%   The size of S is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Type: Discrete, bounded, {0,...,N}
%   Restrictions:
%        N >= 1   (integer)
%
%   See also HMTPDF, HMTCDF, HMTINV, HMTSTAT, HMTFIT, HMTLIKE, 
%            HMTRND, HMTHAZ
%

%   Mike Sheppard
%   Last Modified 3-Jan-2012


if nargin < 2
    error('hmtsf:TooFewInputs',...
          'Requires two input argument.'); 
end


try
    s = 1 - hmtcdf(x,n);
catch
    error('hmtsf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end