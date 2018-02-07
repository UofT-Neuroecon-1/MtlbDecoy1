function h = hmthaz(x,n)
%HMTHAZ Heads-Minus-Tails hazard function
%   H = HMTHAZ(X,N) returns the hazard function of the 
%   Heads-Minus-Tails Distribution of having an absolute difference
%   of heads and tails of 2X given a fair coin is tossed 2N number
%   of times.
%
%   The size of H is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Type: Discrete, bounded, {0,...,N}
%   Restrictions:
%        N >= 1   (integer)
%
%   See also HMTPDF, HMTCDF, HMTINV, HMTSTAT, HMTFIT, HMTLIKE, 
%            HMTRND, HMTSF, HMTHAZ
%

%   Mike Sheppard
%   Last Modified 3-Jan-2012


if nargin < 2
    error('hmthaz:TooFewInputs',...
          'Requires two input argument.'); 
end


try
    yt = hmtpdf(x,n);
    st = hmtsf(x,n);
    h = yt ./ (yt+st);  % +yt term in denominator for discrete r.v.
catch
    error('hmthaz:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


end