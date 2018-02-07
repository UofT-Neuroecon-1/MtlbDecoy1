function h = bnkwblhaz(x,a,b)
%BNKWBLHAZ Benktander-Weibull hazard function
%   H = BNKWBLHAZ(X,A,B) returns the Benktander-Weibull hazard 
%   function with parameters A and B at the values in X.
%
%   The size of H is the common sizes of the input arguments. A scalar input   
%   functions as a constant matrix of the same size as the other input. 
%
%   Type: Continuous, semi-bounded, [1,Inf)
%   Restrictions:
%        A > 0
%        0 < B <= 1
%
%   Note: The Benktander-Weibull Distribution is also known as the
%   Benktander Distribution of Type II
%
%   See also BNKWBLPDF, BNKWBLCDF, BNKWBLINV, BNKWBLSTAT, 
%            BNKWBLFIT, BNKWBLLIKE, BNKWBLRND, BNKWBLSF
%

%   Mike Sheppard
%   Last Modified: 20-Dec-2011

if nargin < 3
    error('bnkwblhaz:TooFewInputs',...
          'Requires at least three input argument.'); 
end

try
    h = bnkwblpdf(x,a,b) ./ bnkwblsf(x,a,b);
catch
    error('bnkwblhaz:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end