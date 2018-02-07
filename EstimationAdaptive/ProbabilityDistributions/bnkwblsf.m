function s = bnkwblsf(x,a,b)
%BNKWBLSF Benktander-Weibull survival function
%   S = BNKWBLSF(X,A,B) returns the survival function of the 
%   Benktander-Weibull distribution with parameters A and B
%   at the values in X.
%
%   The size of S is the common sizes of the input arguments. A scalar input   
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
%   See also BNKWBLPDF, BNKWBLCDF, BNKWBLINV, BNKWBLSTAT, BNKWBLFIT, 
%            BNKWBLLIKE, BNKWBLRND, BNKWBLHAZ
%

%   Mike Sheppard
%   Last Modified: 20-Dec-2011

if nargin < 3
    error('bnkwblsf:TooFewInputs',...
          'Requires at least three input argument.'); 
end

try
    s=1-bnkwblcdf(x,a,b);
catch
    error('bnkwblsf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end