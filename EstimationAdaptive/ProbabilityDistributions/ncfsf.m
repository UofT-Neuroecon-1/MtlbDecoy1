function s = ncfsf(x,nu1,nu2,delta)
%NCFSF Noncentral F survival function
%   S = NCFSF(X,NU1,NU2,DELTA) returns the survival function of the 
%   noncentral F distribution with numerator degrees of freedom (df), 
%   NU1, denominator df, NU2, and noncentrality parameter, DELTA, 
%   at the values in X.
%
%   The size of S is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   See also NCFPDF, NCFCDF, NCFINV, NCFSTAT, NCFFIT, NCFLIKE,
%            NCFRND, NCFHAZ
%

%   Mike Sheppard
%   Last Modified 26-Dec-2011


if nargin <  4,
    error(message('ncfsf:TooFewInputs'));
end

try
    s = 1 - ncfcdf(x,nu1,nu2,delta);
catch
    error('ncfsf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end