function s = nctsf(x,nu,delta)
%NCTSF Noncentral T survival function
%   S = NCTSF(X,NU,DELTA) returns the survival function of the
%   noncentral T distribution with NU degrees of freedom and
%   noncentrality parameter, DELTA, at the values in X.
%
%   The size of S is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   See also NCTPDF, NCTCDF, NCTINV, NCTSTAT, NCTFIT, NCTLIKE,
%            NCTRND, NCTHAZ
%

%   Mike Sheppard
%   Last Modified 26-Dec-2011

if nargin <  3,
    error(message('nctsf:TooFewInputs'));
end

try
    s = 1 - nctcdf(x,nu,delta);
catch
    error('nctsf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end