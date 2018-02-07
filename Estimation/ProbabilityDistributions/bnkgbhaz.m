function h = bnkgbhaz(x,a,b)
%BNKGBHAZ Benktander-Gibrat hazard function
%   H = BNKGBHAZ(X,A,B) returns the Benktander-Gibrat hazard function
%   with parameters A and B at the values in X.
%
%   The size of H is the common sizes of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other input.
%
%   Type: Continuous, semi-bounded, [1,Inf)
%   Restrictions:
%        A(A+1) >= 2B
%        A , B > 0
%
%   Note: The Benktander-Gibrat Distribution is also known as the
%   Benktander Distribution of Type I
%
%   See also BNKGBPDF, BNKGBCDF, BNKGBINV, BNKGBSTAT, BNKGBFIT, BNKGBLIKE,
%            BNKGBRND, BNKGBSF
%

%   Mike Sheppard
%   Last Modified 20-Dec-2011

if nargin < 3
    error('bnkgbhaz:TooFewInputs',...
        'Requires at least three input argument.');
end


try
    h = bnkgbpdf(x,a,b) ./ bnkgbsf(x,a,b);
catch
    error('bnkgbhaz:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end