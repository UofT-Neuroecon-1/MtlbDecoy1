function s = bnkgbsf(x,a,b)
%BNKGBSF Benktander-Gibrat survival function
%   S = BNKGBSF(X,A,B) returns the survival function of the 
%   Benktander-Gibrat distribution with parameters A and B 
%   at the values in X.
%
%   The size of S is the common sizes of the input arguments. A scalar input
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
%            BNKGBRND, BNKGBHAZ
%

%   Mike Sheppard
%   Last Modified 20-Dec-2011

if nargin < 3
    error('bnkgbsf:TooFewInputs',...
        'Requires at least three input argument.');
end


try
    s=1-bnkgbcdf(x,a,b);
catch
    error('bnkgbsf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end