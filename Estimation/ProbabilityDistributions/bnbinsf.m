function s = bnbinsf(x,n,a,b)
%BNBINSF Beta Negative Binomial survival function
%   S = BNBINSF(X,N,A,B) returns the survival function of the 
%   Beta Negative Binomial Distribution with parameters A and B, 
%   with N successes, at the values in X.
%
%   The size of S is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Type: Discrete, semi-bounded, {0,...,Inf}
%   Restrictions:
%        A , B > 0
%            N >= 0 (integer)
%
%   See also BNBINPDF, BNBINCDF, BNBININV, BNBINSTAT, BNBINFIT, BNBINLIKE,
%            BNBINRND, BNBINHAZ
%

%   Mike Sheppard
%   Last Modified 21-Dec-2011


if nargin < 4
    error('bnbinsf:TooFewInputs',...
        'Requires at least four input argument.');
end



try
    s = 1 - bnbincdf(x,n,a,b);
catch
    error('bnbinsf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end