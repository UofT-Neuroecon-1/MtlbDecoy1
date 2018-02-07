function h = bnbinhaz(x,n,a,b)
%BNBINHAZ Beta Negative Binomial hazard function
%   H = BNBINHAZ(X,N,A,B) returns the hazard function of the 
%   Beta Negative Binomial Distribution with parameters A and B, 
%   with N successes, at the values in X.
%
%   The size of H is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Type: Discrete, semi-bounded, {0,...,Inf}
%   Restrictions:
%        A , B > 0
%            N >= 0 (integer)
%
%   See also BNBINPDF, BNBINCDF, BNBININV, BNBINSTAT, BNBINFIT, BNBINLIKE,
%            BNBINRND, BNBINSF,
%

%   Mike Sheppard
%   Last Modified 21-Dec-2011


if nargin < 4
    error('bnbinhaz:TooFewInputs',...
        'Requires at least four input argument.');
end


try
    yt = bnbinpdf(x,n,a,b);
    st = bnbinsf(x,n,a,b);
    h = yt ./ (yt+st);  % +yt term in denominator for discrete r.v.
catch
    error('bnbinhaz:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end