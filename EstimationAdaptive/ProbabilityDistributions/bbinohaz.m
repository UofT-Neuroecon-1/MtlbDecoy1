function h = bbinohaz(x,n,a,b)
%BBINOHAZ Beta Binomial hazard function
%   H = BBINOHAZ(X,N,A,B) returns the hazard function of the 
%   Beta Binomial Distribution with parameters A and B, with N trials,
%   at the values in X.
%
%   The size of H is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   The Beta Binomial Distribution is defined as the distribution of
%   observing X successes in N binomial trials where the probability of
%   success is the Beta Distribution with parameters A and B.
%
%   Type: Discrete, bounded, {0,...,N}
%   Restrictions:
%         A, B > 0
%         N >= 1       (N integer)
%
%   See also BBINOPDF, BBINOCDF, BBINOINV, BBINOSTAT, BBINOFIT, BBINOLIKE,
%            BBINORND, BBINOSF
%

%   Mike Sheppard
%   Last Modified 21-Dec-2011


if nargin < 4
    error('bbinohaz:TooFewInputs',...
          'Requires at least four input argument.'); 
end


try
    yt = bbinopdf(x,n,a,b);
    st = bbinosf(x,n,a,b);
    h = yt ./ (yt+st);  % +yt term in denominator for discrete r.v.
catch
    error('bbinohaz:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


end