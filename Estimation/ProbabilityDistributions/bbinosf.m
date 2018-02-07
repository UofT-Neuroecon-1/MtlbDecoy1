function s = bbinosf(x,n,a,b)
%BBINOSF Beta Binomial survival function
%   S = BBINOSF(X,N,A,B) returns the survival function of the 
%   Beta Binomial Distribution with parameters A and B, with N trials,
%   at the values in X.
%
%   The size of S is the common size of the input arguments. A scalar input
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
%            BBINORND, BBINOHAZ
%

%   Mike Sheppard
%   Last Modified 21-Dec-2011


if nargin < 4
    error('bbinosf:TooFewInputs',...
          'Requires at least four input argument.'); 
end


try
    s = 1 - bbinocdf(x,n,a,b);
catch
    error('bbinosf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


end