function [m,v] = bbinostat(n,a,b)
%BBINOSTAT Mean and variance for the Beta Binomial Distribution
%   [M,V] = BBINOSTAT(N,A,B) returns the mean and variance of the Beta
%   Binomial distribution with parameters A and B, with N trials.
%
%   The sizes of M and V are the common size of the input arguments.  
%   A scalar input functions as a constant matrix of the same size as
%   the other inputs.
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
%   See also BBINOPDF, BBINOCDF, BBINOINV, BBINOFIT, BBINOLIKE, BBINORND,
%            BBINOSF, BBINOHAZ
%

%   Mike Sheppard
%   Last Modified 16-Dec-2011


if nargin < 3
    error('bbinostat:TooFewInputs',...
          'Requires at least three input arguments.'); 
end

[errorcode n a b] = distchck(3,n,a,b);

if errorcode > 0
    error('bbinostat:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (0<a & a<Inf) & (0<b & b<Inf) & (n>=1 & n==round(n));

%Use Beta Distribution first, to catch most errors and NaNs
[m,v]=betastat(a,b);

%Transform to Beta-Binomial
m = m.*n;
v = v.*n.*(a+b+n);

%Correct out of bounds
m(~okparam)=NaN;
v(~okparam)=NaN;

end