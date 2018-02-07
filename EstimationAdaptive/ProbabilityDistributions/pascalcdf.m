function y=pascalcdf(x,n,p)
%PASCALCDF Pascal's  cumulative distribution function
%   Y = PASCALCDF(X,N,P) returns the Pascal cumulative distribution
%   function of the number of trials with success probability
%   P before N successes occur, at the values of X.
%
%   Type: Discrete, Semi-bounded
%   Restrictions:
%      x>=n
%      n>=0 (integer)
%      0<=p<=1
%
%   The size of Y is the size of the input variable X.
%

%   Mike Sheppard
%   Last Modified 18-Jun-2011


if nargin < 3
    error('pascalcdf:TooFewInputs',...
          'Requires three input argument.');
end


[errorcode, x, n,p] = distchck(3,x,n,p);

if errorcode > 0
    error('pascalcdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

%PASCALCDF(X,N,P)=NBINCDF(X-N,N,P)
%Let nbincdf handle all error checking
y=nbincdf(x-n,n,p);

end