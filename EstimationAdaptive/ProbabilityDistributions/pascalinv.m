function x = pascalinv(y,n,p)
%PASCALINV Inverse of the Pascal cumulative distribution function (cdf)
%   X = PASCALINV(Y,N,P) returns the inverse of the Pascal cdf with
%   parameters N and P. Since the Pascal distribution is
%   discrete, PASCALINV returns the least integer X such that 
%   the Pascal cdf evaluated at X, equals or exceeds Y.  
%
%   Note that X takes the values N,N+1,N+2,...
%
%   Type: Discrete, Semi-bounded
%   Restrictions:
%      y>=n
%      n>=0 (integer)
%      0<=p<=1
%
%   The size of X is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs. 
%

%   Mike Sheppard
%   Last Modified 18-Jun-2011


if nargin < 3
    error('pascalinv:TooFewInputs',...
          'Requires three input argument.');
end

[errorcode, y, n,p] = distchck(3,y,n,p);

if errorcode > 0
    error('pascalinv:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

%PASCALCDF(X,N,P)=NBINCDF(X-N,N,P)
%Let nbincdf handle all error checking
x=n+nbininv(y,n,p);

end