function p = normdiffcdf(x,m,v)
%NORMDIFFCDF Normal Difference cumulative distribution function
%   P = NORMDIFFCDF(X,M,V) returns the Normal Difference cumulative
%   distribution function of the difference of two Normal Distributions
%   with means M, and variances V, at the values of X.
%
%   The parameters M and V must be 1x2 vectors.
%
%   Type: Continuous, Unbounded
%   Restrictions:
%      v>0
%
%   The size of Y is the size of X.
%

%   Mike Sheppard
%   Last Modified 19-Jun-2011


if nargin < 3
    error('normdiffcdf:TooFewInputs',...
          'Requires at least three input argument.');
end

%Check sizes
if ~isvector(m) || ~isvector(v) || numel(m)~=2 || numel(v)~=2
    error('normdiffcdf:TooFewInputs',...
          'Mean and variances must be vectors of length 2.');
end

m=m(:)'; v=v(:)';
w=[1 -1]; %Use Normal Sum with weights of +1 and -1
p=normsumcdf(x,m,v,w);

end