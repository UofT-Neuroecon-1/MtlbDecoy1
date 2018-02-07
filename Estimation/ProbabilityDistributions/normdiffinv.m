function x = normdiffinv(p,m,v)
%NORMDIFFINV Inverse of Normal Difference cumulative distribution function
%   X = NORMDIFFINV(P,M,V) returns the inverse of the Normal Difference
%   cumulative distribution function of the difference of two Normal
%   Distributions with means M, and variances V, at the values in X.
%
%   The parameters M and V must be 1x2 vectors.
%
%   Type: Continuous, Unbounded
%   Restrictions:
%      v>0
%
%   The size of X is the size of P.
%

%   Mike Sheppard
%   Last Modified 19-Jun-2011


if nargin < 3
    error('normdiffinv:TooFewInputs',...
          'Requires at least three input argument.');
end

%Check sizes
if ~isvector(m) || ~isvector(v) || numel(m)~=2 || numel(v)~=2
    error('normdiffinv:TooFewInputs',...
          'Mean and variances must be vectors of length 2.');
end

m=m(:)'; v=v(:)';
w=[1 -1]; %Use Normal Sum with weights of +1 and -1
x=normsuminv(p,m,v,w);

end