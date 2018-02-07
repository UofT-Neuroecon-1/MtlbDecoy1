function [mu,vr] = normdiffstat(m,v)
%NORMDIFFSTAT Mean and variance of Normal Difference Distribution
%   [MU,VR] = NORMDIFFSTAT(M,V) returns the mean and variance for the
%   Normal Difference distribution function of the difference of two Normal
%   Distributions with means M, and variances V
%
%   The parameters M and V must be 1x2 vectors.
%
%   Type: Continuous, Unbounded
%   Restrictions:
%      v>0
%

%   Mike Sheppard
%   Last Modified 19-Jun-2011

if nargin < 2
    error('normdiffstat:TooFewInputs',...
          'Requires two input argument.');
end

%Check sizes
if ~isvector(m) || ~isvector(v) || numel(m)~=2 || numel(v)~=2
    error('normdiffstat:TooFewInputs',...
          'Mean and variances must be vectors of length 2.');
end

m=m(:)'; v=v(:)';
w=[1 -1]; %Use Normal Sum with weights of +1 and -1
[mu,vr]=normsumstat(m,v,w);

end
