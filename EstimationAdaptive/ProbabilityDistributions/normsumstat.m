function [mu,vr] = normsumstat(m,v,w)
%NORMSUMSTAT Mean and variance of Normal Sum Distribution
%   [MU,VR] = NORMSUMSTAT(M,V,W) returns the mean and variance for the 
%   Normal Sum distribution function of the sum of Normal Distributions
%   with means M, variances V, and weight W.
%
%   The parameters M,V,W must be vectors.
%   Default values of W are equal weights of 1.
%
%   Type: Continuous, Unbounded
%   Restrictions:
%      v>0
%      w>0
%

%   Mike Sheppard
%   Last Modified 19-Jun-2011


if nargin < 2
    error('normsumstat:TooFewInputs',...
          'Requires at least three input argument.');
end

if nargin==2
   w=1;
end

[errorcode, m,v,w] = distchck(3,m,v,w);

if errorcode > 0
    error('normsumstat:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

mu=dot(m,w);
vr=dot(v,w.^2);

end