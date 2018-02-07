function x = normsuminv(p,m,v,w)
%NORMSUMINV Inverse of Normal Sum cumulative distribution function
%   X = NORMSUMINV(P,M,V,W) returns the inverse of the Normal Sum
%   cumulative distribution function of the sum of Normal Distributions
%   with means M, variances V, and weight W, at the values in P.
%
%   The parameters M,V,W must be vectors
%   Default values of W are equal weights of 1.
%
%   Type: Continuous, Unbounded
%   Restrictions:
%      v>0
%      w>0
%
%   The size of Y is the size of X.
%

%   Mike Sheppard
%   Last Modified 19-Jun-2011


if nargin < 3
    error('normsuminv:TooFewInputs',...
          'Requires at least three input argument.');
end

if nargin ==3
   w=1;
end

[errorcode, m,v,w] = distchck(3,m,v,w);

if errorcode > 0
    error('normsuminv:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end


% Initialize y to zero.
if isa(m,'single') || isa(v,'single') || isa(w,'single')
    x=zeros(size(p),'single');
else
    x=zeros(size(p));
end

mu=dot(m,w);
vr=dot(v,w.^2);
sig=sqrt(vr);
x=norminv(p,mu,sig);



end