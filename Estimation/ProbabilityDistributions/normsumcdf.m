function p = normsumcdf(x,m,v,w)
%NORMSUMCDF Normal Sum cumulative distribution function
%   P = NORMSUMCDF(X,M,V,W) returns the Normal Sum cumulative distribution
%   function of the sum of Normal Distributions with means M, variances V,
%   and weight W, at the values of X.
%
%   The parameters M,V,W must be vectors; and the sum is evaluated at all
%   values in X. Default values of W are equal weights of 1.
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
    error('normsumcdf:TooFewInputs',...
          'Requires at least three input argument.');
end

if nargin ==3
   w=1;
end

[errorcode, m,v,w] = distchck(3,m,v,w);

if errorcode > 0
    error('normsumcdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end


% Initialize y to zero.
if isa(m,'single') || isa(v,'single') || isa(w,'single')
    p=zeros(size(m),'single');
else
    p=zeros(size(m));
end

mu=dot(m,w);
vr=dot(v,w.^2);
sig=sqrt(vr);
p=normcdf(x,mu,sig);


end