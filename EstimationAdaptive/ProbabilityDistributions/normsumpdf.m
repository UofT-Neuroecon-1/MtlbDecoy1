function y = normsumpdf(x,m,v,w)
%NORMSUMPDF Normal Sum probability density function
%   Y = NORMSUMPDF(X,M,V,W) returns the probability density function of the
%   Normal Sum Distribution of the sum of Normal Distributions with 
%   means M, variances V, and weight W, at the values of X.
%
%   The parameters M and V must be 1x2 vectors. W may be a scalar or 1x2
%   vector. If W is a scalar with value w it acts the same as [w w]
%
%   NORMSUM(X,M,V) is the same as NORMSUM(X,M,V,[1 1])
%
%   Type: Continuous, Unbounded
%   Restrictions:
%      V > 0 (each element)
%
%   The size of Y is the size of the input variable X
%

%   Mike Sheppard
%   Last Modified 4-Jul-2011


if nargin < 3
    error('normsumpdf:TooFewInputs',...
          'Requires at least three input arguments.');
end

if nargin ==3
   w=1;
end

[errorcode, m,v,w] = distchck(3,m,v,w);

if errorcode > 0
    error('normsumpdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Return NaN for out of range parameters.
v(v<0)=NaN;

mu=dot(m,w);
sig=sqrt(dot(v,w.^2)); %Standard Deviation
y=normpdf(x,mu,sig);

%Round off
y(y<0)=0;

end