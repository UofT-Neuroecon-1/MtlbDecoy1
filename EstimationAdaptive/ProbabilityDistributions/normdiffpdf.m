function y = normdiffpdf(x,m,v)
%NORMDIFFPDF Normal Difference probability density function
%   Y = NORMDIFFPDF(X,M,V) returns the probability density function of the
%   Normal Difference Distribution of the difference of two 
%   Normal Distributions with means M, and variances V, at the values of X.
%
%   The parameters M and V must be 1x2 vectors.
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
    error('normdiffpdf:TooFewInputs',...
          'Requires at least three input argument.');
end

%Check sizes
if ~isvector(m) || ~isvector(v) || numel(m)~=2 || numel(v)~=2
    error('normdiffpdf:TooFewInputs',...
          'Mean and variances must be vectors of length 2.');
end

% Return NaN for out of range parameters.
v(v<0)=NaN;

m=m(:)'; v=v(:)';

%Use Normal Sum with weights of +1 and -1
w=[1 -1]; 
y=normsumpdf(x,m,v,w);

%Round off
y(y<0)=0;

end