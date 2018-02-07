function y = radpdf(x)
%RADPDF Rademacher probability density function
%   Y = RADPDF(X) returns the probability density function of the 
%   Rademacher Distribution at the values in X.
%
%   The Rademacher Distribution is defined to be 0.50 if X is +1 or -1,
%   and zero otherwise.
%
%   Type: Discrete, Bounded
%   Restrictions:
%        X = +1 or X = -1
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%

%   Mike Sheppard
%   Last Modified 13-Jun-2011


if nargin < 1
    error('radpdf:TooFewInputs',...
          'Requires one input argument.');
end

% Initialize Y to zero.
if isa(x,'single')
    y=zeros(size(x),'single');
else
    y=zeros(size(x));
end

y(x==1 | x==-1) = 1/2;

end