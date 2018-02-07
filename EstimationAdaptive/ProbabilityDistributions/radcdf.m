function p = radcdf(x)
%RADCDF Rademacher cumulative distribution function
%   P = RADCDF(X) returns the Rademacher cumulative distribution function
%   at the values in X.
%
%   Type: Discrete, Bounded
%   Restrictions:
%   X real
%
%   The size of p is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%

%   Mike Sheppard
%   Last Modified 13-Jun-2011


if nargin < 1
    error('radcdf:TooFewInputs',...
          'Requires one input argument.');
end

% Initialize Y to zero.
if isa(x,'single')
    p=zeros(size(x),'single');
else
    p=zeros(size(x));
end

p(x>=1) = 1;
p(x>=-1 & x<1) = .5;

end