function y = mappdf(x)
%MAPPDF Map-Airy probability density function
%   Y = MAPPDF(X) returns the probability density function of the
%   Map-Airy Distribution at the values in X
%
%   Type: Continuous, unbounded
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Mike Sheppard
%   Last Modified 23-Jun-2011


if nargin < 1
    error('mappdf:TooFewInputs',...
          'Requires one input argument.');
end

%Valid for all real X
y = 2.*exp(-2.*(x.^3)./3).*(x.*airy(x.^2)-airy(1,x.^2));

%Round off
y(y<0)=0;

end