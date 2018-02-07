function h = exphaz(x,mu)
%EXPHAZ Exponential hazard function.
%   H = EXPHAZ(X,MU) returns the hazard function of the exponential
%   distribution with mean parameter MU, evaluated at the values in X.
%
%   The size of H is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   The default value for MU is 1.
%
%   See also EXPPDF, EXPCDF, EXPINV, EXPSTAT, EXPFIT, EXPLIKE,
%            EXPRND, EXPSF
%

%   Mike Sheppard
%   Last Modified 26-Dec-2011


if nargin < 1
    error('exphaz:TooFewInputs',...
        'Requires at least one input argument.');
end
if nargin < 2
    mu = 1;
end

try
    mu = mu + zeros(size(x));
    h = 1 ./ mu;
catch
   error('exphaz:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end
