function h = gamhaz(x,a,b)
%GAMHAZ Gamma hazard function.
%   H = GAMHAZ(X,A,B) returns the hazard function of the gamma distribution
%   with shape and scale parameters A and B, respectively, at the values in
%   X.  The size of H is the common size of the input arguments.  A scalar
%   input functions as a constant matrix of the same size as the other
%   inputs.
%
%   Some references refer to the gamma distribution with a single
%   parameter.  This corresponds to the default of B = 1.
%
%   See also GAMPDF, GAMCDF, GAMINV, GAMSTAT, GAMFIT, GAMLIKE,
%            GAMRND, GAMSF
%

%   Mike Sheppard
%   Last Modified 26-Dec-2011

if nargin < 2
    error('gamhaz:TooFewInputs',...
        'Requires at least two input arguments.');
elseif nargin < 3
    b = 1;
end

try
    h = gampdf(x,a,b) ./ gamsf(x,a,b);
catch
    error('gamhaz:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end

