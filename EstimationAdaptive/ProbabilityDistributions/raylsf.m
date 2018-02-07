function s = raylsf(x,b)
%RAYLSF  Rayleigh survival function
%   S = RAYLSF(X,B) returns the survival function of the Rayleigh
%   distribution function with parameter B at the values in X.
%
%   The size of S is the common size of X and B. A scalar input
%   functions as a constant matrix of the same size as the other input.
%
%   See also RAYLPDF, RAYLCDF, RAYLINV, RAYLSTAT, RAYLFIT, RAYLLIKE,
%            RAYLRND, RAYLHAZ
%

%   Mike Sheppard
%   Last Modified 26-Dec-2011

if nargin < 1
    error(message('raylsf:TooFewInputs'));
end
if nargin<2
    b = 1;
end

try
    s = 1 - raylcdf(x,b);
catch
    error('raylsf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end
