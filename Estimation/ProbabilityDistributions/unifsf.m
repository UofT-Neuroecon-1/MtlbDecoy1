function s = unifsf(x,a,b)
%UNIFSF Uniform (continuous) survival function
%   S = UNIFSF(X,A,B) returns the survival function of the uniform
%   distribution on the interval [A,B] at the values in X.
%
%   The size of S is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   By default, A = 0 and B = 1.
%
%   See also UNIFPDF, UNIFCDF, UNIFINV, UNIFSTAT, UNIFFIT, UNIFLIKE, 
%            UNIFRND, UNIFHAZ
%

%   Mike Sheppard
%   Last Modified 26-Dec-2011


if nargin < 1,
    error(message('unifsf:TooFewInputs'));
end

if nargin == 1
    a = 0;
    b = 1;
end

try
    s = 1 - unifcdf(x,a,b);
catch
    error('unifsf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end