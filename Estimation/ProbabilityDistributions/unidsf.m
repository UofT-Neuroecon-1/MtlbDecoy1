function s = unidsf(x,n)
%UNIDSF Uniform (discrete) survival function
%   S = UNIDSF(X,N) returns the survival function for a random variable
%   uniform on (1,2,...,N), at the values in X.
%
%   The size of S is the common size of X and N. A scalar input   
%   functions as a constant matrix of the same size as the other input.    
%
%   See also UNIDPDF, UNIDCDF, UNIDINV, UNIDSTAT, UNIDFIT, UNIDLIKE, 
%            UNIDRND, UNIDHAZ
%

%   Mike Sheppard
%   Last Modified 26-Dec-2011


if nargin < 2, 
    error(message('unidsf:TooFewInputs')); 
end

try
    s = 1 - unidcdf(x,n);
catch
    error('unidsf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end
