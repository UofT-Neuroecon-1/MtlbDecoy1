function s = hygesf(x,m,k,n)
%HYGESF Hypergeometric survival function.
%   S = HYGESF(X,M,K,N) returns the survival function of the hypergeometric
%   distribution with parameters M, K, and N at the values in X.
%
%   The size of S is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   See also HYGEPDF, HYGECDF, HYGEINV, HYGESTAT, HYGEFIT, HYGELIKE, 
%            HYGERND, HYGEHAZ
%

%   Mike Sheppard
%   Last Modified 26-Dec-2011


if nargin < 4,
    error(message('hygesf:TooFewInputs'));
end

try
    s = 1 - hygecdf(x,m,k,n);
catch
    error('hygesf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end
