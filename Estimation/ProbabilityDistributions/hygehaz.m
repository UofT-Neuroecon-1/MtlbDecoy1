function h = hygehaz(x,m,k,n)
%HYGEHAZ Hypergeometric hazard function.
%   H = HYGEHAZ(X,M,K,N) returns the hazard function of the hypergeometric
%   distribution with parameters M, K, and N at the values in X.
%
%   The size of H is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   See also HYGEPDF, HYGECDF, HYGEINV, HYGESTAT, HYGEFIT, HYGELIKE, 
%            HYGERND, HYGESF
%

%   Mike Sheppard
%   Last Modified 26-Dec-2011


if nargin < 4,
    error(message('hygehaz:TooFewInputs'));
end

try
    yt = hygepdf(x,m,k,n);
    st = hygesf(x,m,k,n);
    h = yt ./ (yt+st);  % +yt term in denominator for discrete r.v.
catch
    error('hygehaz:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end