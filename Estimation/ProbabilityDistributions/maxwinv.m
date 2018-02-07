function x = maxwinv(p,a)
%MAXWINV Inverse of the Maxwell-Boltzmann distribution function (cdf)
%   X = MAXWINV(P,A) returns the inverse of the Maxwell-Boltzmann CDF
%   with scale parameter A, at the values in P.
%

%   Mike Sheppard
%   Last Modified 24-Mar-2011

if nargin < 2,
    error('maxwinv:TooFewInputs','Requires two input arguments.');
end

try
    x=sqrt(2*(a.^2).*gammaincinv(p,3/2));
catch
    error('maxwinv:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end
