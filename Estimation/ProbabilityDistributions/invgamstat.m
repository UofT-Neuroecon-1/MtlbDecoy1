function [m,v] = invgamstat(a,b)
%INVGAMSTAT Mean and variance for the Inverse-Gamma Distribution
%   [M,V] = INVGAMSTAT(A,B) returns the mean and variance for the
%   Inverse-Gamma distribution with shape parameter A and scale
%   parameter B.
%
%   The sizes of M and V are the common size of the input arguments.
%   A scalar input functions as a constant matrix of the same size as
%   the other inputs.
%
%   Type: Continuous, semi-bounded
%   Restrictions:
%     A, B > 0
%
%   See also INVGAMPDF, INVGAMCDF, INVGAMINV, INVGAMFIT, 
%            INVGAMLIKE, INVGAMRND, INVGAMSF, INVGAMHAZ
%

%   Mike Sheppard
%   Last Modified 25-Jun-2011


if nargin < 2
    error('invgamstat:TooFewInputs','Requires two input arguments.');
end

% Return NaN for out of range parameters.
b(b<0)=NaN;

try
    a(a<1)=NaN;
    m = b ./ (-1+a);
    a(a<2)=NaN;
    v = (b.^2) ./ ((-2+a).*((-1+a).^2));
catch
    error('invgamstat:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


end