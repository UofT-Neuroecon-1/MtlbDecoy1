function [m,v] = rcosstat(u,s)
%RCOSSTAT Mean and variance for the raised cosine distribution
%   [M,V] = RCOSSTAT(U,S) returns the mean and variance of the raised
%   cosine distribution with mean U and scale parameter S.
%
%   Type: Continuous, Bounded
%   Restrictions:
%   S>0
%   U-S<=X<=U+S
%
%   The size of P is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%

%   Mike Sheppard
%   Last Modified 13-Jun-2011


if nargin < 2
    error('rcosstat:TooFewInputs',...
        'Requires two input arguments.');
end


% Return NaN for out of range parameters.
s(s<=0)=NaN;

try
    m=u+zeros(size(s));
    const=(1/3)-(2/(pi^2));
    v=(const.*(s.^2))+zeros(size(u));
catch
    error('rcosstat:InputSizeMismatch',...
        'Non-scalar arguments must match in size.');
end



end