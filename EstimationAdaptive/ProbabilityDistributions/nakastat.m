function [m,v] = nakastat(mu,omega)
%NAKASTAT Mean and variance of the Nakagami Distribution
%   [M,V] = NAKASTAT(MU,OMEGA) returns the mean and variance of the
%   Nakagami Distribution with shape parameter MU and scale parameter
%   OMEGA.
%
%   The sizes of M and V are the common size of the input arguments.
%   A scalar input functions as a constant matrix of the same size as
%   the other inputs.
%
%   Type: Continuous, semi-bounded
%   Restrictions:
%     MU, OMEGA > 0
%
%   See also NAKAPDF, NAKACDF, NAKAINV, NAKAFIT, NAKALIKE, NAKARND
%

%   Mike Sheppard
%   Last Modified 12-Dec-2011


if nargin < 2
    error('nakastat:TooFewInputs',...
        'Requires two input arguments.');
end


% Return NaN for out of range parameters.
mu(mu <= 0) = NaN;
omega(omega <= 0) = NaN;

try
    %From addnaka.m in statistics toolbox
    %-----
    gamratio = exp(gammaln(mu+.5) - gammaln(mu));
    m = gamratio .* sqrt(omega./mu);
    v = omega .* (1 - gamratio.^2 ./ mu);
    %-----
catch
    error('nakastat:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


end