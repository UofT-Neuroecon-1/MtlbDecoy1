function [m,v] = logistat(mu,sigma)
%LOGISTAT Mean and variance for the Logistic Distribution
%   [M,V] = LOGISTAT(MU,SIGMA) returns the mean and variance for the
%   Logistic Distribution with mean MU and scale parameter SIGMA.
%
%   The sizes of M and V are the common size of the input arguments.
%   A scalar input functions as a constant matrix of the same size as
%   the other inputs.
%
%   Default values for MU and SIGMA are 0 and 1, respectively.
%
%   Type: Continuous, unbounded
%   Restrictions:
%     SIGMA>0
%
%   See also LOGIPDF, LOGICDF, LOGIINV, LOGIFIT, LOGILIKE, LOGIRND
%

%   Mike Sheppard
%   Last Modified 12-Dec-2011


if nargin < 1, mu = 0; end
if nargin < 2, sigma = 1; end


% Return NaN for out of range parameters.
sigma(sigma <= 0) = NaN;

try
    %From addlogi.m in statistics toolbox
    %-----
    m = mu + zeros(size(sigma));
    v = (sigma.^2 .* pi.^2 ./ 3) + zeros(size(mu));
    %-----
catch
    error('logistat:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


end