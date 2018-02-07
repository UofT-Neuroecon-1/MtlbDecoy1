function [m,v] = ricestat(s,sigma)
%RICESTAT Mean and variance for the Rician Distribution.
%   [M,V] = RICESTAT(S,SIGMA) returns the mean and variance of the
%   Rician Distribution with noncentrality parameter S and scale
%   parameter SIGMA
%
%   The sizes of M and V are the common size of the input arguments.
%   A scalar input functions as a constant matrix of the same size as
%   the other inputs.
%
%   Type: Continuous, semi-bounded
%   Restrictions:
%     S >= 0
%     SIGMA > 0
%
%   See also RICEPDF, RICECDF, RICEINV, RICEFIT, RICELIKE, RICERND
%

%   Mike Sheppard
%   Last Modified 10-Dec-2011



if nargin < 2
    error('ricestat:TooFewInputs',...
        'Requires two input arguments.');
end

% Return NaN for out of range parameters.
s(s < 0) = NaN;
sigma(sigma <= 0) = NaN;

try
    %From addrice.m in statistics toolbox
    %-----
    t = .5 .* (s./sigma).^2;
    m = sigma.*sqrt(.5.*pi) .* ((1+t).*besseli(0,.5.*t,1) + t.*besseli(1,.5.*t,1));
    v = 2.*sigma.^2 + s.^2 - m.^2;
    %-----
catch
    error('ricestat:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


end