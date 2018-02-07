function [m,v] = loglstat(mu,sigma)
%LOGLSTAT Mean and variance for the Log-Logistic Distribution
%   [M,V] = LOGLSTAT(MU,SIGMA) returns the mean and variance for the
%   Log-Logistic Distribution with log-mean MU and log-scale SIGMA.
%
%   The sizes of M and V are the common size of the input arguments.
%   A scalar input functions as a constant matrix of the same size as
%   the other inputs.
%
%   Default values for MU and SIGMA are 0 and 1, respectively.
%
%   Type: Continuous, semi-bounded
%   Restrictions:
%     SIGMA>0
%
%   See also LOGLPDF, LOGLCDF, LOGLINV, LOGLFIT, LOGLLIKE, LOGLRND
%

%   Mike Sheppard
%   Last Modified 12-Dec-2011


if nargin < 1, mu = 0; end
if nargin < 2, sigma = 1; end

try
    %Expand size if necessary
    mu=mu+zeros(size(sigma));
    sigma=sigma+zeros(size(mu));
catch
    error('loglstat:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

%pre-allocate memory
m=zeros(size(mu)); v=m;

% Return NaN for out of range parameters.
sigma(sigma <= 0) = NaN;


%Modified from addlogi.m in statistics toolbox
%-----

k=(sigma < 1);
m(k) = exp(mu(k) + gammaln(1+sigma(k)) + gammaln(1-sigma(k)));
m(~k) = Inf;

k=(sigma < .5);
v(k) = exp(2.*mu(k) + gammaln(1+2.*sigma(k)) + gammaln(1-2.*sigma(k))) - m(k).^2;
v(~k) = Inf;

%-----


end