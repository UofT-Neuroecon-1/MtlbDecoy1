function [m,v] = laplstat(mu,sigma)
%LAPLSTAT Mean and variance for the Laplace Distribution
%   [M,V] = LAPLSTAT(MU,SIGMA) returns the mean and variance for the
%   Laplace Distribution with location MU and scale SIGMA.
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
%   See also LAPLPDF, LAPLCDF, LAPLINV, LAPLFIT, LAPLLIKE, LAPLRND
%

%   Mike Sheppard
%   Last Modified 12-Dec-2011



if nargin < 1, mu = 0; end
if nargin < 2, sigma = 1; end


% Return NaN for out of range parameters.
sigma(sigma <= 0) = NaN;

try
    %Modified from dftoolinittemplate.m in statistics toolbox
    %-----
    m = mu + zeros(size(sigma));
    v = 2 * sigma.^2 + zeros(size(mu));
    %-----
catch
    error('laplstat:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


end