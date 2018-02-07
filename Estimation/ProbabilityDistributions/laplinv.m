function x = laplinv(p,mu,sigma)
%LAPLINV Inverse of the Laplace cumulative distribution function
%   X = LOGLINV(P,MU,SIGMA) returns the inverse cumulative distribution 
%   function of the Laplace distribution with location MU and scale SIGMA, 
%   evaluated at the values in P.
%
%   The size of X is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Default values for MU and SIGMA are 0 and 1, respectively.
%
%   Type: Continuous, unbounded
%   Restrictions:
%     SIGMA>0
%
%   See also LAPLPDF, LAPLCDF, LAPLSTAT, LAPLFIT, LAPLLIKE, LAPLRND
%

%   Mike Sheppard
%   Last Modified 10-Dec-2011



if nargin < 1
    error('laplinv:TooFewInputs',...
          'Requires at least one input argument.'); 
end
if nargin < 2, mu = 0; end
if nargin < 3, sigma = 1; end


[errorcode, p,mu,sigma] = distchck(3,p,mu,sigma);

if errorcode > 0
    error('laplinv:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize y to zero.
if isa(p,'single') || isa(mu,'single') || isa(sigma,'single')
   x = zeros(size(p),'single');
else
   x = zeros(size(p));
end

% Return NaN for out of range parameters.
sigma(sigma <= 0) = NaN;


%Modified from dftoolinittemplate.m in statistics toolbox
%-----
z = zeros(size(p));
t = (p<=.5);
z(t) = log(2*p(t));
z(~t) = -log(1-2*(p(~t)-.5));
x = mu + sigma.*z;
%-----



end