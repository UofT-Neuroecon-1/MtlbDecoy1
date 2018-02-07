function p = laplcdf(x,mu,sigma)
%LAPLPDF Laplace cumulative distribution function
%   P = LOGLCDF(X,MU,SIGMA) returns the cumulative distribution function
%   of the Laplace Distribution with location MU and scale SIGMA, evaluated
%   at the values in X.
%
%   The size of P is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Default values for MU and SIGMA are 0 and 1, respectively.
%
%   Distribution: Continuous, unbounded, (-Inf, Inf)
%   Restrictions:
%      SIGMA > 0
%
%   See also LAPLPDF, LAPLINV, LAPLSTAT, LAPLFIT,
%            LAPLLIKE, LAPLRND, LAPLSF, LAPLHAZ
%


%   Mike Sheppard
%   Last Modified 10-Dec-2011


if nargin < 1
    error('laplcdf:TooFewInputs',...
          'Requires at least one input argument.'); 
end
if nargin < 2, mu = 0; end
if nargin < 3, sigma = 1; end


[errorcode, x,mu,sigma] = distchck(3,x,mu,sigma);

if errorcode > 0
    error('laplcdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize p to zero.
if isa(x,'single') || isa(mu,'single') || isa(sigma,'single')
   p = zeros(size(x),'single');
else
   p = zeros(size(x));
end

% Return NaN for out of range parameters.
sigma(sigma <= 0) = NaN;

%Modified from dftoolinittemplate.m in statistics toolbox
%-----
z = (x-mu)./sigma;
t = (z<=0);
p(t) = exp(z(t))/2;
p(~t) = 1 - exp(-z(~t))/2;
%-----

%Catch round off
p(p<0)=0; p(p>1)=1;

end


