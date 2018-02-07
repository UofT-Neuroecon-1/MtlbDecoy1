function y = invgpdf(x,mu,lambda)
%INVGPDF Inverse Gaussian probability density function
%   Y=INVGPDF(X,MU,LAMBDA) returns the probability density function of the
%   Inverse Gaussian Distribution with mean MU and scale parameter LAMBDA,
%   evaluated at the values in X. 
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Distribution: Continuous, semi-bounded, (0,Inf)
%   Restrictions:
%      MU, LAMBDA > 0
%
%   NOTE: The Inverse Gaussian Distribution is also known as the 
%   Wald Distribution.
%
%   See also INVGCDF, INVGINV, INVGSTAT, INVGFIT, 
%            INVGLIKE, INVGRND, INVGSF, INVGHAZ
%

%   Mike Sheppard
%   Last Modified: 7-Dec-2011

if nargin ~= 3
    error('invgpdf:TooFewInputs',...
          'Requires three input arguments.'); 
end

[errorcode x mu lambda] = distchck(3,x,mu,lambda);

if errorcode > 0
    error('invgpdf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

% Initialize Y to zero.
if isa(x,'single') || isa(mu,'single') || isa(lambda,'single')
    y = zeros(size(x),'single');
else
    y = zeros(size(x));
end

% Return NaN for out of range parameters.
mu(mu <= 0) = NaN;
lambda(lambda <= 0) = NaN;


%From addinvg.m in statistics toolbox
%-----
nonpos = (x <= 0);
x(nonpos)= realmin;
y = sqrt(lambda./(2.*pi.*x.^3)) .* exp(-0.5.*lambda.*(x./mu - 2 + mu./x)./mu);
% this would happen automatically for x==0, but generates DivideByZero warnings
y(nonpos) = 0;
%-----

%Catch round off
y(y<0)=0;


end