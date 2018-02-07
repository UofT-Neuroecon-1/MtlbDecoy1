function p = invgcdf(x,mu,lambda)
%INVGCDF Inverse Gaussian cumulative distribution function
%   P=INVGCDF(X,MU,LAMBDA) returns the cumulative distribution function
%   of the Inverse Gaussian Distribution with mean MU and scale parameter
%   LAMBDA, evaluated at the values in X. 
%
%   The size of P is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Distribution: Continuous, semi-bounded, (0,Inf)
%   Restrictions:
%      MU, LAMBDA > 0
%
%   NOTE: The Inverse Gaussian Distribution is also known as the 
%   Wald Distribution.
%
%   See also INVGPDF, INVGINV, INVGSTAT, INVGFIT, 
%            INVGLIKE, INVGRND, INVGSF, INVGHAZ
%

%   Mike Sheppard
%   Last Modified: 7-Dec-2011

if nargin ~= 3
    error('invgcdf:TooFewInputs',...
          'Requires three input arguments.'); 
end

[errorcode x mu lambda] = distchck(3,x,mu,lambda);

if errorcode > 0
    error('invgcdf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

% Initialize Y to zero.
if isa(x,'single') || isa(mu,'single') || isa(lambda,'single')
    p = zeros(size(x),'single');
else
    p = zeros(size(x));
end

% Return NaN for out of range parameters.
mu(mu <= 0) = NaN;
lambda(lambda <= 0) = NaN;


%From addinvg.m in statistics toolbox
%-----
nonpos = (x <= 0);
posinf = (x == Inf);
x(nonpos)= realmin;
z1 = (x./mu - 1).*sqrt(lambda./x);
z2 = -(x./mu + 1).*sqrt(lambda./x);
p = 0.5.*erfc(-z1./sqrt(2)) + exp(2.*lambda./mu) .* 0.5.*erfc(-z2./sqrt(2));
% this would happen automatically for x==0, but generates DivideByZero warnings
p(nonpos) = 0;
p(posinf) = 1;
%-----

%Catch round off
p(p<0)=0; p(p>1)=1;

end