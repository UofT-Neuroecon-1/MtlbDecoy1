function h = lognhaz(x,mu,sigma)
%LOGNHAZ Lognormal hazard function
%   H = LOGNHAZ(X,MU,SIGMA) returns the hazard function of the lognormal
%   distribution with parameters MU and SIGMA, at the values in X.
%   MU and SIGMA are the mean and standard deviation, respectively, of the
%   associated normal distribution.
%
%   The size of H is the common size of X, MU and SIGMA.  A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Default values for MU and SIGMA are 0 and 1, respectively.
%
%   See also LOGNPDF, LOGNCDF, LOGNINV, LOGNSTAT, LOGNFIT, LOGNLIKE,
%            LOGNRND, LOGNSF
%

%   Mike Sheppard
%   Last Modified 26-Dec-2011


if nargin<1
    error('lognhaz:TooFewInputs','Input argument X is undefined.');
end
if nargin < 2
    mu = 0;
end
if nargin < 3
    sigma = 1;
end


try
    h = lognpdf(x,mu,sigma) ./ lognsf(x,mu,sigma);
catch
    error('lognhaz:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end