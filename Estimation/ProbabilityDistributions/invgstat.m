function [m,v] = invgstat(mu,lambda)
%INVGSTAT Mean and variance for the Inverse Gaussian Distribution
%   [M,V] = INVGSTAT(MU,LAMBDA) returns the mean and variance of the
%   Inverse Gaussian Distribution with mean MU and scale parameter
%   LAMBDA.
%
%   The sizes of M and V are the common size of the input arguments.  
%   A scalar input functions as a constant matrix of the same size as
%   the other inputs.
%
%   Type: Continuous, semi-bounded, (0,Inf)
%   Restrictions:
%        MU, LAMBDA > 0
%
%   Note: The Inverse Gaussian Distribution is also known as the
%   Wald Distribution.
%
%   See also INVGPDF, INVGCDF, INVGINV, INVGFIT, INVGLIKE, INVGRND, 
%            INVGSF, INVGHAZ
%

%   Mike Sheppard
%   Last Modified: 7-Dec-2011

if nargin < 2
    error('invgstat:TooFewInputs',...
          'Requires two input arguments.'); 
end


% Return NaN for out of range parameters.
mu(mu <= 0) = NaN;
lambda(lambda <= 0) = NaN;

try
    m = mu;
    v = mu.^3 ./ lambda;
catch
    error('invgstat:InputSizeMismatch',...
          'Non-scalar arguments must match in size.');
end

end