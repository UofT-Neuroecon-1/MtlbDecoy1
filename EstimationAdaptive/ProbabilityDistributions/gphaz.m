function h = gphaz(x,k,sigma,theta)
%GPHAZ Generalized Pareto hazard function
%   H = GPHAZ(X,K,SIGMA,THETA) returns the hazard function of the 
%   generalized Pareto (GP) distribution with tail index (shape) 
%   parameter K, scale parameter SIGMA, and threshold (location) 
%   parameter THETA, evaluated at the values in X.
%
%   The size of H is the common size of the input arguments.  A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Default values for K, SIGMA, and THETA are 0, 1, and 0, respectively.
%
%   When K = 0 and THETA = 0, the GP is equivalent to the exponential
%   distribution.  When K > 0 and THETA = SIGMA/K, the GP is equivalent to the
%   Pareto distribution.  The mean of the GP is not finite when K >= 1, and the
%   variance is not finite when K >= 1/2.  When K >= 0, the GP has positive
%   density for X>THETA, or, when K < 0, for 0 <= (X-THETA)/SIGMA <= -1/K.
%
%   See also GPPDF, GPCDF, GPINV, GPSTAT, GPFIT, GPLIKE, GPRND, GPSF
%

%   Mike Sheppard
%   Last Modified 26-Dec-2011

if nargin < 1
    error(message('gphaz:TooFewInputs'));
end
if nargin < 2 || isempty(k), k = 0;     end
if nargin < 3 || isempty(sigma), sigma = 1; end
if nargin < 4 || isempty(theta), theta = 0; end


try
    h = gppdf(x,k,sigma,theta) ./ gpsf(x,k,sigma,theta);
catch
    error('gphaz:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end