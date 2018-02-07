function h = evhaz(x,mu,sigma)
%EVHAZ Extreme value hazard function
%   H = EVHAZ(X,MU,SIGMA) returns the hazard function of the type 1
%   extreme value distribution with location parameter MU and scale
%   parameter SIGMA, evaluated at the values in X.
%
%   The size of H is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Default values for MU and SIGMA are 0 and 1, respectively.
%
%   See also EVPDF, EVCDF, EVINV, EVSTAT, EVFIT, EVLIKE, EVRND, EVSF
%

%   Mike Sheppard
%   Last Modified 26-Dec-2011

if nargin < 1
    error('evhaz:TooFewInputs',...
        'Requires at least one input argument.');
end
if nargin < 2
    mu = 0;
end
if nargin < 3
    sigma = 1;
end

try
    h = evpdf(x,mu,sigma) ./ evsf(x,mu,sigma);
catch
    error('evhaz:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


end
