function [s,slo,sup] = gamsf(x,a,b,pcov,alpha)
%GAMSF Gamma survival function.
%   S = GAMSF(X,A,B) returns the survival function of the gamma distribution
%   with shape and scale parameters A and B, respectively, at the values in
%   X.  The size of S is the common size of the input arguments.  A scalar
%   input functions as a constant matrix of the same size as the other
%   inputs.
%
%   Some references refer to the gamma distribution with a single
%   parameter.  This corresponds to the default of B = 1.
%
%   [S,SLO,SUP] = GAMSF(X,A,B,PCOV,ALPHA) produces confidence bounds for
%   S when the input parameters A and B are estimates.  PCOV is a 2-by-2
%   matrix containing the covariance matrix of the estimated parameters.
%   ALPHA has a default value of 0.05, and specifies 100*(1-ALPHA)%
%   confidence bounds.  SLO and SUP are arrays of the same size as S
%   containing the lower and upper confidence bounds.
%
%   See also GAMPDF, GAMCDF, GAMINV, GAMSTAT, GAMFIT, GAMLIKE, 
%            GAMRND, GAMHAZ
%

%   Mike Sheppard
%   Last Modified 26-Dec-2011

if nargin < 2
    error('gamsf:TooFewInputs',...
        'Requires at least two input arguments.');
elseif nargin < 3
    b = 1;
end

try
    s = 1 - gamcdf(x,a,b);
catch
    error('gamsf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

if nargout>=2
    switch nargin
        case 3
            [p,plo,pup] = gamcdf(x,a,b);
        case 4
            [p,plo,pup] = gamcdf(x,a,b,pcov);
        case 5
            [p,plo,pup] = gamcdf(x,a,b,pcov,alpha);
    end
   [s,sup,slo] = deal(1-p,1-plo,1-pup);    %switch low/high
end

end