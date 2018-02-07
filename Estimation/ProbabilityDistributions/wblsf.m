function [s,slo,sup] = wblsf(x,A,B,pcov,alpha)
%WBLSF Weibull survival function
%   S = WBLSF(X,A,B) returns the survival function of the Weibull distribution
%   with scale parameter A and shape parameter B, evaluated at the
%   values in X.  The size of S is the common size of the input arguments.
%   A scalar input functions as a constant matrix of the same size as the
%   other inputs.
%
%   Default values for A and B are 1 and 1, respectively.
%
%   [S,SLO,SUP] = WBLSF(X,A,B,PCOV,ALPHA) produces confidence
%   bounds for S when the input parameters A and B are estimates.
%   PCOV is a 2-by-2 matrix containing the covariance matrix of the estimated
%   parameters.  ALPHA has a default value of 0.05, and specifies
%   100*(1-ALPHA)% confidence bounds.  SLO and SUP are arrays of the same
%   size as P containing the lower and upper confidence bounds.
%
%   See also WBLPDF, WBLCDF, WBLINV, WBLSTAT, WBLFIT, WBLLIKE,
%            WBLRND, WBLHAZ
%

%   Mike Sheppard
%   Last Modified 26-Dec-2011

if nargin<1
    error('wblsf:TooFewInputs','Input argument X is undefined.');
end
if nargin < 2
    A = 1;
end
if nargin < 3
    B = 1;
end

try
    s = 1 - wblcdf(x,A,B);
catch
    error('wblsf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

if nargout>=2
    switch nargin
        case 3
            [p,plo,pup] = wblcdf(x,A,B);
        case 4
            [p,plo,pup] = wblcdf(x,A,B,pcov);
        case 5
            [p,plo,pup] = wblcdf(x,A,B,pcov,alpha);
    end
    [s,sup,slo] = deal(1-p,1-plo,1-pup);    %switch low/high
end

end