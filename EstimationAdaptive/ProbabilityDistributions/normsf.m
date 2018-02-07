function [s,slo,sup] = normsf(x,mu,sigma,pcov,alpha)
%NORMSF Normal survival function
%   S = NORMSF(X,MU,SIGMA) returns the survival function of the normal 
%   distribution with mean MU and standard deviation SIGMA, evaluated at 
%   the values in X.
%
%   The size of S is the common size of X, MU and SIGMA.  A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Default values for MU and SIGMA are 0 and 1, respectively.
%
%   [S,SLO,SUP] = NORMSF(X,MU,SIGMA,PCOV,ALPHA) produces confidence bounds
%   for S when the input parameters MU and SIGMA are estimates.  PCOV is a
%   2-by-2 matrix containing the covariance matrix of the estimated parameters.
%   ALPHA has a default value of 0.05, and specifies 100*(1-ALPHA)% confidence
%   bounds.  SLO and SUP are arrays of the same size as S containing the lower
%   and upper confidence bounds.
%
%   See also NORMPDF, NORMCDF, NORMINV, NORMSTAT, NORMFIT, NORMLIKE, 
%            NORMRND, NORMHAZ
%

%   Mike Sheppard
%   Last Modified 26-Dec-2011


if nargin<1
    error('normsf:TooFewInputs','Input argument X is undefined.');
end
if nargin < 2
    mu = 0;
end
if nargin < 3
    sigma = 1;
end


try
    s = 1 - normcdf(x,mu,sigma);
catch
    error('normsf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

if nargout>=2
    switch nargin
        case 3
            [p,plo,pup] = normcdf(x,mu,sigma);
        case 4
            [p,plo,pup] = normcdf(x,mu,sigma,pcov);
        case 5
            [p,plo,pup] = normcdf(x,mu,sigma,pcov,alpha);
    end
   [s,sup,slo] = deal(1-p,1-plo,1-pup);    %switch low/high
end

end
