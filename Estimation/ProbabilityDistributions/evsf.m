function [s,slo,sup] = evsf(x,mu,sigma,pcov,alpha)
%EVSF Extreme value survival function
%   S = EVSF(X,MU,SIGMA) returns the survival function of the type 1
%   extreme value distribution with location parameter MU and scale
%   parameter SIGMA, evaluated at the values in X.
%
%   The size of S is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Default values for MU and SIGMA are 0 and 1, respectively.
%
%   [S,SLO,SUP] = EVSF(X,MU,SIGMA,PCOV,ALPHA) produces confidence bounds
%   for S when the input parameters MU and SIGMA are estimates.  PCOV is a
%   2-by-2 matrix containing the covariance matrix of the estimated parameters.
%   ALPHA has a default value of 0.05, and specifies 100*(1-ALPHA)% confidence
%   bounds.  SLO and SUP are arrays of the same size as S containing the lower
%   and upper confidence bounds.
%
%   The type 1 extreme value distribution is also known as the Gumbel
%   distribution.  The version used here is suitable for modeling minima; the
%   mirror image of this distribution can be used to model maxima by negating
%   X.  If Y has a Weibull distribution, then X=log(Y) has the type 1 extreme
%   value distribution.
%
%   See also EVPDF, EVCDF, EVINV, EVSTAT, EVFIT, EVLIKE, EVRND, EVHAZ
%

%   Mike Sheppard
%   Last Modified 26-Dec-2011

if nargin < 1
    error('evsf:TooFewInputs',...
        'Requires at least one input argument.');
end
if nargin < 2
    mu = 0;
end
if nargin < 3
    sigma = 1;
end

% More checking if we need to compute confidence bounds.
if nargout>1
    if nargin<4
        error('evsf:TooFewInputs',...
            'Must provide covariance matrix to compute confidence bounds.');
    end
    if ~isequal(size(pcov),[2 2])
        error('evsf:BadCovariance',...
            'Covariance matrix must have 2 rows and columns.');
    end
    if nargin<5
        alpha = 0.05;
    elseif ~isnumeric(alpha) || numel(alpha)~=1 || alpha<=0 || alpha>=1
        error(message('evsf:BadAlpha'));
    end
    if nargout>=2
        try
            z = (x-mu)./sigma;
        catch
            error(message('evsf:InputSizeMismatch'));
        end
        zvar = (pcov(1,1) + 2*pcov(1,2)*z + pcov(2,2)*z.^2) ./ (sigma.^2);
        if any(zvar<0)
            error('evsf:BadCovariance',...
                'PCOV must be a positive semi-definite matrix.');
        end
    end
end


try
    if nargout==1
        s = 1 - evcdf(x,mu,sigma);
    else
        [p,plo,pup] = evcdf(x,mu,sigma,pcov,alpha);
        [s,sup,slo] = deal(1-p,1-plo,1-pup);    %switch low/high
    end
catch
    error('evsf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end
