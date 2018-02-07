function [s,slo,sup] = expsf(x,mu,pcov,alpha)
%EXPSF Exponential survival function.
%   S = EXPSF(X,MU) returns the survival function of the exponential
%   distribution with mean parameter MU, evaluated at the values in X.
%
%   The size of S is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   The default value for MU is 1.
%
%   [S,SLO,SUP] = EXPSF(X,MU,PCOV,ALPHA) produces confidence bounds
%   for S when the input parameter MU is an estimate.  PCOV is the
%   variance of the estimated MU.  ALPHA has a default value of 0.05, and
%   specifies 100*(1-ALPHA)% confidence bounds.  SLO and SUP are arrays of
%   the same size as S containing the lower and upper confidence bounds.
%   The bounds are based on a normal approximation for the distribution of
%   the log of the estimate of MU.  You can get a more accurate set of
%   bounds simply by using EXPFIT to get a confidence interval for MU,
%   and evaluating EXPSF at the lower and upper end points of that interval.
%
%   See also EXPPDF, EXPCDF, EXPINV, EXPSTAT, EXPFIT, EXPLIKE,
%            EXPRND, EXPHAZ
%

%   Mike Sheppard
%   Last Modified 26-Dec-2011


if nargin < 1
    error('expsf:TooFewInputs',...
        'Requires at least one input argument.');
end
if nargin < 2
    mu = 1;
end

        
% More checking if we need to compute confidence bounds.
if nargout>1
    if nargin<3
        error('expsf:TooFewInputs',...
            'Must provide parameter variance to compute confidence bounds.');
    end
    if numel(pcov)~=1
        error('expsf:BadVariance','Variance must be a scalar.');
    end
    if nargin<4
        alpha = 0.05;
    elseif ~isnumeric(alpha) || numel(alpha)~=1 || alpha<=0 || alpha>=1
        error(message('expsf:BadAlpha'));
    end
    if (nargout>=2)&&(pcov<0)
      error('expsf:BadVariance','PCOV must be non-negative.');
    end
end


try
    if nargout==1
        s = 1 - expcdf(x,mu);
    else
        [p,plo,pup] = expcdf(x,mu,pcov,alpha);
        [s,sup,slo] = deal(1-p,1-plo,1-pup);    %switch low/high
    end
catch
    error('expsf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end
