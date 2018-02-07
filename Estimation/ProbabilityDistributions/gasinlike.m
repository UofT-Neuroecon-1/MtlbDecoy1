function logL = gasinlike(params,data)
%GASINLIKE Negative generalized arcsine log-likelihood function.
%   LOGL = GASINLIKE(PARAMS,DATA) returns the negative of the generalized 
%   arcsine log-likehood function for the parameters PARAMS(1) = ALPHA,
%   PARAMS(2) = A, and PARAMS(3)=B, given DATA.
%
%   Type: Continuous, bounded
%
%   See also GASINPDF, GASINCDF, GASININV, GASINSTAT, GASINFIT, GASINRND
%

%   Mike Sheppard
%   Last Modified 26-Apr-2011

if nargin < 2
    error('gasinlike:TooFewInputs','Requires at least two input arguments.'); 
end

if isscalar(params)
    params=[params 0 1];
end

data_normalized = (data-params(2))./(params(3)-params(2));
%Special case BETALIKE
logL = betalike([1-params(1) params(1)],data_normalized);


end
