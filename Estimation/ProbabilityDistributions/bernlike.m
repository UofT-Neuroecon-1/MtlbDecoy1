function [nlogL,avar] = bernlike(param,data)
%BERNLIKE Negative log-likelihood for the Bernoulli distribution.
%   NLOGL = BERNLIKE(PARAM,DATA) returns the negative of the log-likelihood
%   for the Bernoulli distribution, evaluated at the success parameter
%   PARAM, given DATA.  NLOGL is a scalar.
%
%   [NLOGL, AVAR] = BERNLIKE(PARAM,DATA) returns the inverse of Fisher's
%   information, AVAR, a scalar.  If the input parameter value in PARAM is
%   the maximum likelihood estimate, AVAR is its asymptotic variance.  AVAR
%   is based on the observed Fisher's information, not the expected
%   information.
%
%   Distribution: Discrete, bounded, {0,1}
%   Restrictions:
%        0 <= P <= 1
%
%   See also BERNPDF, BERNCDF, BERNINV, BERNSTAT,
%            BERNFIT, BERNRND, BERNSF, BERNHAZ
%

%   Mike Sheppard
%   Last Modified 15-May-2012

if nargin < 2
    error('bernlike:TooFewInputs',...
        'Requires two input arguments.');
elseif ~isvector(data)
    error('bernlike:VectorRequired',...
        'Requires DATA to be a vector');
end
x = data(:);

p = param;

% Return NaN for out of range parameters or data
p(p<0 | p>1)=NaN;
x(~(x==0 | x==1)) = NaN;

% Sum up the individual log-likelihood terms, and return the negative
% log-likelihood.
sx=sum(x); n=numel(x);
nlogL = - (sx.*log(p) + (n-sx).*log(1-p));

% Compute the negative hessian at the parameter value, and invert to get
% the observed information.
if nargout == 2
    nhess = (sx./p.^2) + ((n-sx)./(1-p).^2);
    avar = (1./nhess);
end

end
