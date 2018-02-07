function [logL, avar] = pascallike(params, data)
%PASCALLIKE Negative of the Pascal Distribution log-likelihood function.
%   L = PASCALLIKE(PARAMS,DATA) returns the negative of the Pascal
%   Distribution log-likelihood function for the parameters PARAMS(1) = N
%   and PARAMS(2) = P, given DATA.
%
%   [LOGL, AVAR] = PASCALLIKE(PARAMS,DATA) adds the inverse of Fisher's
%   information matrix, AVAR. If the input parameter values in PARAMS
%   are the maximum likelihood estimates, the diagonal elements of AVAR
%   are their asymptotic variances.  AVAR is based on the observed
%   Fisher's information, not the expected information.
%
%   PASCALLIKE is a utility function for maximum likelihood estimation. 

%   Mike Sheppard
%   Last Modified 18-Jun-2011

if nargin < 2, 
    error('pascallike:TooFewInputs',...
          'Requires two input argument.'); 
end

[n, m] = size(data);
if min(n,m) > 1
    error('pascallike:InvalidData',...
          'Requires data to be more than one data point');
end

if nargout == 2 & max(m,n) == 1
    error('pascallike:NotEnoughData',...
        'Requires data to be more than one data point');
end

if n == 1
   data = data';
   n = m;
end

r = params(1);
p = params(2);

%Transform to Negative Binomial Distribution
%PASCALPDF(X,R,P)=NBINPDF(X-R,R,P)
data=data-r;
[logL, avar] = nbinlike(params, data);

end
