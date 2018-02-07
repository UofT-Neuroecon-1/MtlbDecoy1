function nlogL = asinlike(params,data)
%ASINLIKE Negative arcsine log-likelihood function.
%   NLOGL = ASINLIKE(PARAMS,DATA) returns the negative of the log-likelihood
%   for the Arcsine Distribution, evaluated at parameters PARAMS(1) = A and
%   PARAMS(2) = B, given DATA.  NLOGL is a scalar.
%
%   Distribution: Continuous, bounded, [A,B]
%   Restrictions:
%         A < B
%
%   See also ASINPDF, ASINCDF, ASININV, ASINSTAT, 
%            ASINFIT, ASINRND, ASINSF, ASINHAZ
%

%   Mike Sheppard
%   Last Modified 25-May-2012


if nargin ~=2
    error('asinlike:TooFewInputs','Requires two input arguments.'); 
end

data_n = (data-params(1))./(params(2)-params(1));
%Special case BETALIKE
nlogL = betalike([1 1]/2,data_n);

%d2= 1/ (2*(a-x)^2)
%d2 = 1/ (2*(b-x)^2)


end
