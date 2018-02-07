function h = gevhaz(x,k,sigma,mu)
%GEVHAZ Generalized extreme value hazard function
%   H = GEVHAZ(X,K,SIGMA,MU) returns the hazard function of the 
%   generalized extreme value (GEV) distribution with shape parameter K, 
%   scale parameter SIGMA, and location parameter MU, evaluated at the
%   values in X.  The size of S is the common size of the input arguments. 
%   A scalar input functions as a constant matrix of the same size as the
%   other inputs.
%
%   Default values for K, SIGMA, and MU are 0, 1, and 0, respectively.
%
%   When K < 0, the GEV is the type III extreme value distribution.  When K >
%   0, the GEV distribution is the type II, or Frechet, extreme value
%   distribution.  If W has a Weibull distribution as computed by the WBLCDF
%   function, then -W has a type III extreme value distribution and 1/W has a
%   type II extreme value distribution.  In the limit as K approaches 0, the
%   GEV is the mirror image of the type I extreme value distribution as
%   computed by the EVCDF function.
%
%   The mean of the GEV distribution is not finite when K >= 1, and the
%   variance is not finite when K >= 1/2.  The GEV distribution has positive
%   density only for values of X such that K*(X-MU)/SIGMA > -1.
%
%   See also GEVPDF, GEVCDF, GEVINV, GEVSTAT, GEVFIT, GEVLIKE, 
%            GEVRND, GEVSF
%

%   Mike Sheppard
%   Last Modified 26-Dec-2011

if nargin < 1
    error(message('gevhaz:TooFewInputs'));
end
if nargin < 2, k = 0;     end
if nargin < 3, sigma = 1; end
if nargin < 4, mu = 0;    end


try
    h = gevpdf(x,k,sigma,mu) ./ gevsf(x,k,sigma,mu);
catch
    error('gevhaz:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end

