function h = expgamhaz(x,k,t,mu)
%EXPGAMHAZ Exp-Gamma hazard function
%   H = EXPGAMHAZ(X,K,T,MU) returns the hazard function
%   for the Exp-Gamma Distribution with shape parameter K, scale 
%   parameter T, and location parameter MU, at the values of X.
%
%   The size of S is the common sizes of the inputs. A scalar input   
%   functions as a constant matrix of the same size as the other input.    
%
%   Type: Continuous, unbounded, (-Inf,Inf)
%   Restrictions:
%        K, T > 0
%
%   Note: The Exp-Gamma Distribution is also known as the Generalized
%   Extreme Value Distribution, and also should not be confused with the
%   Log-Gamma Distribution.
%
%   See also EXPGAMPDF, EXPGAMCDF, EXPGAMINV, EXPGAMSTAT, 
%            EXPGAMFIT, EXPGAMLIKE, EXPGAMRND, EXPGAMSF
%

%   Mike Sheppard
%   Last Modified 20-Dec-2011

if nargin < 4
    error('expgamhaz:TooFewInputs',...
          'Requires at least four input argument.'); 
end


try
    h = expgampdf(x,k,t,mu) ./ expgamsf(x,k,t,mu);
catch
    error('expgamhaz:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end