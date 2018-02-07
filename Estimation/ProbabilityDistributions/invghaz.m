function h = invghaz(x,mu,lambda)
%INVGHAZ Inverse Gaussian hazard function
%   H=INVGHAZ(X,MU,LAMBDA) returns the hazard function
%   of the Inverse Gaussian Distribution with mean MU and scale parameter
%   LAMBDA, evaluated at the values in X. 
%
%   The size of H is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Type: Continuous, semi-bounded, (0,Inf)
%   Restrictions:
%        MU, LAMBDA > 0
%
%   Note: The Inverse Gaussian Distribution is also known as the 
%   Wald Distribution.
%
%   See also INVGPDF, INVGCDF, INVGINV, INVGSTAT, INVGFIT, INVGLIKE, 
%            INVGRND, INVGSF
%

%   Mike Sheppard
%   Last Modified: 3-Jan-2012


if nargin < 3
    error('invghaz:TooFewInputs',...
          'Requires three input arguments.'); 
end


try
    h = invgpdf(x,mu,lambda) ./ invgsf(x,mu,lambda);
catch
    error('invghaz:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end