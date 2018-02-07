function s = invgsf(x,mu,lambda)
%INVGSF Inverse Gaussian survival function
%   S=INVGSF(X,MU,LAMBDA) returns the survival function
%   of the Inverse Gaussian Distribution with mean MU and scale parameter
%   LAMBDA, evaluated at the values in X. 
%
%   The size of S is the common size of the input arguments. A scalar input
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
%            INVGRND, INVGHAZ
%

%   Mike Sheppard
%   Last Modified: 3-Jan-2012


if nargin < 3
    error('invgsf:TooFewInputs',...
          'Requires three input arguments.'); 
end


try
    s = 1 - invgcdf(x,mu,lambda);
catch
    error('invgsf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end