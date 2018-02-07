function h = invgamhaz(x,a,b)
%INVGAMHAZ Inverse-Gamma hazard function
%   H = INVGAMHAZ(X,A,B) returns the hazard function of the Inverse-Gamma
%   Distribution with shape parameter A and scale parameter B, evaluated
%   at the values in X.
%
%   The size of H is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Type: Continuous, semi-bounded, (0,Inf)
%   Restrictions:
%     A, B > 0
%
%   See also INVGAMPDF, INVGAMCDF, INVGAMINV, INVGAMSTAT, 
%            INVGAMFIT, INVGAMLIKE, INVGAMRND, INVGAMSF
%

%   Mike Sheppard
%   Last Modified: 13-May-2012


if nargin < 3
   error('invgamhaz:TooFewInputs','Requires three input arguments.');
end


try
    h = invgampdf(x,a,b) ./ invgamsf(x,a,b);
catch
    error('invgamhaz:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


end